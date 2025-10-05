#!/usr/bin/env python3
import os
import re
import subprocess
import numpy as np
import threading
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed

# -----------------------------
# 입력 디렉토리
# -----------------------------
input_dirs = {
    "ITS": "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_its_trimmed",
    "16S": "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_trimmed"
}

# -----------------------------
# 출력/리포트/통계 경로
# -----------------------------
base_tr_dir = "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastp_qc_tr"
report_dir = os.path.join(base_tr_dir, "reports")
os.makedirs(report_dir, exist_ok=True)

stats_qs_file = "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_ad_qc_tr_final.csv"
os.makedirs(os.path.dirname(stats_qs_file), exist_ok=True)

dada2_dirs = {
    "ITS": "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_its_trimmed_dada2",
    "16S": "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_trimmed_dada2"
}

# -----------------------------
# CSV 실시간 기록 준비
# -----------------------------
csv_lock = threading.Lock()
csv_header = ["File_name", "BioProject", "Seq", "Layout",
              "mean_qs20_bf", "mean_qs25_bf", "mean_qs30_bf",
              "mean_qs20_af", "mean_qs25_af", "mean_qs30_af"]

def append_qs_result_realtime(file_path, record):
    with csv_lock:
        file_exists = os.path.isfile(file_path)
        with open(file_path, "a", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            if not file_exists:
                writer.writeheader()
            writer.writerow(record)

# -----------------------------
# 유틸: JS 배열 파싱
# -----------------------------
_js_number_re = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")

def parse_js_array(array_text):
    nums = _js_number_re.findall(array_text)
    return [float(n) for n in nums]

# -----------------------------
# HTML에서 마지막 mean trace 추출
# -----------------------------
def extract_mean_from_html(html_path, plot_id):
    try:
        with open(html_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except Exception as e:
        print(f"⚠️ HTML 파일 열기 실패: {html_path} — {e}")
        return None, None

    marker = f"Plotly.newPlot('{plot_id}'"
    pos = content.find(marker)
    if pos == -1:
        return None, None

    data_matches = list(re.finditer(r"var\s+data\s*=\s*(\[[\s\S]*?\]);", content))
    if not data_matches:
        return None, None

    chosen_block = None
    for start, end, m in [(m.start(), m.end(), m) for m in data_matches]:
        if end < pos:
            chosen_block = m.group(1)
    if not chosen_block:
        chosen_block = data_matches[-1].group(1)

    mean_matches = list(re.finditer(r"{[^}]*name\s*:\s*['\"]mean['\"][^}]*}", chosen_block))
    if not mean_matches:
        return None, None

    mean_block = mean_matches[-1].group(0)

    x_match = re.search(r"x\s*:\s*\[([^\]]+)\]", mean_block)
    y_match = re.search(r"y\s*:\s*\[([^\]]+)\]", mean_block)
    if not x_match or not y_match:
        return None, None

    try:
        xs = parse_js_array(x_match.group(1))
        ys = parse_js_array(y_match.group(1))
    except Exception as e:
        print(f"⚠️ HTML 숫자 변환 실패: {e}")
        return None, None

    return xs, ys

# -----------------------------
# Q score threshold 위치 찾기 (numpy 최적화)
# -----------------------------
def first_position_below_threshold(xs, ys, threshold, min_pos=50):
    if not xs or not ys:
        return None
    xs_arr = np.array(xs)
    ys_arr = np.array(ys)
    mask = (xs_arr >= min_pos) & (ys_arr <= threshold)
    if np.any(mask):
        return int(xs_arr[mask][0])
    return int(xs_arr[-1])

# -----------------------------
# HTML 분석
# -----------------------------
def analyze_html(html_path, read_num):
    xs_bf, ys_bf = extract_mean_from_html(html_path, f"plot_Before_filtering__read{read_num}__quality")
    xs_af, ys_af = extract_mean_from_html(html_path, f"plot_After_filtering__read{read_num}__quality")

    pos20_bf = first_position_below_threshold(xs_bf, ys_bf, 20)
    pos25_bf = first_position_below_threshold(xs_bf, ys_bf, 25)
    pos30_bf = first_position_below_threshold(xs_bf, ys_bf, 30)

    pos20_af = first_position_below_threshold(xs_af, ys_af, 20)
    pos25_af = first_position_below_threshold(xs_af, ys_af, 25)
    pos30_af = first_position_below_threshold(xs_af, ys_af, 30)

    return pos20_bf, pos25_bf, pos30_bf, pos20_af, pos25_af, pos30_af

# -----------------------------
# fastp 실행
# -----------------------------
def run_fastp(r1, r2, out_r1, out_r2, html_report_path):
    cmd = ["fastp",
           "-i", r1,
           "-o", out_r1,
           "-q", "25",
           "-l", "50",
           "-w", "16",
           "--detect_adapter_for_pe",
           "--thread", "16",
           "--html", html_report_path]
    if r2:
        cmd += ["-I", r2, "-O", out_r2]
    subprocess.run(cmd, check=True)

# -----------------------------
# 파일 단위 처리
# -----------------------------
def process_file(seq_type, root, f):
    dada2_base = dada2_dirs[seq_type]
    rel_path = os.path.relpath(root, input_dirs[seq_type])
    out_dir = os.path.join(dada2_base, rel_path)
    os.makedirs(out_dir, exist_ok=True)

    bioproject = os.path.basename(root)

    if f.endswith("_1_trim.fastq"):
        # PE
        r1 = os.path.join(root, f)
        r2 = r1.replace("_1_trim.fastq", "_2_trim.fastq")
        if not os.path.exists(r2):
            return

        out_r1 = os.path.join(out_dir, f)
        out_r2 = os.path.join(out_dir, os.path.basename(r2))
        html_report = os.path.join(report_dir, f.replace("_1_trim.fastq", "_fastp.html"))

        run_fastp(r1, r2, out_r1, out_r2, html_report)

        # read1
        pos20_bf, pos25_bf, pos30_bf, pos20_af, pos25_af, pos30_af = analyze_html(html_report, 1)
        record1 = {
            "File_name": f,
            "BioProject": bioproject,
            "Seq": seq_type,
            "Layout": "PE",
            "mean_qs20_bf": pos20_bf,
            "mean_qs25_bf": pos25_bf,
            "mean_qs30_bf": pos30_bf,
            "mean_qs20_af": pos20_af,
            "mean_qs25_af": pos25_af,
            "mean_qs30_af": pos30_af
        }
        append_qs_result_realtime(stats_qs_file, record1)

        # read2
        pos20_bf, pos25_bf, pos30_bf, pos20_af, pos25_af, pos30_af = analyze_html(html_report, 2)
        record2 = {
            "File_name": os.path.basename(r2),
            "BioProject": bioproject,
            "Seq": seq_type,
            "Layout": "PE",
            "mean_qs20_bf": pos20_bf,
            "mean_qs25_bf": pos25_bf,
            "mean_qs30_bf": pos30_bf,
            "mean_qs20_af": pos20_af,
            "mean_qs25_af": pos25_af,
            "mean_qs30_af": pos30_af
        }
        append_qs_result_realtime(stats_qs_file, record2)

        print(f"[PE] {f}, {os.path.basename(r2)} QS분석 완료")

    elif f.endswith("_s_trim.fastq"):
        # SE
        r1 = os.path.join(root, f)
        out_r1 = os.path.join(out_dir, f)
        html_report = os.path.join(report_dir, f.replace("_s_trim.fastq", "_fastp.html"))

        run_fastp(r1, None, out_r1, None, html_report)

        pos20_bf, pos25_bf, pos30_bf, pos20_af, pos25_af, pos30_af = analyze_html(html_report, 1)
        record = {
            "File_name": f,
            "BioProject": bioproject,
            "Seq": seq_type,
            "Layout": "SE",
            "mean_qs20_bf": pos20_bf,
            "mean_qs25_bf": pos25_bf,
            "mean_qs30_bf": pos30_bf,
            "mean_qs20_af": pos20_af,
            "mean_qs25_af": pos25_af,
            "mean_qs30_af": pos30_af
        }
        append_qs_result_realtime(stats_qs_file, record)

        print(f"[SE] {f} QS분석 완료")

# -----------------------------
# 메인 (병렬 처리 + 실시간 CSV 기록)
# -----------------------------
def main():
    with ThreadPoolExecutor(max_workers=4) as executor:  # 동시에 4개 파일 처리
        futures = []
        for seq_type, in_dir in input_dirs.items():
            for root, dirs, files in os.walk(in_dir):
                for f in files:
                    if f.endswith("_1_trim.fastq") or f.endswith("_s_trim.fastq"):
                        futures.append(executor.submit(process_file, seq_type, root, f))

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"❌ 처리 중 오류 발생: {e}")

    print(f"✅ 모든 결과가 실시간으로 {stats_qs_file} 에 기록되었습니다.")

if __name__ == "__main__":
    main()
