import pandas as pd
import numpy as np
import os
import shutil
from pathlib import Path

# ==============================
# 경로 설정
# ==============================
input_file = "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_ad_qc_tr_final.csv"
amplicon_file = "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_final_read_cutoff.csv"
output_file_step2 = "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_ad_qc_tr_final_meta_2.csv"
output_file_final = "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_ad_qc_tr_final_meta_final.csv"

fastq_dir_16s = "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_trimmed_dada2"
fastq_dir_its = "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_its_trimmed_dada2"
fastq_removed = "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_dada2_removed"

# ==============================
# 1. CSV 읽기
# ==============================
df = pd.read_csv(input_file, sep=",", dtype=str)
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
df.columns = df.columns.str.strip()

# SRR 추출
df['SRR'] = df['File_name'].str.extract(r'(.+?)(_s|_1|_2)')[0]

# ==============================
# 2. F/R 컬럼 초기화
# ==============================
F_cols = ["mean_qs20_bf_F", "mean_qs25_bf_F", "mean_qs30_bf_F",
          "mean_qs20_af_F", "mean_qs25_af_F", "mean_qs30_af_F"]
R_cols = ["mean_qs20_bf_R", "mean_qs25_bf_R", "mean_qs30_bf_R",
          "mean_qs20_af_R", "mean_qs25_af_R", "mean_qs30_af_R"]

for col_f, col_r, src_col in zip(F_cols, R_cols,
                                 ['mean_qs20_bf','mean_qs25_bf','mean_qs30_bf',
                                  'mean_qs20_af','mean_qs25_af','mean_qs30_af']):
    df[col_f] = df.apply(lambda x: x[src_col] if ("_s" in x['File_name'] or "_1" in x['File_name']) else "", axis=1)
    df[col_r] = df.apply(lambda x: x[src_col] if "_2" in x['File_name'] else "", axis=1)

cols_out = ["SRR", "BioProject", "Seq", "Layout"] + F_cols + R_cols
df = df[cols_out]

# ==============================
# 3. PE 병합
# ==============================
def merge_pe(group):
    if len(group) == 2:
        merged = group.iloc[0].copy()
        for col in R_cols:
            merged[col] = group.iloc[1][col]
        return merged
    else:
        return group.iloc[0]

df = df.groupby("SRR", group_keys=False).apply(merge_pe)

# ==============================
# 4. F/R 합산 → mean_qs*_af
# ==============================
sum_cols = ["mean_qs20_af", "mean_qs25_af", "mean_qs30_af"]
F_af_cols = ["mean_qs20_af_F", "mean_qs25_af_F", "mean_qs30_af_F"]
R_af_cols = ["mean_qs20_af_R", "mean_qs25_af_R", "mean_qs30_af_R"]

for f_col, r_col, new_col in zip(F_af_cols, R_af_cols, sum_cols):
    df[f_col] = pd.to_numeric(df[f_col], errors='coerce').fillna(0)
    df[r_col] = pd.to_numeric(df[r_col], errors='coerce').fillna(0)
    df[new_col] = df[f_col] + df[r_col]

df_subset = df[["SRR", "BioProject", "Seq", "Layout"] + sum_cols]

# ==============================
# 5. 그룹 평균 → mean_qs*_af_final
# ==============================
df_grouped = df_subset.groupby(["BioProject", "Seq"], as_index=False)[sum_cols].mean()
df_grouped = df_grouped.apply(lambda x: np.floor(x) if x.name in sum_cols else x)
df_grouped = df_grouped.rename(columns={
    "mean_qs20_af": "mean_qs20_af_final",
    "mean_qs25_af": "mean_qs25_af_final",
    "mean_qs30_af": "mean_qs30_af_final"
})

df_final = df_subset.merge(df_grouped, on=["BioProject", "Seq"], how="left")

# ==============================
# 6. Amplicon length 병합
# ==============================
df_amp = pd.read_csv(amplicon_file, dtype=str)
df_amp = df_amp.applymap(lambda x: x.strip() if isinstance(x, str) else x)
df_amp.columns = df_amp.columns.str.strip()

df_final = df_final.merge(
    df_amp[["BioProject", "Seq", "Expected_amplicon_length", "Expected_amplicon_length_cutoff"]],
    on=["BioProject", "Seq"],
    how="left"
)

df_final['mean_qs20_af'] = pd.to_numeric(df_final['mean_qs20_af'], errors='coerce')
df_final['Expected_amplicon_length_cutoff'] = pd.to_numeric(df_final['Expected_amplicon_length_cutoff'], errors='coerce')

# ==============================
# 7. Test 열 생성
# ==============================
df_final['Test'] = np.where(df_final['mean_qs20_af'] >= df_final['Expected_amplicon_length_cutoff']*0.9, "pass", "fail")

# ==============================
# 8. Trunc1, Trunc2 (fail도 포함)
# ==============================
df_final['Trunc1'] = ""
df_final['Trunc2'] = ""

mask_pe = (df_final['Layout'] == "PE")
df_final.loc[mask_pe, 'Trunc1'] = ((df_final.loc[mask_pe, 'Expected_amplicon_length_cutoff'] + 20) / 2).astype(int)
df_final.loc[mask_pe, 'Trunc2'] = ((df_final.loc[mask_pe, 'Expected_amplicon_length_cutoff'] - 20) / 2).astype(int)

mask_se = (df_final['Layout'] == "SE")
df_final.loc[mask_se, 'Trunc1'] = df_final.loc[mask_se, 'Expected_amplicon_length_cutoff'].astype(int)
df_final.loc[mask_se, 'Trunc2'] = ""

# ==============================
# 9. 1차 결과 저장
# ==============================
df_final.to_csv(output_file_step2, index=False)
print(f"1차 완료! 저장된 파일: {output_file_step2}")

# ==============================
# 10. fail SRR 추출
# ==============================
fail_srrs = df_final.loc[df_final['Test'] == "fail", 'SRR'].unique().tolist()

# ==============================
# 11. fail 파일 이동 + 정보 기록
# ==============================
def move_fail_files(src_dir, dst_dir, fail_srrs):
    moved_files = []  # (SRR, file_name, src_path, dst_path)
    for root, dirs, files in os.walk(src_dir):
        for file in files:
            for srr in fail_srrs:
                if srr in file:
                    src_path = Path(root) / file
                    rel_path = Path(root).relative_to(src_dir)
                    dst_path = Path(dst_dir) / rel_path / file
                    dst_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(src_path), str(dst_path))
                    moved_files.append((srr, file, str(src_path), str(dst_path)))
                    break
    return moved_files

moved_16s = move_fail_files(fastq_dir_16s, fastq_removed, fail_srrs)
moved_its = move_fail_files(fastq_dir_its, fastq_removed, fail_srrs)
moved_all = moved_16s + moved_its

print(f"\n옮겨진 파일 수: {len(moved_all)}")

# ==============================
# 12. fail 제거 후 최종 메타데이터 저장
# ==============================
df_final_filtered = df_final[df_final['Test'] != "fail"].copy()
df_final_filtered.to_csv(output_file_final, index=False)
print(f"최종 완료! 저장된 파일: {output_file_final}")

# ==============================
# 13. 옮겨진 파일 정보 출력
# ==============================
if moved_all:
    print("\n[FAIL된 SRR 및 옮겨진 파일 정보]")
    fail_info = df_final[df_final['Test'] == "fail"][["SRR","BioProject","Seq","Layout"]].drop_duplicates()
    fail_info = fail_info.set_index("SRR")

    for srr, fname, src_path, dst_path in moved_all:
        if srr in fail_info.index:
            meta = fail_info.loc[srr]
            print(f"SRR: {srr}\tBioProject: {meta['BioProject']}\tSeq: {meta['Seq']}\tLayout: {meta['Layout']}\tFile: {fname}")

# ==============================
# 14. 파일 개수 출력
# ==============================
def count_fastq_files(base_dir):
    count = 0
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".fastq") or file.endswith(".fastq.gz"):
                count += 1
    return count

count_16s = count_fastq_files(fastq_dir_16s)
count_its = count_fastq_files(fastq_dir_its)
count_removed = count_fastq_files(fastq_removed)

print("\n파일 개수 요약:")
print(f"- {fastq_dir_16s}: {count_16s}개")
print(f"- {fastq_dir_its}: {count_its}개")
print(f"- {fastq_removed}: {count_removed}개 (옮겨진 fail fastq 파일)")
