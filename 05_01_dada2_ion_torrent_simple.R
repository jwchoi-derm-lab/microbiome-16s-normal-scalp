# DADA2 Ion Torrent Analysis - Stringent Parameters to Reduce ASV Inflation (Clean Output)
library(dada2)
library(data.table)
library(dplyr)

# -----------------------------
# Paths
# -----------------------------
fastq_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_trimmed_dada2/se/PRJNA1115970" # 바이오프로젝트 시퀀스 파일 위치 
filtered_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/IonT/filtered"
results_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/IonT/results"
summary_file <- file.path(results_dir, "summary.csv")

dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Collect FASTQ files
# -----------------------------
fastq_files <- list.files(fastq_dir, pattern = "\\.fastq|\\.fastq\\.gz$", full.names = TRUE)
sample_names <- gsub("_s.*fastq.*$", "", basename(fastq_files))
names(fastq_files) <- sample_names

cat(sprintf("[INFO] Found %d FASTQ files.\n", length(fastq_files)))

# -----------------------------
# Quality assessment (first 6 samples)
# -----------------------------
pdf(file.path(results_dir, "quality_profiles.pdf"))
plotQualityProfile(fastq_files[1:min(6, length(fastq_files))])
dev.off()

# -----------------------------
# Filter & Trim (stringent)
# -----------------------------
filtered_files <- file.path(filtered_dir, paste0(sample_names, "_filt.fastq.gz"))
out <- filterAndTrim(fastq_files, filtered_files,
                     truncLen = 150, maxN = 0, maxEE = 1, truncQ = 20,
                     minLen = 145, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

keep <- out[, "reads.out"] >= 1000
filtered_files <- filtered_files[keep]
sample_names <- sample_names[keep]
names(filtered_files) <- sample_names

if(length(filtered_files) == 0) stop("[ERROR] No samples passed filtering.")

cat(sprintf("[INFO] Retained %d samples after filtering.\n", length(filtered_files)))

# -----------------------------
# Learn error rates
# -----------------------------
errF <- learnErrors(filtered_files, multithread = TRUE, nbases = 5e7,
                    randomize = TRUE, MAX_CONSIST = 15, qualityType = "FastqQuality", OMEGA_C = 0)

pdf(file.path(results_dir, "error_rates.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

# -----------------------------
# Dereplication & DADA2 inference
# -----------------------------
dadaFs <- list()
for (i in seq_along(filtered_files)) {
  drp <- derepFastq(filtered_files[i], verbose = FALSE)
  dadaFs[[i]] <- dada(drp, err = errF, multithread = FALSE,
                      HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 16,
                      OMEGA_A = 1e-20, OMEGA_C = 1e-40,
                      DETECT_SINGLETONS = FALSE, GAPLESS = TRUE,
                      MAX_CLUST = 0, MIN_FOLD = 2, MIN_HAMMING = 1, pool = FALSE)
}
names(dadaFs) <- names(filtered_files)

# -----------------------------
# Make sequence table
# -----------------------------
seqtab <- makeSequenceTable(dadaFs)

seq_lengths <- nchar(getSequences(seqtab))
length_filter <- seq_lengths >= 100 & seq_lengths <= 300
seqtab_length <- seqtab[, length_filter]

# -----------------------------
# Chimera removal
# -----------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab_length, method = "consensus", multithread = TRUE,
                                    verbose = FALSE, minFoldParentOverAbundance = 2.0,
                                    allowOneOff = FALSE, minOneOffParentDistance = 2)

# -----------------------------
# Abundance & prevalence filtering
# -----------------------------
asv_counts <- colSums(seqtab.nochim)
total_reads <- sum(seqtab.nochim)
min_reads_relative <- max(10, round(total_reads * 0.00001))
min_reads_absolute <- 50
min_abundance <- max(min_reads_relative, min_reads_absolute)

seqtab.filtered <- seqtab.nochim[, asv_counts >= min_abundance]

asv_prevalence <- colSums(seqtab.filtered > 0)
prev_threshold <- 2
seqtab.final <- seqtab.filtered[, asv_prevalence >= prev_threshold]

final_asv_counts <- colSums(seqtab.final)
asvs_per_sample <- rowSums(seqtab.final > 0)

# -----------------------------
# Track reads
# -----------------------------
getN <- function(x) sum(getUniques(x))
track <- data.frame(Sample = sample_names,
                    input = out[keep, "reads.in"],
                    filtered = out[keep, "reads.out"],
                    denoised = sapply(dadaFs, getN),
                    length_filt = rowSums(seqtab_length),
                    nonchim = rowSums(seqtab.nochim),
                    abundance_filt = rowSums(seqtab.filtered),
                    final = rowSums(seqtab.final),
                    pct_retained = round(rowSums(seqtab.final)/out[keep, "reads.in"]*100, 1))
write.csv(track, file.path(results_dir, "read_tracking.csv"), row.names = FALSE)

# -----------------------------
# Save ASV sequences (FASTA)
# -----------------------------
sequences <- getSequences(seqtab.final)
asv_names <- paste0("ASV_", sprintf("%04d", 1:length(sequences)))

fasta_file <- file.path(results_dir, "ASV_sequences_final.fasta")
cat("", file = fasta_file)
for(i in seq_along(sequences)) {
  cat(paste0(">", asv_names[i], "\n"), file = fasta_file, append = TRUE)
  cat(paste0(sequences[i], "\n"), file = fasta_file, append = TRUE)
}

# ASV count table
seqtab_named <- seqtab.final
colnames(seqtab_named) <- asv_names
write.csv(data.frame(Sample = rownames(seqtab_named), seqtab_named),
          file.path(results_dir, "ASV_count_table.csv"), row.names = FALSE)

# RDS
saveRDS(seqtab.final, file.path(results_dir, "seqtab_final.rds"))

# ASV summary
asv_summary <- data.frame(
  ASV_ID = asv_names,
  Sequence = sequences,
  Length = nchar(sequences),
  Total_reads = final_asv_counts,
  Prevalence = colSums(seqtab.final > 0),
  Max_abundance = apply(seqtab.final, 2, max),
  stringsAsFactors = FALSE
)
write.csv(asv_summary, file.path(results_dir, "ASV_summary.csv"), row.names = FALSE)

# -----------------------------
# Final summary
# -----------------------------
summary_df <- data.frame(
  Dataset = "PRJNA1115970",
  Analysis = "Ion_Torrent_stringent",
  Total_samples = length(fastq_files),
  Samples_retained = nrow(seqtab.final),
  Input_reads = sum(out[keep, "reads.in"]),
  Final_reads = sum(seqtab.final),
  Retention_rate_pct = round(sum(seqtab.final)/sum(out[keep, "reads.in"])*100, 1),
  Total_ASVs = ncol(seqtab.final),
  Median_ASVs_per_sample = median(asvs_per_sample),
  Mean_ASVs_per_sample = round(mean(asvs_per_sample), 1),
  Min_abundance_threshold = min_abundance,
  Min_prevalence_threshold = prev_threshold,
  Chimera_rate_pct = round(100*(1-sum(seqtab.nochim)/sum(seqtab_length)),1),
  stringsAsFactors = FALSE
)
write.csv(summary_df, summary_file, row.names = FALSE)

# -----------------------------
# Clean up unnecessary files
# -----------------------------
# Keep only filtered fastq (compressed), results CSV, FASTA, RDS, PDF
# Remove intermediate variables if needed
cat("[INFO] Analysis complete. Final outputs saved to results directory.\n")
cat(sprintf("  Samples: %d, ASVs: %d\n", nrow(seqtab.final), ncol(seqtab.final)))


# 파일 경로
count_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/test/results/ASV_count_table.csv"
summary_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/test/results/ASV_summary.csv"
out_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/test/results/16s_SE_PRJNA1115970.csv"

# 데이터 불러오기
asv_counts <- fread(count_file)
asv_summary <- fread(summary_file)

# ASV ID -> Sequence 매핑 벡터 만들기
asv_map <- setNames(asv_summary$Sequence, asv_summary$ASV_ID)

# count table의 ASV 컬럼 이름을 Sequence로 치환
new_colnames <- colnames(asv_counts)
new_colnames[-1] <- asv_map[new_colnames[-1]]   # 첫 컬럼(Sample)은 그대로 두고 나머지만 치환
colnames(asv_counts) <- new_colnames

# 결과 저장
fwrite(asv_counts, out_file)

cat("✅ 변환된 파일이 저장되었습니다:", out_file, "\n")
