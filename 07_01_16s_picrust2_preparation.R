# -----------------------------
# 16S ASV-level Filtering + PICRUSt2 Input + Merge Stats
# -----------------------------

library(dplyr)
library(readr)
library(Biostrings)

# -----------------------------
# Input / Output directories
# -----------------------------
input_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq/"
output_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/picrust2_input_merged/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

stats_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/picrust2_merge_stats.csv"

# -----------------------------
# biom executable path (절대경로)
# -----------------------------
biom_path <- "/home/jwchoi/miniconda3/bin/biom"

# -----------------------------
# List files
# -----------------------------
files <- list.files(input_dir, pattern = "16s.*species.*\\.csv$", full.names = TRUE)

# -----------------------------
# User-defined filtering thresholds
# -----------------------------
min_reads_total <- 50
min_prevalence <- 0.01

# -----------------------------
# Initialize
# -----------------------------
all_abundance <- list()
all_seqs <- DNAStringSet()
stats_list <- list()

# -----------------------------
# Loop through files
# -----------------------------
for(f in files){
  fname <- basename(f)
  message("Processing: ", fname)
  
  dat <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  
  # -----------------------------
  # 중복 Species 컬럼 처리
  # -----------------------------
  species_idx <- which(colnames(dat) == "Species")
  if(length(species_idx) > 1){
    colnames(dat)[species_idx[1:2]] <- c("Species1","Species2")
    dat$Species <- ifelse(dat$Species1 != "NA", dat$Species1,
                          ifelse(dat$Species2 != "NA", dat$Species2, "NA"))
    dat <- dat[, !(colnames(dat) %in% c("Species1","Species2"))]
  }
  
  # -----------------------------
  # Identify sample columns
  # -----------------------------
  sample_cols <- setdiff(colnames(dat), c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  
  # -----------------------------
  # Pre-filter unwanted taxa
  # -----------------------------
  dat_filtered <- dat %>%
    filter(!(Class %in% c("Chloroplast","Mitochondria") |
               Order %in% c("Chloroplast","Mitochondria") |
               Family %in% c("Chloroplast","Mitochondria") |
               Genus %in% c("Chloroplast","Mitochondria"))) %>%
    filter(Kingdom == "Bacteria") %>%
    filter(!(Phylum=="NA" & Class=="NA" & Order=="NA" & Family=="NA" & Genus=="NA" & Species=="NA"))
  
  # -----------------------------
  # Low abundance / low prevalence filtering
  # -----------------------------
  dat_filtered <- dat_filtered %>%
    rowwise() %>%
    mutate(
      total_reads = sum(c_across(all_of(sample_cols)), na.rm=TRUE),
      prevalence = sum(c_across(all_of(sample_cols)) > 0)/length(sample_cols)
    ) %>%
    ungroup() %>%
    filter(total_reads >= min_reads_total, prevalence >= min_prevalence) %>%
    dplyr::select(-total_reads, -prevalence)
  
  
  # -----------------------------
  # Statistics
  # -----------------------------
  stats_list[[fname]] <- data.frame(
    file = fname,
    n_samples = length(sample_cols),
    n_ASVs = nrow(dat_filtered),
    total_reads = sum(dat_filtered[, sample_cols], na.rm = TRUE)
  )
  
  # -----------------------------
  # Append to global lists
  # -----------------------------
  all_abundance[[fname]] <- dat_filtered[, c("ASV", sample_cols)]
  
  seqs <- DNAStringSet(dat_filtered$ASV)
  names(seqs) <- dat_filtered$ASV
  all_seqs <- append(all_seqs, seqs)
  
  message("File processed: ", fname,
          " | ASVs: ", nrow(dat_filtered),
          " | Samples: ", length(sample_cols))
}

# -----------------------------
# Merge all abundance tables
# -----------------------------
merged_abundance <- bind_rows(all_abundance)

# -----------------------------
# ASV별 합산 (중복 제거)
# -----------------------------
merged_abundance <- merged_abundance %>%
  group_by(ASV) %>%
  summarise(across(everything(), sum, na.rm=TRUE)) %>%
  ungroup()

# -----------------------------
# FASTA 중복 제거
# -----------------------------
all_seqs <- all_seqs[!duplicated(names(all_seqs))]

# -----------------------------
# Replace NA with 0 and force numeric
# -----------------------------
merged_abundance[is.na(merged_abundance)] <- 0
merged_abundance[, -1] <- lapply(merged_abundance[, -1], as.numeric)

# -----------------------------
# Save merged CSV and TSV
# -----------------------------
out_csv <- file.path(output_dir, "ASV_filtered_all.csv")
write.csv(merged_abundance, out_csv, row.names=FALSE)

out_tsv <- file.path(output_dir, "ASV_filtered_all.tsv")
write.table(merged_abundance, out_tsv, sep="\t", quote=FALSE, row.names=FALSE)

# -----------------------------
# Save merged FASTA
# -----------------------------
out_fasta <- file.path(output_dir, "ASV_filtered_all.fasta")
writeXStringSet(all_seqs, out_fasta)

# -----------------------------
# BIOM 생성
# -----------------------------
out_biom <- file.path(output_dir, "ASV_filtered_all.biom")
cmd <- paste(biom_path, "convert -i", out_tsv, "-o", out_biom,
             "--to-hdf5 --table-type='OTU table'")
system(cmd)
message("BIOM file created: ", out_biom)

# -----------------------------
# Merge statistics and save
# -----------------------------
merge_stats <- do.call(rbind, stats_list)
merge_stats$merged_ASVs <- nrow(merged_abundance)
merge_stats$merged_samples <- ncol(merged_abundance) - 1
write.csv(merge_stats, stats_file, row.names = FALSE)
message("Merge statistics saved to: ", stats_file)

# -----------------------------
# PICRUSt2 실행 예제
# -----------------------------
message("\nYou can now run PICRUSt2 with the following command (adjust threads -p):")
message("picrust2_pipeline.py -s ", out_fasta, 
        " -i ", out_biom,
        " -o ", file.path(output_dir, "picrust2_out"),
        " -p 4")
