# -----------------------------
# Make phyloseq objects from Genus / Species filtered CSV files (Dynamic min_reads)
# -----------------------------
library(phyloseq)
library(dplyr)
library(readr)
library(stringr)
# -----------------------------
# Directories
# -----------------------------
genus_dir   <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq_genus_16s_filtered/"
species_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq_species_16s_filtered/"
# -----------------------------
# Helper function: CSV -> phyloseq with dynamic min_reads
# -----------------------------
make_phyloseq_from_csv <- function(file, level=c("genus","species"), min_reads) {
  level <- match.arg(level)
  dat <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  # Identify taxonomy vs sample columns
  if (level == "genus") {
    tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")
  } else {
    tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  }
  sample_cols <- setdiff(colnames(dat), tax_cols)
  # -----------------------------
  # ① taxonomy 불완전(NA, "NA", 빈 값) ASV 제거
  # -----------------------------
  dat <- dat %>%
    filter(
      !if_any("Phylum",
              ~ is.na(.) | . == "" | . == "NA")
    )
  # -----------------------------
  # ② OTU table / taxonomy table 생성
  # -----------------------------
  otu_mat <- as.matrix(dat[, sample_cols])
  rownames(otu_mat) <- apply(dat[, tax_cols], 1, function(x) paste(x, collapse="|"))
  otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)
  tax_mat <- as.matrix(dat[, tax_cols])
  rownames(tax_mat) <- rownames(otu_mat)
  tax_tab <- tax_table(tax_mat)
  phy <- phyloseq(otu_tab, tax_tab)
  # -----------------------------
  # ③ 샘플 read count < min_reads 제거 (dynamic threshold)
  # -----------------------------
  phy <- prune_samples(sample_sums(phy) >= min_reads, phy)
  return(phy)
}
# -----------------------------
# Load metadata
# -----------------------------
meta_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/scalp_16s.csv"
metadata <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta_cols <- c("SRR","BioSample","BioProject","Subject_id","Dx","Lesion",
               "Age","Sex","Ethnicity","Sample","Site","Scalp","Time_points")
metadata <- metadata[, meta_cols]
# -----------------------------
# Function: add metadata to phyloseq object
# -----------------------------
add_metadata_to_phyloseq <- function(phy, metadata){
  sample_names_phy <- sample_names(phy)
  meta_sub <- metadata %>%
    filter(SRR %in% sample_names_phy) %>%
    distinct(SRR, .keep_all = TRUE) %>%
    tibble::column_to_rownames("SRR")
  meta_sub <- meta_sub[sample_names_phy, , drop = FALSE]
  samp <- sample_data(meta_sub)
  phy <- merge_phyloseq(phy, samp)
  return(phy)
}
# -----------------------------
# Helper: convert phyloseq object to relative abundance
# -----------------------------
make_relabund <- function(phy) {
  transform_sample_counts(phy, function(x) x / sum(x))
}
# -----------------------------
# Helper: filter empty samples (모든 OTU=0인 샘플 제거)
# -----------------------------
filter_empty_samples <- function(ps_obj, output_csv) {
  otu_mat <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) {
    otu_mat <- t(otu_mat)
  }
  is_all_zero_or_empty <- function(x) {
    all(is.na(x) | x == 0 | x == "")
  }
  bad_samples <- colnames(otu_mat)[apply(otu_mat, 2, is_all_zero_or_empty)]
  if (length(bad_samples) > 0) {
    meta_bad <- as.data.frame(sample_data(ps_obj)[bad_samples, ])
    meta_bad <- tibble::rownames_to_column(meta_bad, var = "SampleID")
    write_csv(meta_bad, output_csv)
  } else {
    message("No problematic samples found for ", output_csv)
  }
  ps_filtered <- prune_samples(!(sample_names(ps_obj) %in% bad_samples), ps_obj)
  return(ps_filtered)
}
# -----------------------------
# Pipeline 실행 함수: now min_reads is dynamic
# -----------------------------
process_phyloseq_dir <- function(input_dir, level, metadata, output_csv, min_reads) {
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  objs <- list()
  for (f in files) {
    objname <- tools::file_path_sans_ext(basename(f))
    message("Processing: ", objname)
    phy <- make_phyloseq_from_csv(f, level=level, min_reads=min_reads)
    if (nsamples(phy) == 0) next
    phy <- add_metadata_to_phyloseq(phy, metadata)
    phy <- make_relabund(phy)
    objs[[objname]] <- phy
  }
  if (length(objs) > 1) {
    merged <- do.call(merge_phyloseq, objs)
  } else {
    merged <- objs[[1]]
  }
  merged_filtered <- filter_empty_samples(merged, output_csv)
  return(merged_filtered)
}
# -----------------------------
# 실행: set min_reads as needed!
# -----------------------------
my_min_reads <- 50  # Set to your desired value

genus_csv   <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/no_read_count_sample_genus.csv"
species_csv <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/no_read_count_sample_species.csv"
merged_16s_genus_rel_filtered <- process_phyloseq_dir(genus_dir, "genus", metadata, genus_csv, min_reads=my_min_reads)
merged_16s_species_rel_filtered <- process_phyloseq_dir(species_dir, "species", metadata, species_csv, min_reads=my_min_reads)

# 결과 확인
merged_16s_genus_rel_filtered
merged_16s_species_rel_filtered
