# -----------------------------
# 16S Species & Genus Level Filtering + Aggregation with additional filtering
# -----------------------------

library(dplyr)
library(readr)

# -----------------------------
# Input / Output directories
# -----------------------------
input_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq/"

output_species <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq_species_16s_filtered/"
output_genus   <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/taxo_phyloseq_genus_16s_filtered/"

if(!dir.exists(output_species)) dir.create(output_species, recursive = TRUE)
if(!dir.exists(output_genus))   dir.create(output_genus, recursive = TRUE)

stats_species_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_asv_16s_species_filtered.csv"
stats_genus_file   <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_asv_16s_genus_filtered.csv"

# -----------------------------
# List only files containing "16s" and "species"
# -----------------------------
files <- list.files(input_dir, pattern = "16s.*species.*\\.csv$", full.names = TRUE)

stats_species_list <- list()
stats_genus_list <- list()

# -----------------------------
# User-defined filtering thresholds
# -----------------------------
min_reads_total <-  50     # 최소 read count
min_prevalence <- 0.01       # 최소 존재 샘플 비율 (5%)

for(f in files){
  message("Processing: ", f)
  fname <- basename(f)
  bioproject <- sub("_.*", "", fname)   # 파일명에서 BioProject 추출
  
  dat <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  
  # -----------------------------
  # Species1(8), Species2(9) -> Species 하나로 통합
  # -----------------------------
  colnames(dat)[8:9] <- c("Species1", "Species2")
  dat$Species <- ifelse(dat$Species1 != "NA", dat$Species1,
                        ifelse(dat$Species2 != "NA", dat$Species2, "NA"))
  dat <- dat[, !(colnames(dat) %in% c("Species1","Species2"))]
  
  # -----------------------------
  # Identify sample columns
  # -----------------------------
  sample_cols <- setdiff(colnames(dat), c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  
  # -----------------------------
  # Pre-filter unwanted taxa
  # 1. Chloroplast / Mitochondria 제거
  # 2. Kingdom == "Bacteria"만 남기기
  # 3. Phylum~Species 모두 NA면 제거
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
    mutate(total_reads = sum(c_across(all_of(sample_cols)), na.rm=TRUE),
           prevalence = sum(c_across(all_of(sample_cols)) > 0)/length(sample_cols)) %>%
    ungroup() %>%
    filter(total_reads >= min_reads_total, prevalence >= min_prevalence) %>%
    select(-total_reads, -prevalence)
  
  # -----------------------------
  # Species-level aggregation
  # -----------------------------
  dat_species <- dat_filtered %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    filter(rowSums(across(all_of(sample_cols))) > 0)
  
  total_reads_species <- sum(dat_species[, sample_cols], na.rm = TRUE)
  total_asv_species <- nrow(dat_species)
  
  fname_species <- sub("_tax_species\\.csv$", "_tax_species.csv", fname)
  out_species <- file.path(output_species, fname_species)
  write.csv(dat_species, out_species, row.names = FALSE)
  
  stats_species_list[[fname]] <- data.frame(
    file = fname_species,
    BioProject = bioproject,
    total_asv_before = nrow(dat),
    total_reads_before = sum(dat[, sample_cols], na.rm=TRUE),
    total_asv_after = total_asv_species,
    total_reads_after = total_reads_species
  )
  
  # -----------------------------
  # Genus-level aggregation
  # -----------------------------
  dat_genus <- dat_filtered %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>%
    summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    filter(rowSums(across(all_of(sample_cols))) > 0)
  
  total_reads_genus <- sum(dat_genus[, sample_cols], na.rm = TRUE)
  total_asv_genus <- nrow(dat_genus)
  
  fname_genus <- sub("_tax_species\\.csv$", "_tax_genus.csv", fname)
  out_genus <- file.path(output_genus, fname_genus)
  write.csv(dat_genus, out_genus, row.names = FALSE)
  
  stats_genus_list[[fname]] <- data.frame(
    file = fname_genus,
    BioProject = bioproject,
    total_asv_before = nrow(dat),
    total_reads_before = sum(dat[, sample_cols], na.rm=TRUE),
    total_asv_after = total_asv_genus,
    total_reads_after = total_reads_genus
  )
  
  message("Finished: ", fname,
          " | Species ASV: ", total_asv_species,
          " | Genus ASV: ", total_asv_genus,
          " | Reads(Species): ", total_reads_species,
          " | Reads(Genus): ", total_reads_genus)
}

# -----------------------------
# Save statistics files
# -----------------------------
stats_species_df <- do.call(rbind, stats_species_list)
write.csv(stats_species_df, stats_species_file, row.names = FALSE)
message("Species statistics saved to: ", stats_species_file)

stats_genus_df <- do.call(rbind, stats_genus_list)
write.csv(stats_genus_df, stats_genus_file, row.names = FALSE)
message("Genus statistics saved to: ", stats_genus_file)
