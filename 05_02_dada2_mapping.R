# Taxonomy Mapping + Phyloseq (Error-safe)
library(dada2)
library(stringr)
library(phyloseq)

# Input/Output directories
input_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/"
output_dir <- file.path(input_dir, "taxo_phyloseq")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Database paths
silva_genus <- "/home/jwchoi/projects/microbiome/databases/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species <- "/home/jwchoi/projects/microbiome/databases/silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
silva_species_assign <- "/home/jwchoi/projects/microbiome/databases/silva/silva_v138.2_assignSpecies.fa.gz"
unite_species <- "/home/jwchoi/projects/microbiome/databases/unite/sh_general_release_dynamic_19.02.2025.fasta"

# List input files
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  message("Processing: ", f)
  fname <- basename(f)
  is_16s <- str_detect(tolower(fname), "16s")
  is_its <- str_detect(tolower(fname), "its")
  
  if (!is_16s && !is_its) {
    warning("Skipping file (not 16S or ITS): ", fname)
    next
  }
  
  # Load ASV table (samples = rows, ASVs = columns)
  asv_tab <- read.csv(f, check.names = FALSE, row.names = 1)
  
  if (ncol(asv_tab) == 0) {
    warning("No ASVs found in file: ", fname)
    next
  }
  
  asv_seqs <- colnames(asv_tab)
  otu_mat <- t(as.matrix(asv_tab)) # ASVs = rows
  
  # Safe taxonomy assignment
  safe_assign <- function(seqs, refFasta) {
    tryCatch({
      assignTaxonomy(seqs, refFasta = refFasta, multithread = TRUE)
    }, error = function(e) {
      warning("Taxonomy assignment failed: ", e$message)
      return(matrix(NA, nrow = length(seqs), ncol = 1))
    })
  }
  
  if (is_16s) {
    tax_genus <- safe_assign(asv_seqs, silva_genus)
    tax_species <- safe_assign(asv_seqs, silva_species)
    tax_species <- tryCatch({
      addSpecies(tax_species, silva_species_assign, verbose = TRUE)
    }, error = function(e) {
      warning("Species assignment failed: ", e$message)
      return(tax_species)
    })
    
    # Save taxonomy + sample counts
    genus_out <- file.path(output_dir, gsub("\\.csv$", "_tax_genus.csv", fname))
    species_out <- file.path(output_dir, gsub("\\.csv$", "_tax_species.csv", fname))
    
    write.csv(cbind(ASV = asv_seqs, tax_genus, t(asv_tab)), genus_out, row.names = FALSE)
    write.csv(cbind(ASV = asv_seqs, tax_species, t(asv_tab)), species_out, row.names = FALSE)
    
    # Create phyloseq objects
    ps_genus <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
                         tax_table(as.matrix(tax_genus)))
    ps_species <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
                           tax_table(as.matrix(tax_species)))
    
    saveRDS(ps_genus, file = file.path(output_dir, gsub("\\.csv$", "_phyloseq_genus.rds", fname)))
    saveRDS(ps_species, file = file.path(output_dir, gsub("\\.csv$", "_phyloseq_species.rds", fname)))
    
  } else if (is_its) {
    tax_species <- safe_assign(asv_seqs, unite_species)
    
    species_out <- file.path(output_dir, gsub("\\.csv$", "_tax_species.csv", fname))
    write.csv(cbind(ASV = asv_seqs, tax_species, t(asv_tab)), species_out, row.names = FALSE)
    
    ps_species <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE),
                           tax_table(as.matrix(tax_species)))
    saveRDS(ps_species, file = file.path(output_dir, gsub("\\.csv$", "_phyloseq_species.rds", fname)))
  }
  
  message("Finished: ", fname)
}
