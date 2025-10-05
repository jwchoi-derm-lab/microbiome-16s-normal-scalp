library(ggforce)
library(phyloseq)
library(tidyr)
library(scales)
library(microbiome)
library(microbiomeutilities)
library(tidyverse)
library(reshape2)
library(ggalluvial)
library(forcats)
library(ggforce)

setwd("/media/jwchoi/ssd2/projects/microbiome/16S_1/R")

metadata_finder = function(df, metadata){
  BioProject = NULL
  Sample = NULL
  Scalp = NULL
  for (i in as.character(df$variable)){
    BioProject = append(BioProject, metadata$BioProject[metadata$BioSample == i])
    Sample <- append(Sample, metadata$Sample[metadata$BioSample == i])
    Scalp = append(Scalp, metadata$Scalp[metadata$BioSample== i])
  }
  df$BioProject <- BioProject
  df$Sample <- Sample
  df$Scalp <- Scalp
  return(df)
}
'%ni%' <- Negate('%in%')

TargetLvl <- function(df, mytab){
  joinby <- 'BioSample'
  spn<-3
  # original legend order preserved:
  target_taxa <- c(
    'Actinomycetota', 
    paste(paste(rep(" ", spn), collapse = ""),  'Cutibacterium'),
    paste(paste(rep(" ", spn), collapse = ""),  'Corynebacterium'),
    paste(paste(rep(" ", spn), collapse = ""),  'Lawsonella'),
    'Bacillota', 
    paste(paste(rep(" ", spn), collapse = ""),  'Staphylococcus'),
    paste(paste(rep(" ", spn), collapse = ""),  'Streptococcus'), 
    'Pseudomonadota',
    paste(paste(rep(" ", spn), collapse = ""),  'Acinetobacter'),
    paste(paste(rep(" ", spn), collapse = ""),  'Pseudomonas'), 
    paste(paste(rep(" ", spn), collapse = ""),  'Enhydrobacter'), 
    'Bacteroidota',
    'Others'
  )
  df_origin <- df
  df$Species <- gsub('_', ' ', df$Species)
  df <- df[!(df$Kingdom %in% c('Archaea','Viruses','Eukaryota')),]
  taxa_only <- df[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
  df <- melt(df, id=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  df$Target <- df$Kingdom
  df$Target[df$Kingdom == 'Bacteria'] <- 'Others'
  df$Target[df$Phylum == 'Actinomycetota'] <- 'Actinomycetota'
  df$Target[df$Phylum == 'Bacillota'] <- 'Bacillota'
  df$Target[df$Phylum == 'Pseudomonadota'] <- 'Pseudomonadota'
  df$Target[df$Phylum == 'Bacteroidota'] <- 'Bacteroidota'
  df$Target[df$Genus == 'Cutibacterium'] <- paste(paste(rep(" ", spn), collapse = ""), 'Cutibacterium')
  df$Target[df$Genus == 'Corynebacterium'] <- paste(paste(rep(" ", spn), collapse = ""),  'Corynebacterium')
  df$Target[df$Genus == 'Lawsonella'] <- paste(paste(rep(" ", spn), collapse = ""), 'Lawsonella')
  df$Target[df$Genus == 'Staphylococcus'] <- paste(paste(rep(" ", spn), collapse = ""),  'Staphylococcus')
  df$Target[df$Genus == 'Streptococcus'] <- paste(paste(rep(" ", spn), collapse = ""),  'Streptococcus')
  df$Target[df$Genus == 'Acinetobacter'] <- paste(paste(rep(" ", spn), collapse = ""),  'Acinetobacter')
  df$Target[df$Genus == 'Pseudomonas'] <- paste(paste(rep(" ", spn), collapse = ""),  'Pseudomonas')
  df$Target[df$Genus == 'Enhydrobacter'] <- paste(paste(rep(" ", spn), collapse = ""),  'Enhydrobacter')
  df$Target <- factor(df$Target, levels=target_taxa)
  df$Kingdom <- NULL
  df$Phylum <- NULL
  df$Class <- NULL
  df$Order <- NULL
  df$Family <- NULL
  df$Genus <- NULL
  df$Species <- NULL
  df <- data.frame(df %>% group_by(variable, Target) %>% summarise_all(list(sum)))
  df <- metadata_finder(df, mytab)
  return(df)
}

# 1. Select, clean, normalize data
merged_16s_species_rel_filtered_subset <- merged_16s_species_rel_filtered %>%
  subset_samples(Dx == "CTRL" & Time_points == "t=1 (Day 0)") %>%
  subset_taxa(Genus != "Methylobacterium") %>%
  subset_taxa(Genus != "Xanthomonas") %>%
    subset_taxa(!(Species %in% c("algae",
                                 "alkaliphila",
                                 "aurantiacum",
                                 "equorum",
                                 "grossensis",
                                 "ureicelerivorans",
                                 "vesicularis",
                                 "yamanorum",
                                 "rhizophila")))

ra_rel.bac <- prune_taxa(taxa_sums(merged_16s_species_rel_filtered_subset)>0, merged_16s_species_rel_filtered_subset)
ra_rel.bac <- transform_sample_counts(ra_rel.bac, function(x) x/sum(x))
df <- cbind(data.frame(tax_table(ra_rel.bac)), data.frame(otu_table(ra_rel.bac)))
my_meta <- sample_data(ra_rel.bac)
my_meta <- as.data.frame(my_meta)
my_meta[,"BioSample"] <- rownames(my_meta) 
rownames(my_meta) <- NULL 
df_sorted <- df[order(rowSums(df[, 8:ncol(df)]), decreasing=TRUE), ]
df <- TargetLvl(df_sorted, my_meta)

# 2. Color palette for taxa (match above order!)
color_pall <- c(
  "#000080",  # Actinomycetota
  "#1560BD",  # Cutibacterium
  "#1E90FF",  # Corynebacterium
  "#87CEEB",  # Lawsonella
  "#8B0000",  # Bacillota
  "#FF0000",  # Staphylococcus
  "#F08080",  # Streptococcus
  "#6B8E23",  # Pseudomonadota 
  "#228B22",  # Acinetobacter
  "#00FF00",  # Pseudomonas
  "#7FFF00",  # Enhydrobacter
  "#FFD700",  # Bacteroidota
  "#D3D3D3"   # Others
)

# 3. Reorder Sample and BioProject factors
factor_levels <- unique(df$Sample)  
other_factors <- c("Hair", "Swab", "Tissue")
Sample_order <- c("Hair", "Swab", "Tissue")
df$Sample <- factor(df$Sample, levels = Sample_order)
df$BioProject <- factor(df$BioProject, levels = unique(df$BioProject)) 
df$BioProject <- as.character(df$BioProject)

# 4. Summarize for stack bar plot
df_by_BioProject <- df %>%
  dplyr::group_by(BioProject, Sample, Target) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(BioProject, Sample) %>%
  dplyr::mutate(scaling_factor = 0.9999999 / sum(value, na.rm = TRUE),
                value = value * scaling_factor) %>%
  dplyr::ungroup() %>%
  dplyr::select(-scaling_factor)

# 5. Compute Cutibacterium abundance for ordering
cuti_abundance <- df_by_BioProject %>%
  filter(Target == "    Cutibacterium") %>%
  select(BioProject, Sample, cuti_value = value)

# 6. Join abundance back and create compound key for unique x-axis
df_by_BioProject <- df_by_BioProject %>%
  left_join(cuti_abundance, by = c("BioProject", "Sample")) %>%
  mutate(cuti_value = ifelse(is.na(cuti_value), 0, cuti_value),
         Sample_BioProject = paste0(Sample, "_", BioProject))

# 7. Reorder within each Sample group descending by Cutibacterium abundance
df_by_BioProject <- df_by_BioProject %>%
  group_by(Sample) %>%
  mutate(
    Sample_BioProject = forcats::fct_reorder(Sample_BioProject, cuti_value, .desc = TRUE)
  ) %>%
  ungroup()

# 7b. Create new factor just for stacking order (Cutibacterium at bottom)
df_by_BioProject <- df_by_BioProject %>%
  mutate(Target_plot = forcats::fct_relevel(Target, "    Cutibacterium", after = 0))


# 8. Plot
p.all_ge_bac_rel_combined = ggplot(df_by_BioProject) + 
  aes(x = Sample_BioProject, y = value, fill = Target_plot) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.85, alpha = 0.95) +
  facet_grid(~Sample, scales = 'free_x', space = 'free_x') + 
  theme_bw(base_size = 12) + 
  ylab('Relative Abundance (%)') + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c('0', '', '50', '', '100'),
                     limits = c(0, 1),
                     oob = scales::squish) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), 
    axis.title.y = element_text(size = 15), 
    axis.text.y = element_text(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    legend.position = 'right', 
    legend.box = "vertical", 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic")   # <--- italic legend text
  ) + 
  scale_fill_manual(name = "TAXA", values = color_pall, drop = FALSE,
                    breaks = levels(df_by_BioProject$Target)) + 
  theme(strip.text.x = element_text(color = "white", face = "bold"), 
        strip.background.x = element_rect(fill = "black")) + 
  guides(fill = guide_legend(reverse = FALSE, ncol = 1)) +
  scale_x_discrete(labels = function(x) gsub('.*_', '', x))

print(p.all_ge_bac_rel_combined)

# Save the ggplot object to a PDF file with the desired aspect ratio
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_phylum_genus_level_relative_abundance_barplot_by_project.pdf", plot = p.all_ge_bac_rel_combined, width = 9, height = 4)


# By Sample

# 4. Summarize for stack bar plot
df_by_BioProject <- df %>%
  dplyr::group_by(Sample, Target) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(scaling_factor = 0.999999 / sum(value, na.rm = TRUE),
                value = value * scaling_factor) %>%
  dplyr::ungroup() %>%
  dplyr::select(-scaling_factor)

# 5. Compute Cutibacterium abundance for ordering
cuti_abundance <- df_by_BioProject %>%
  filter(Target == "    Cutibacterium") %>%
  select(Sample, cuti_value = value)

# 6. Join abundance back and create compound key for unique x-axis
df_by_BioProject <- df_by_BioProject %>%
  left_join(cuti_abundance, by = c("Sample")) %>%
  mutate(cuti_value = ifelse(is.na(cuti_value), 0, cuti_value),
         Sample_BioProject = paste0(Sample))

# 7. Reorder within each Sample group descending by Cutibacterium abundance
df_by_BioProject <- df_by_BioProject %>%
  group_by(Sample) %>%
  mutate(
    Sample_BioProject = forcats::fct_reorder(Sample_BioProject, cuti_value, .desc = TRUE)
  ) %>%
  ungroup()


# 7b. Create new factor just for stacking order (Cutibacterium at bottom)
df_by_BioProject <- df_by_BioProject %>%
  mutate(Target_plot = forcats::fct_relevel(Target, "    Cutibacterium", after = 0))


df_alluvial <- df_by_BioProject

df_alluvial$Target_plot <- factor(df_alluvial$Target_plot, levels = rev(levels(df_alluvial$Target_plot)))

p.alluvial <- ggplot(df_alluvial,
                     aes(x = Sample_BioProject,
                         stratum = Target_plot,
                         alluvium = Target_plot,
                         y = value,
                         fill = Target)) +
  geom_flow(alpha = 0.6, knot.pos = 0.5) +
  geom_stratum(width = 0.6, color = "white") +
  theme_bw(base_size = 12) +
  ylab("Relative Abundance (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, size = 10, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 15),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "italic")
  ) +
  scale_fill_manual(name = "TAXA", values = color_pall, drop = FALSE,
                    breaks = levels(df_by_BioProject$Target)) +  # legend 순서 고정
  theme(strip.text.x = element_text(color = "white", face = "bold"), 
        strip.background.x = element_rect(fill = "black")) + 
  guides(fill = guide_legend(reverse = FALSE, ncol = 1)) +
  scale_x_discrete(labels = function(x) gsub('.*_', '', x))


print(p.alluvial)

# Save the ggplot object to a PDF file with the desired aspect ratio
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_phylum_genus_level_relative_abundance_barplot_by_sample_streamlines.pdf", plot = p.alluvial, width = 6, height = 4)


# 8. Plot, normal bar charts
p.all_ge_bac_rel_combined_sample = ggplot(df_by_BioProject) + 
  aes(x = Sample_BioProject, y = value, fill = Target_plot) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.85, alpha = 0.95) +
  # facet_grid(~Sample, scales = 'free_x', space = 'free_x') + 
  theme_bw(base_size = 12) + 
  ylab('Relative Abundance (%)') + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c('0', '', '50', '', '100'),
                     limits = c(0, 1),
                     oob = scales::squish) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), 
    axis.title.y = element_text(size = 15), 
    axis.text.y = element_text(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    legend.position = 'right', 
    legend.box = "vertical", 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic")   # <--- italic legend text
  ) + 
  scale_fill_manual(name = "TAXA", values = color_pall, drop = FALSE,
                    breaks = levels(df_by_BioProject$Target)) + 
  theme(strip.text.x = element_text(color = "white", face = "bold"), 
        strip.background.x = element_rect(fill = "black")) + 
  guides(fill = guide_legend(reverse = FALSE, ncol = 1)) +
  scale_x_discrete(labels = function(x) gsub('.*_', '', x))

print(p.all_ge_bac_rel_combined_sample)

# Save the ggplot object to a PDF file with the desired aspect ratio
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_phylum_genus_level_relative_abundance_barplot_by_sample.pdf", plot = p.all_ge_bac_rel_combined_sample, width = 5, height = 4)


# As a whole

# Genus

# Helper: Metadata finder
metadata_finder = function(df, metadata) {
  BioProject <- NULL
  Sample <- NULL
  Scalp <- NULL
  
  for (i in as.character(df$variable)) {
    BioProject <- append(BioProject, metadata$BioProject[metadata$BioSample == i])
    Sample <- append(Sample, metadata$Sample[metadata$BioSample == i])
    Scalp <- append(Scalp, metadata$Scalp[metadata$BioSample == i])
  }
  
  df$BioProject <- BioProject
  df$Sample <- Sample
  df$Scalp <- Scalp
  return(df)
}

# Negate %in%
'%ni%' <- Negate('%in%')

# Helper: Custom Target Level assignment (legend order preserved!)
TargetLvl <- function(df, mytab) {
  joinby <- 'BioSample'
  spn <- 3
  
  # Original legend order preserved
  target_taxa <- c(
    'Cutibacterium (A)', 'Corynebacterium (A)', 'Lawsonella (A)',
    'Staphylococcus (B)', 'Streptococcus (B)',
    'Acinetobacter (P)', 'Pseudomonas (P)', 'Enhydrobacter (P)'
  )
  
  df_origin <- df
  df$Species <- gsub('_', ' ', df$Species)
  
  # Filter out non-bacterial kingdoms
  df <- df[!(df$Kingdom %in% c('Archaea', 'Viruses', 'Eukaryota')), ]
  
  # Extract taxonomic columns
  taxa_only <- df[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
  
  # Melt data
  df <- melt(df, id = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  
  # Assign target levels
  df$Target <- df$Kingdom
  df$Target[df$Kingdom == 'Bacteria'] <- 'Others'
  df$Target[df$Phylum == 'Actinomycetota'] <- 'Actinomycetota'
  df$Target[df$Phylum == 'Bacillota'] <- 'Bacillota'
  df$Target[df$Phylum == 'Pseudomonadota'] <- 'Pseudomonadota'
  df$Target[df$Phylum == 'Bacteroidota'] <- 'Bacteroidota'
  df$Target[df$Genus == 'Cutibacterium'] <- 'Cutibacterium (A)'
  df$Target[df$Genus == 'Corynebacterium'] <- 'Corynebacterium (A)'
  df$Target[df$Genus == 'Lawsonella'] <- 'Lawsonella (A)'
  df$Target[df$Genus == 'Staphylococcus'] <- 'Staphylococcus (B)'
  df$Target[df$Genus == 'Streptococcus'] <- 'Streptococcus (B)'
  df$Target[df$Genus == 'Acinetobacter'] <- 'Acinetobacter (P)'
  df$Target[df$Genus == 'Pseudomonas'] <- 'Pseudomonas (P)'
  df$Target[df$Genus == 'Enhydrobacter'] <- 'Enhydrobacter (P)'
  
  # Factor levels (legend order)
  df$Target <- factor(df$Target, levels = target_taxa)
  
  # Remove original taxonomic columns
  df$Kingdom <- NULL
  df$Phylum <- NULL
  df$Class <- NULL
  df$Order <- NULL
  df$Family <- NULL
  df$Genus <- NULL
  df$Species <- NULL
  
  # Summarize by variable and target
  df <- data.frame(df %>% group_by(variable, Target) %>% summarise_all(list(sum)))
  
  # Add metadata
  df <- metadata_finder(df, mytab)
  return(df)
}

# 1. Select, clean, normalize data
merged_16s_species_rel_filtered_subset <- merged_16s_species_rel_filtered %>%
  subset_samples(Dx == "CTRL" & Time_points == "t=1 (Day 0)") %>%
  subset_taxa(Genus != "Methylobacterium") %>%
  subset_taxa(Genus != "Xanthomonas")

ra_rel.bac <- transform_sample_counts(
  merged_16s_species_rel_filtered_subset, 
  function(x) x / sum(x)
)
ra_rel.bac <- prune_taxa(taxa_sums(ra_rel.bac) > 0, ra_rel.bac)

df <- cbind(
  data.frame(tax_table(ra_rel.bac)),
  data.frame(otu_table(ra_rel.bac))
)

my_meta <- sample_data(ra_rel.bac)
my_meta <- as.data.frame(my_meta)
my_meta[,"BioSample"] <- rownames(my_meta)
rownames(my_meta) <- NULL

df_sorted <- df[order(rowSums(df[, 8:ncol(df)]), decreasing = TRUE), ]

df <- TargetLvl(df_sorted, my_meta)

# 2. Summarize for stack bar plot
df_by_BioProject <- df %>%
  dplyr::group_by(Target) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(scaling_factor = 0.999999 / sum(value, na.rm = TRUE),
                value = value * scaling_factor) %>%
  dplyr::select(-scaling_factor)

# 🔹 NA → "Others" 치환
df_by_BioProject$Target <- as.character(df_by_BioProject$Target)
df_by_BioProject$Target[is.na(df_by_BioProject$Target)] <- "Others"

# 🔹 Target factor 레벨 재정의 (Others 맨 마지막에 위치)
target_levels <- c(
  'Cutibacterium (A)','Corynebacterium (A)','Lawsonella (A)',
  'Staphylococcus (B)','Streptococcus (B)','Acinetobacter (P)',
  'Pseudomonas (P)','Enhydrobacter (P)','Others'
)
df_by_BioProject$Target <- factor(df_by_BioProject$Target, levels = target_levels)

# stacking order용 Target_plot 도 같이 정의
df_by_BioProject <- df_by_BioProject %>%
  mutate(Target_plot = Target)

# 3. Plot
# Summarize again per Target (already done above)
df_bar_labels <- df_by_BioProject %>%
  dplyr::group_by(Target) %>%
  dplyr::summarise(total = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Plot
p.all_ge_bac_rel_combined_all = ggplot(df_by_BioProject) + 
  aes(x = Target, y = value, fill = Target_plot) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), 
           width = 0.85, alpha = 0.95) +
  # 🔹 Add one total number above each bar
  geom_text(data = df_bar_labels,
            aes(x = Target, y = total, label = round(total*100, 1)),
            vjust = -0.5, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  theme_bw(base_size = 12) + 
  ylab('Relative Abundance (%)') + 
  scale_y_continuous(
    breaks = c(0, 0.1, 0.2, 0.3, 0.4),
    labels = c('0', '10', '20', '30','40'),
    limits = c(0, 0.45),   # extend for labels
    oob = scales::squish
  ) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), 
    axis.title.y = element_text(size = 15), 
    axis.text.y = element_text(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    legend.position = 'right', 
    legend.box = "vertical", 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic")
  ) + 
  scale_fill_manual(
    name = "Top 8 genera",
    values = c(
      "#1560BD",  # Cutibacterium
      "#1E90FF",  # Corynebacterium
      "#87CEEB",  # Lawsonella
      "#FF0000",  # Staphylococcus
      "#F08080",  # Streptococcus
      "#228B22",  # Acinetobacter
      "#00FF00",  # Pseudomonas
      "#7FFF00",  # Enhydrobacter
      "#D3D3D3"   # Others
    ), 
    drop = FALSE,
    breaks = levels(df_by_BioProject$Target)
  ) + 
  guides(fill = guide_legend(reverse = FALSE, ncol = 1))

print(p.all_ge_bac_rel_combined_all)

# Save the ggplot object to a PDF file with the desired aspect ratio
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_genus_level_relative_abundance_barplot_all.pdf", plot = p.all_ge_bac_rel_combined_all, width = 7, height = 4)

# Phylum

metadata_finder = function(df, metadata){
  BioProject = NULL
  Sample = NULL
  Scalp = NULL
  for (i in as.character(df$variable)){
    BioProject = append(BioProject, metadata$BioProject[metadata$BioSample == i])
    Sample <- append(Sample, metadata$Sample[metadata$BioSample == i])
    Scalp = append(Scalp, metadata$Scalp[metadata$BioSample== i])
  }
  df$BioProject <- BioProject
  df$Sample <- Sample
  df$Scalp <- Scalp
  return(df)
}
'%ni%' <- Negate('%in%')

# Helper: Custom Target Level assignment (legend order preserved!)
TargetLvl <- function(df, mytab){
  joinby <- 'BioSample'
  spn<-3
  # original legend order preserved:
  target_taxa <- c(
    'Actinomycetota', 
    'Bacillota', 
    'Pseudomonadota',
    'Bacteroidota',
    'Others'
  )
  df_origin <- df
  df$Species <- gsub('_', ' ', df$Species)
  df <- df[!(df$Kingdom %in% c('Archaea','Viruses','Eukaryota')),]
  taxa_only <- df[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
  df <- melt(df, id=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  df$Target <- df$Kingdom
  df$Target[df$Kingdom == 'Bacteria'] <- 'Others'
  df$Target[df$Phylum == 'Actinomycetota'] <- 'Actinomycetota'
  df$Target[df$Phylum == 'Bacillota'] <- 'Bacillota'
  df$Target[df$Phylum == 'Pseudomonadota'] <- 'Pseudomonadota'
  df$Target[df$Phylum == 'Bacteroidota'] <- 'Bacteroidota'
  df$Target <- factor(df$Target, levels=target_taxa)
  df$Kingdom <- NULL
  df$Phylum <- NULL
  df$Class <- NULL
  df$Order <- NULL
  df$Family <- NULL
  df$Genus <- NULL
  df$Species <- NULL
  df <- data.frame(df %>% group_by(variable, Target) %>% summarise_all(list(sum)))
  df <- metadata_finder(df, mytab)
  return(df)
}


# 1. Select, clean, normalize data
merged_16s_species_rel_filtered_subset <- merged_16s_species_rel_filtered %>%
  subset_samples(Dx == "CTRL" & Time_points == "t=1 (Day 0)") %>%
  subset_taxa(Genus != "Methylobacterium") %>%
  subset_taxa(Genus != "Xanthomonas")

ra_rel.bac <- transform_sample_counts(merged_16s_species_rel_filtered_subset, function(x) x/sum(x))
ra_rel.bac <- prune_taxa(taxa_sums(ra_rel.bac)>0, ra_rel.bac)
df <- cbind(data.frame(tax_table(ra_rel.bac)), data.frame(otu_table(ra_rel.bac)))
my_meta <- sample_data(ra_rel.bac)
my_meta <- as.data.frame(my_meta)
my_meta[,"BioSample"] <- rownames(my_meta) 
rownames(my_meta) <- NULL 
df_sorted <- df[order(rowSums(df[, 8:ncol(df)]), decreasing=TRUE), ]
df <- TargetLvl(df_sorted, my_meta)

# 2. Summarize for stack bar plot
df_by_BioProject <- df %>%
  dplyr::group_by(Target) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(scaling_factor = 0.999999 / sum(value, na.rm = TRUE),
                value = value * scaling_factor) %>%
  dplyr::select(-scaling_factor)

# 🔹 NA → "Others" 치환
df_by_BioProject$Target <- as.character(df_by_BioProject$Target)
df_by_BioProject$Target[is.na(df_by_BioProject$Target)] <- "Others"

# 🔹 Target factor 레벨 재정의 (Others 맨 마지막에 위치)
target_levels <- df_by_BioProject$Target 
df_by_BioProject$Target <- factor(df_by_BioProject$Target, levels = target_levels)

# stacking order용 Target_plot 도 같이 정의
df_by_BioProject <- df_by_BioProject %>%
  mutate(Target_plot = Target)

# 3. Plot
# Summarize again per Target (already done above)
df_bar_labels <- df_by_BioProject %>%
  dplyr::group_by(Target) %>%
  dplyr::summarise(total = sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Plot
p.all_ph_bac_rel_combined_all = ggplot(df_by_BioProject) + 
  aes(x = Target, y = value, fill = Target_plot) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), 
           width = 0.85, alpha = 0.95) +
  # 🔹 Add one total number above each bar
  geom_text(data = df_bar_labels,
            aes(x = Target, y = total, label = round(total*100, 1)),
            vjust = -0.5, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  theme_bw(base_size = 12) + 
  ylab('Relative Abundance (%)') + 
  scale_y_continuous(
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    labels = c('0', '10', '20', '30','40', '50'),
    limits = c(0, 0.5),   # extend for labels
    oob = scales::squish
  ) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1), 
    axis.title.y = element_text(size = 15), 
    axis.text.y = element_text(),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    legend.position = 'right', 
    legend.box = "vertical", 
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic")
  ) + 
  scale_fill_manual(
    name = "Top 4 phyla",
    values = c(
      "#000080",  # Actinomycetota
      "#8B0000",  # Bacillota
      "#6B8E23",  # Pseudomonadota 
      "#FFD700",  # Bacteroidota
      "#D3D3D3"   # Others
    ), 
    drop = FALSE,
    breaks = levels(df_by_BioProject$Target)
  ) + 
  guides(fill = guide_legend(reverse = FALSE, ncol = 1))

print(p.all_ph_bac_rel_combined_all)


# Save the ggplot object to a PDF file with the desired aspect ratio
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_phyla_level_relative_abundance_barplot_all.pdf", plot = p.all_ph_bac_rel_combined_all, width = 5, height = 4)





