# ----------------------------------------------------------
# 패키지 로드
# ----------------------------------------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(stringr)
library(rstatix)

setwd("/media/jwchoi/ssd2/projects/microbiome/16S_1/R")

# ----------------------------------------------------------
# 1. Subset / Genus level aggregation
# ----------------------------------------------------------
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

ps <- merged_16s_species_rel_filtered_subset

if (!exists("ps_genus")) {
  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
}

sample_data(ps_genus)$Sample <- factor(sample_data(ps_genus)$Sample,
                                       levels = c("Hair", "Swab", "Tissue"))

# ----------------------------------------------------------
# 2. Relative abundance → CLR transform
# ----------------------------------------------------------
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x/sum(x))
ps_genus_rel <- prune_taxa(taxa_sums(ps_genus_rel) > 0, ps_genus_rel)

mat <- as(otu_table(ps_genus_rel), "matrix")
if (!taxa_are_rows(ps_genus_rel)) {
  mat <- t(mat)
}

pseudocount <- 1
log_mat <- log(mat + pseudocount)
clr_mat <- sweep(log_mat, 2, colMeans(log_mat), FUN = "-")

# ----------------------------------------------------------
# 3. 관심 Genus 추출
# ----------------------------------------------------------
target_genera <- c("Cutibacterium","Corynebacterium","Lawsonella",
                   "Staphylococcus")

tax <- as.data.frame(tax_table(ps_genus_rel))
tax$TaxaID <- rownames(tax)
genus_names <- setNames(as.character(tax$Genus), tax$TaxaID)

keep_taxa <- names(genus_names)[genus_names %in% target_genera]
clr_sub <- clr_mat[keep_taxa, , drop = FALSE]

# ----------------------------------------------------------
# 4. 데이터 프레임 변환
# ----------------------------------------------------------
df <- as.data.frame(clr_sub) %>%
  rownames_to_column("TaxaID") %>%
  pivot_longer(-TaxaID, names_to = "SampleID", values_to = "CLR") %>%
  mutate(Genus = genus_names[TaxaID])

meta <- as(sample_data(ps_genus_rel), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, Sample)

df <- df %>% inner_join(meta, by = "SampleID")
df$Genus  <- factor(df$Genus,  levels = target_genera)
df$Sample <- factor(df$Sample, levels = c("Hair","Swab","Tissue"))

# ----------------------------------------------------------
# 5. Outlier trimming (상하위 0.1% 잘라냄)
# ----------------------------------------------------------
df_trim <- df %>%
  group_by(Genus, Sample) %>%
  mutate(
    q_low  = quantile(CLR, 0.001),
    q_high = quantile(CLR, 0.999),
    CLR_trim = pmin(pmax(CLR, q_low), q_high)
  ) %>%
  ungroup()

# ----------------------------------------------------------
# 6. Wilcoxon test + 보정된 p<0.05만 추출
# ----------------------------------------------------------
comparisons <- list(c("Hair","Swab"), c("Hair","Tissue"), c("Swab","Tissue"))

stat_test <- df_trim %>%
  group_by(Genus) %>%
  wilcox_test(CLR_trim ~ Sample, comparisons = comparisons) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# facet별 y.position 고정 (0.9에서 시작, 0.05 간격 아래로)
stat_sig <- stat_test %>%
  filter(p.adj < 0.05) %>%
  group_by(Genus) %>%
  arrange(desc(p.adj)) %>%   # 역순 정렬
  mutate(y.position = 0.8 - (row_number() - 1) * 0.06) %>%
  ungroup()

# ----------------------------------------------------------
# 7. 시각화
# ----------------------------------------------------------
p <- ggplot(df_trim, aes(x = Sample, y = CLR_trim, color = Sample)) +
  geom_boxplot(fill = NA, width = 0.6, linewidth = 1,
               outlier.shape = NA, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 4) +
  stat_pvalue_manual(stat_sig, label = "p.adj.signif", size = 6, tip.length = 0.02) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = NULL, y = "CLR abundance", color = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "black", color = "black"),
        strip.text = element_text(color = "white", face = "bold.italic", size = 12),
        axis.text.x  = element_text(size = 12),   # x축 눈금 크기
        axis.text.y  = element_text(size = 12),
        legend.text  = element_text(size = 13), 
        legend.position = "right") +
  scale_y_continuous(limits = c(-0.025, 0.85),   # 안전하게 0.95까지 확장
                     breaks = seq(0, 0.75, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  guides(color = guide_legend(override.aes = list(size = 4)))

print(p)

ggsave("16s_scalp_final_genus_level_clr_diff_analysis.pdf", plot = p,
       width = 12, height = 4, units = "in", dpi = 600)
