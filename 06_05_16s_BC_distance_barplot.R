# ----------------------------------------------------------
# 패키지
# ----------------------------------------------------------
library(phyloseq)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
library(readr)

setwd("/media/jwchoi/ssd2/projects/microbiome/16S_1/R")

merged_16s_species_rel_filtered_subset <- merged_16s_species_rel_filtered %>%
  subset_samples(Dx == "CTRL" & Time_points == "t=1 (Day 0)") %>%
  subset_taxa(Genus != "Methylobacterium") %>%
  subset_taxa(Genus != "Xanthomonas")

merged_16s_species_rel_filtered_subset_rel <- transform_sample_counts(
  merged_16s_species_rel_filtered_subset, function(x) x/sum(x)
)

# ----------------------------------------------------------
# 2. Bray-Curtis 거리행렬
# ----------------------------------------------------------
pid_bray <- phyloseq::distance(merged_16s_species_rel_filtered_subset_rel, method = "bray")
pid_bray <- as.matrix(pid_bray)

# ----------------------------------------------------------
# 3. 그룹별 distance 정리 (Sample factor: Swab, Hair, Tissue)
# ----------------------------------------------------------
groups_all <- sample_data(merged_16s_species_rel_filtered_subset_rel)$Sample
sub_dist <- list()

for (group in groups_all) {
  row_group <- which(groups_all == group)
  sample_group <- sample_names(merged_16s_species_rel_filtered_subset_rel)[row_group]
  sub_dist[[group]] <- pid_bray[sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

# 데이터프레임 변환
braygroups <- melt(sub_dist)
df_bray_combined <- braygroups[complete.cases(braygroups), ]
df_bray_combined$L1 <- factor(df_bray_combined$L1,
                              levels = c("Hair", "Swab", "Tissue"))
rownames(df_bray_combined) <- NULL

# ----------------------------------------------------------
# 4. Tukey HSD (Sex 무시, 전체 그룹 간 비교)
# ----------------------------------------------------------
tukey_all <- df_bray_combined %>%
  tukey_hsd(value ~ L1) %>%
  add_significance() %>%
  add_xy_position(x = "L1") %>%
  filter(p.adj < 0.05) 

# ----------------------------------------------------------
# 5. 박스플롯 + 브라켓
# ----------------------------------------------------------
bray_pid_all <- ggboxplot(df_bray_combined, x = "L1", y = "value",
                          color = "black",
                          fill = "L1",
                          outlier.shape = NA) +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label = "Adjusted p < 0.001"),
    data = tukey_all, 
    y.position = max(df_bray_combined$value, na.rm = TRUE) * 1.05,  # ← 더 위로
    step.increase = 0.2,
    tip.length = 0.01,
    vjust = -1
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(vjust=0.5, angle = 30, size = 10),
    panel.grid = element_blank()   # ← 그리드 제거
  ) +
  ylab("Bray-Curtis dissimilarity") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.0, by = 0.5)) +
  scale_x_discrete(limits = c("Hair", "Swab", "Tissue")) +
  scale_fill_manual(
    limits = c("Hair", "Swab", "Tissue"), 
    values = c("#D5AFF5", "#8AB59C", "#65B8F8")
  )


# ----------------------------------------------------------
# 6. 결과 출력 및 저장
# ----------------------------------------------------------
print(bray_pid_all)

write_csv(tukey_all, "16s_scalp_final_beta_bray_curtis_statistics_all.csv")
ggsave(filename = "/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_beta_bray_curtis_boxplot.pdf",
       plot = bray_pid_all, width = 4, height = 4)
