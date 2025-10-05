# -----------------------------
# 라이브러리
# -----------------------------
library(phyloseq)
library(ggplot2)
library(ggtext)
library(ggExtra)
library(dplyr)
library(glue)

# -----------------------------
# 데이터 subset
# -----------------------------
setwd("/media/jwchoi/ssd2/projects/microbiome/16S_1/R")

merged_16s_species_rel_filtered_subset <- merged_16s_species_rel_filtered %>%
  subset_samples(Dx == "CTRL" & Time_points == "t=1 (Day 0)") %>%
  subset_taxa(Genus != "Methylobacterium") %>%
  subset_taxa(Genus != "Xanthomonas")

merged_16s_species_rel_filtered_subset_rel <- transform_sample_counts(
  merged_16s_species_rel_filtered_subset, function(x) x/sum(x)
)

# -----------------------------
# PERMANOVA
# -----------------------------
pid_div_rel_pcoa <- merged_16s_species_rel_filtered_subset_rel
index = "bray"
seed = 42
type = "Sample"
plot = "PCoA"

set.seed(seed)
x.dist <- phyloseq::distance(pid_div_rel_pcoa, method = index)
ord <- phyloseq::ordinate(pid_div_rel_pcoa, plot, index)

myfunc <- function(v1) deparse(substitute(v1))
dist <- myfunc(x.dist)
a <- adonis2(as.formula(glue("{dist} ~ {type}")),
             data=data.frame(sample_data(pid_div_rel_pcoa)),
             permutations=9999, method=index)
Perm.p <- a$`Pr(>F)`[1]
Perm.p <- ifelse(Perm.p < 0.001, "<0.001", Perm.p)

# -----------------------------
# PCoA 데이터 변환
# -----------------------------
mat <- ord$vectors[, 1:2] %>% as.data.frame() %>%
  mutate(SampleID = rownames(.)) %>%
  arrange(SampleID)

meta <- phyloseq::sample_data(pid_div_rel_pcoa) %>% data.frame() %>%
  mutate(SampleID = rownames(.)) %>%
  arrange(SampleID)

pid_div_rel_pcoa_df <- inner_join(meta, mat, by = "SampleID")

PC1 <- round(ord$values["Relative_eig"][1,]*100, 1)
PC2 <- round(ord$values["Relative_eig"][2,]*100, 1)

# -----------------------------
# 색상/순서 세팅
# -----------------------------
# geom_point 그려지는 순서: Swab → Hair → Tissue
# factor 레벨 설정
pid_div_rel_pcoa_df$Sample <- factor(pid_div_rel_pcoa_df$Sample,
                                     levels = c("Swab", "Hair", "Tissue"))

# 그려지는 순서를 강제하기 위해 데이터 정렬
pid_div_rel_pcoa_df <- pid_div_rel_pcoa_df %>%
  arrange(Sample)   # Swab → Hair → Tissue 순서로 정렬됨

# -----------------------------
# 메인 플롯
# -----------------------------
main.plot <- pid_div_rel_pcoa_df %>% 
  ggplot(aes(x = Axis.1, y = Axis.2)) + 
  geom_vline(xintercept = 0, colour = "grey80") + 
  geom_hline(yintercept = 0, colour = "grey80") + 
  geom_point(aes(color = Sample), alpha = 0.5, size = 2) +
  lims(x = c(-1, 1), y = c(-1, 1)) + 
  stat_ellipse(aes(color = Sample), alpha = 0.8) +
  theme_test() +
  scale_color_manual(
    values = c("Hair"="red", "Swab"="green", "Tissue"="blue"),
    breaks = c("Hair", "Swab", "Tissue")  # 레전드 순서 고정
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +  # 레전드 점 크기
  labs(
    y = paste0("PCoA2 (", PC2, "%)"),
    x = paste0("PCoA1 (", PC1, "%)")
  ) +
  theme(
    plot.caption = element_markdown(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12), 
    legend.position = "bottom",
    plot.margin = unit(c(0,0,0,0),"points"),
    plot.title = element_text(vjust=0, size=12,
                              margin=margin(t=10,b=10))
  ) +
  ggtext::geom_richtext(label.color = NA, size = 3, fill = NA, 
                        hjust = 0, vjust =0, x = -Inf, y = -Inf, 
                        label= paste0("**PERMANOVA** *p*-value ", Perm.p))

# -----------------------------
# ggMarginal 추가 (alpha 적용)
# -----------------------------
final_beta_Plot <- ggExtra::ggMarginal(
  main.plot,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE,
  alpha = 0.5
)

# -----------------------------
# 출력
# -----------------------------
final_beta_Plot

ggsave("16s_scalp_final_beta_pcoa.pdf",
       plot = final_beta_Plot,
       width = 4,
       height = 4,
       units = "in")
