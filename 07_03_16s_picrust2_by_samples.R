# -----------------------------
# 패키지 로드
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(ggpicrust2)



# -----------------------------
# 디렉토리 설정
# -----------------------------
picrust2_out <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/picrust2_input_merged/picrust2_out"
output_dir <- file.path(picrust2_out, "R_figures")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Pathway table 불러오기 (descrip 파일)
# -----------------------------
path_file_desc <- file.path(picrust2_out, "pathways_out/path_abun_unstrat_descrip.tsv")

# 첫 두 컬럼은 id, description
path_annot <- read_tsv(path_file_desc, guess_max = 10000)
colnames(path_annot)[1:2] <- c("id","description")

# 샘플 abundance 행렬
sample_cols <- setdiff(colnames(path_annot), c("id","description"))
path_mat <- path_annot %>%
  dplyr::select(id, all_of(sample_cols)) %>%
  column_to_rownames("id") %>%
  as.matrix()

# 라벨 매핑
desc_map <- setNames(path_annot$description, path_annot$id)
label_map <- ifelse(is.na(desc_map) | desc_map == "", names(desc_map), desc_map)
dup <- duplicated(label_map)
label_map[dup] <- paste0(label_map[dup], " (", names(label_map)[dup], ")")

# -----------------------------
# 상위 20개 Pathway 막대그래프
# -----------------------------
top_n <- 20
top_path <- rowSums(path_mat) %>%
  sort(decreasing = TRUE) %>%
  head(top_n) %>%
  names()

path_top_df <- path_mat[top_path, , drop=FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("Pathway_ID") %>%
  mutate(Pathway_Label = label_map[Pathway_ID])

path_top_long <- path_top_df %>%
  pivot_longer(-c(Pathway_ID, Pathway_Label), names_to="Sample", values_to="Abundance")

pal_path <- colorRampPalette(brewer.pal(12, "Paired"))(top_n)
ggplot(path_top_long, aes(x=Sample, y=Abundance, fill=Pathway_Label)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = pal_path) +
  ggtitle("Top 20 MetaCyc Pathways (with description)")
ggsave(file.path(output_dir, "Top20_Pathway_barplot_named.pdf"), width=12, height=6)

# -----------------------------
# Heatmap
# -----------------------------
pheatmap(path_mat[top_path, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         labels_row = label_map[top_path],
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(50),
         filename = file.path(output_dir, "Top20_Pathway_heatmap_named.pdf"),
         width = 10, height = 8)

# -----------------------------
# PCA
# -----------------------------
path_scaled <- t(scale(t(path_mat)))
pca <- prcomp(t(path_scaled), center=TRUE, scale.=FALSE)
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=colnames(path_mat))
ggplot(pca_df, aes(x=PC1, y=PC2, label=Sample)) +
  geom_point(size=4) +
  geom_text(vjust=-0.5) +
  theme_bw() +
  ggtitle("PCA of MetaCyc Pathway Abundance")
ggsave(file.path(output_dir, "Pathway_PCA.pdf"), width=8, height=6)