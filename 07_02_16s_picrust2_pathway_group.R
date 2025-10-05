# -----------------------------
# 패키지 로드
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(tibble)
library(forcats)

setwd("/media/jwchoi/ssd2/projects/microbiome/16S_1/R")

# -----------------------------
# 디렉토리 설정
# -----------------------------
picrust2_out <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results/picrust2_input_merged/picrust2_out"
output_dir <- file.path(picrust2_out, "R_figures")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Pathway table 불러오기
# -----------------------------
path_file_desc <- file.path(picrust2_out, "pathways_out/path_abun_unstrat_descrip.tsv")
path_annot <- readr::read_tsv(path_file_desc, guess_max = 10000)
colnames(path_annot)[1:2] <- c("id","description")

# -----------------------------
# 샘플 abundance 행렬
# -----------------------------
sample_cols <- setdiff(colnames(path_annot), c("id","description"))
path_mat <- path_annot %>%
  dplyr::select(id, all_of(sample_cols)) %>%
  tibble::column_to_rownames("id") %>%
  as.matrix()

# -----------------------------
# 메타데이터 불러오기
# -----------------------------
meta <- readr::read_csv("/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/scalp_16s.csv") %>%
  dplyr::select(SRR, Sample)

# -----------------------------
# 공통 샘플만 선택
# -----------------------------
common_srr <- intersect(colnames(path_mat), meta$SRR)
path_mat <- path_mat[, common_srr, drop=FALSE]
group_info <- meta %>% filter(SRR %in% common_srr)

# -----------------------------
# Functional category mapping
# -----------------------------
pathway_category <- sapply(path_annot$description, function(desc){
  desc <- tolower(desc)
  if(grepl("glycolysis|gluconeogenesis|pyruvate|pentose|tca|starch|sugar", desc)) return("Carbohydrate metabolism")
  if(grepl("arginine|lysine|methionine|tryptophan|amino acid|valine|leucine|isoleucine", desc)) return("Amino acid metabolism")
  if(grepl("respiration|electron transport|oxidative|atp|energy", desc)) return("Energy metabolism")
  if(grepl("nucleotide|purine|pyrimidine|rna|dna", desc)) return("Nucleotide metabolism")
  if(grepl("cofactor|vitamin|folate|riboflavin|coenzyme|biotin|pantothenate", desc)) return("Cofactors & Vitamins")
  if(grepl("lipid|fatty acid|glycerolipid|sterol|cholesterol", desc)) return("Lipid metabolism")
  if(grepl("glycan|glycosyl|cell wall|peptidoglycan", desc)) return("Glycan biosynthesis & degradation")
  if(grepl("xenobiotic|degradation|environmental", desc)) return("Xenobiotics & Environmental processing")
  if(grepl("signal|quorum sensing|autoinducer", desc)) return("Signal molecules / Quorum sensing")
  if(grepl("vitamin d|cholecalciferol|ergocalciferol", desc)) return("Vitamin D biosynthesis")
  if(grepl("steroid|androgen", desc)) return("Hair / Follicle-related metabolism")
  if(grepl("folate|pantothenate|biotin", desc)) return("Hair-related vitamins")
  if(grepl("growth factor|fgf|igf|tgf|bmp|wnt", desc)) return("Hair follicle signaling")
  return("Others")
})
names(pathway_category) <- rownames(path_mat)

# -----------------------------
# Functional category 합산 및 상대 abundance 계산
# -----------------------------
func_mat <- rowsum(path_mat, group = pathway_category)
func_rel <- sweep(func_mat, 2, colSums(func_mat), FUN="/")

func_rel_df <- as.data.frame(func_rel)
func_rel_df$FunctionalCategory <- rownames(func_rel_df)

group_means <- func_rel_df %>%
  tidyr::pivot_longer(-FunctionalCategory, names_to="SRR", values_to="RelAbundance") %>%
  dplyr::left_join(group_info, by="SRR") %>%
  dplyr::group_by(Sample, FunctionalCategory) %>%
  dplyr::summarise(RelAbundance = mean(RelAbundance, na.rm=TRUE), .groups="drop")

# -----------------------------
# Factor & 색상 설정
# -----------------------------
cats <- unique(group_means$FunctionalCategory)
cats <- c(setdiff(cats, "Others"), "Others")   # Legend 순서용
n_cat <- length(cats)

# 원색 계열 팔레트
pal_func <- c(
  "#66C2A5",  # teal green
  "#FC8D62",  # soft orange
  "#8DA0CB",  # muted blue-purple
  "#E78AC3",  # muted pink
  "#A6D854",  # lime green
  "#FFD92F",  # golden yellow
  "#E5C494",  # tan/brown
  "skyblue",  # deeper teal
  "gray80"    # 밝은 회색, Others
)

pal_func <- pal_func[1:n_cat]
names(pal_func) <- cats

group_means$FunctionalCategory <- factor(group_means$FunctionalCategory, levels = rev(cats))

# -----------------------------
# Alluvial-style Stacked Bar Plot
# -----------------------------
df_labels <- df_alluvial %>%
  group_by(Sample) %>%
  arrange(Sample, desc(Target_plot)) %>%
  mutate(ymax = cumsum(value),
         ymin = ymax - value,
         y_center = (ymin + ymax) / 2) %>%
  ungroup() %>%
  filter(value >= 0.03)   # 3% 이상만 표시

# ----------------------------------------------------------
# 3. Alluvial Plot 생성
# ----------------------------------------------------------
p.alluvial.lab <- ggplot(df_alluvial,
                         aes(x = Sample,
                             stratum = Target_plot,
                             alluvium = Target_plot,
                             y = value,
                             fill = Target_plot)) +
  geom_flow(alpha = 0.7, knot.pos = 0.5) +   # 흐름
  geom_stratum(width = 0.6, color = "white") +  # 각 stratum block
  geom_text(data = df_labels,                 # 라벨 중앙 정렬
            aes(x = Sample, y = y_center,
                label = scales::percent(value, accuracy = 0.1)),
            color = "black", size = 3) +
  scale_fill_manual(values = pal_func, name="Functional category", breaks=cats) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, size = 12, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 14),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "italic")
  ) +
  ylab("Relative Abundance (%)")

print(p.alluvial.lab)

# ----------------------------------------------------------
# 4. PDF 저장
# ----------------------------------------------------------
ggsave(file.path("/media/jwchoi/ssd2/projects/microbiome/16S_1/R/16s_scalp_final_functional_barplot_by_sample_streamlines.pdf"),
       plot = p.alluvial.lab, width=6, height=4)