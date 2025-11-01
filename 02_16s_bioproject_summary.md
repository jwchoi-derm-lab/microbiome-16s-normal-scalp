Microbiome 16S Normal Scalp Metadata Summary
================
2025-10-31

``` r
library(dplyr)
library(readr)
library(stringr)
library(moonBook)
library(phyloseq)
library(knitr)

# Load data
meta <- read_csv("data/scalp_16s_o.csv")
merged_16s_species_rel_filtered_subset <- readRDS("data/merged_16s_species_rel_filtered_subset.rds")

filtered <- meta %>%
  filter(Seq == "16S",
         Site == "Scalp",
         Dx == "CTRL",
         Time_points == "t=1 (Day 0)" & !grepl("_u", Subject_id) 
         )

n_bioproject <- filtered %>% distinct(BioProject) %>% nrow()
cat("Number of BioProjects :", n_bioproject, "\n")
```

    ## Number of BioProjects : 17

``` r
cat("Number of samples :", nrow(filtered), "\n")
```

    ## Number of samples : 1190

``` r
filtered <- filtered %>%
  mutate(Subject_core = str_replace(Subject_id, "_.*$", ""))

subject_summary <- filtered %>%
  group_by(BioProject) %>%
  summarise(
    Subject_ids = list(unique(Subject_core)),
    n_unique_subjects = n_distinct(Subject_core),
    .groups = "drop"
  )

cat("Number of subjects:", sum(subject_summary$n_unique_subjects), "\n")
```

    ## Number of subjects: 1112

``` r
# BioProject-wise Metadata Summary
bioproject_summary <- filtered %>%
  group_by(BioProject) %>%
  summarise(
    Meta_sex  = if_else(any(Sex != "NS" & !is.na(Sex)), "Yes", "No"),
    Meta_age  = if_else(any(Age != "NS" & !is.na(Age)), "Yes", "No"),
    Meta_race = if_else(any(Ethnicity != "NS" & !is.na(Ethnicity)), "Yes", "No"),
    Sample_no  = n(),
    Subject_no = n_distinct(Subject_core),
    Layout      = paste(unique(Layout), collapse = ", "),
    Type = paste(sort(unique(Sample)), collapse = ", "),
    File_no = n() * if_else(any(Layout == "PAIRED"), 2, 1),  # ← 여기 수정
    .groups = "drop"
  ) %>%
  ungroup()

cat("Number of sequence files :", sum(bioproject_summary$File_no), "\n\n")
```

    ## Number of sequence files : 2142

``` r
options(width = 2000)
print(bioproject_summary)
```

    ## # A tibble: 17 × 9
    ##    BioProject   Meta_sex Meta_age Meta_race Sample_no Subject_no Layout Type         File_no
    ##    <chr>        <chr>    <chr>    <chr>         <int>      <int> <chr>  <chr>          <dbl>
    ##  1 PRJDB5064    Yes      Yes      Yes              37         37 SINGLE Swab              37
    ##  2 PRJEB16723   Yes      Yes      No               11         11 SINGLE Swab              11
    ##  3 PRJEB25915   No       No       No               25         25 PAIRED Swab              50
    ##  4 PRJEB26870   Yes      No       No               61         61 PAIRED Hair             122
    ##  5 PRJEB62089   Yes      Yes      No               48         48 PAIRED Swab              96
    ##  6 PRJNA1115970 Yes      Yes      No                4          4 SINGLE Swab               4
    ##  7 PRJNA1189034 Yes      No       No               52         52 PAIRED Swab             104
    ##  8 PRJNA1223116 Yes      Yes      Yes             114         57 PAIRED Swab             228
    ##  9 PRJNA1268597 No       No       No               20         20 PAIRED Swab              40
    ## 10 PRJNA314604  Yes      Yes      Yes             110        110 SINGLE Swab             110
    ## 11 PRJNA415710  Yes      Yes      Yes              70         70 PAIRED Swab             140
    ## 12 PRJNA417700  No       No       No               42         21 PAIRED Hair              84
    ## 13 PRJNA510206  No       Yes      No               21         21 PAIRED Swab, Tissue      42
    ## 14 PRJNA542898  No       No       No              494        494 PAIRED Swab             988
    ## 15 PRJNA788988  No       No       No               30         30 SINGLE Swab              30
    ## 16 PRJNA891901  Yes      Yes      Yes              46         46 SINGLE Swab              46
    ## 17 PRJNA953653  Yes      Yes      Yes               5          5 PAIRED Tissue            10

``` r
sample_df <- data.frame(sample_data(merged_16s_species_rel_filtered_subset), 
                        check.names = FALSE, 
                        stringsAsFactors = FALSE)

tbl <- mytable(~ Sex + Ethnicity + Age_group, data = sample_df)
kable(tbl, caption = "Demographic summary of 975 samples used in the final analysis \n\n * NS : Not specified")  
```

Table: Demographic summary of 975 samples used in the final analysis

- NS : Not specified

| name         | stats       | N   | missing | rate    | class       |
|:-------------|:------------|:----|:--------|:--------|:------------|
| Sex          |             | 975 | 0       | ( 0.0%) | categorical |
| \- F         | 225 (23.1%) |     |         |         |             |
| \- M         | 230 (23.6%) |     |         |         |             |
| \- NS        | 520 (53.3%) |     |         |         |             |
| Ethnicity    |             | 975 | 0       | ( 0.0%) | categorical |
| \- African   | 29 (3.0%)   |     |         |         |             |
| \- Asian     | 236 (24.2%) |     |         |         |             |
| \- Caucasian | 16 (1.6%)   |     |         |         |             |
| \- Hispanic  | 20 (2.1%)   |     |         |         |             |
| \- NS        | 674 (69.1%) |     |         |         |             |
| Age_group    |             | 975 | 0       | ( 0.0%) | categorical |
| \- 16-30     | 142 (14.6%) |     |         |         |             |
| \- 31-45     | 115 (11.8%) |     |         |         |             |
| \- 46-60     | 37 (3.8%)   |     |         |         |             |
| \- 61-       | 22 (2.3%)   |     |         |         |             |
| \- NS        | 659 (67.6%) |     |         |         |             |
