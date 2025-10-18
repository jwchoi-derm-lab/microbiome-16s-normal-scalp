Microbiome 16S Normal Scalp Metadata Summary
================
2025-10-22

``` r
library(dplyr)
library(readr)
library(stringr)

# Load Metadata
meta <- read_csv("metadata/scalp_16s_o.csv")

filtered <- meta %>%
filter(Seq == "16S",
Site == "Scalp",
Dx == "CTRL",
Time_points == "t=1 (Day 0)")

n_bioproject <- filtered %>% distinct(BioProject) %>% nrow()
cat("Number of BioProjects :", n_bioproject, "\n")
```

    ## Number of BioProjects : 17

``` r
cat("Number of samples :", nrow(filtered), "\n")
```

    ## Number of samples : 1239

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

    ## Number of subjects: 1114

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
.groups = "drop"
) %>%
ungroup()

options(width = 2000)
print(bioproject_summary)
```

    ## # A tibble: 17 × 8
    ##    BioProject   Meta_sex Meta_age Meta_race Sample_no Subject_no Layout Type        
    ##    <chr>        <chr>    <chr>    <chr>         <int>      <int> <chr>  <chr>       
    ##  1 PRJDB5064    Yes      Yes      Yes              37         37 SINGLE Swab        
    ##  2 PRJEB16723   Yes      Yes      No               11         11 SINGLE Swab        
    ##  3 PRJEB25915   No       No       No               25         25 PAIRED Swab        
    ##  4 PRJEB26870   Yes      No       No               61         61 PAIRED Hair        
    ##  5 PRJEB62089   Yes      Yes      No               48         48 PAIRED Swab        
    ##  6 PRJNA1115970 Yes      Yes      No                4          4 SINGLE Swab        
    ##  7 PRJNA1189034 Yes      No       No              101         54 PAIRED Swab        
    ##  8 PRJNA1223116 Yes      Yes      Yes             114         57 PAIRED Swab        
    ##  9 PRJNA1268597 No       No       No               20         20 PAIRED Swab        
    ## 10 PRJNA314604  Yes      Yes      Yes             110        110 SINGLE Swab        
    ## 11 PRJNA415710  Yes      Yes      Yes              70         70 PAIRED Swab        
    ## 12 PRJNA417700  No       No       No               42         21 PAIRED Hair        
    ## 13 PRJNA510206  No       Yes      No               21         21 PAIRED Swab, Tissue
    ## 14 PRJNA542898  No       No       No              494        494 PAIRED Swab        
    ## 15 PRJNA788988  No       No       No               30         30 SINGLE Swab        
    ## 16 PRJNA891901  Yes      Yes      Yes              46         46 SINGLE Swab        
    ## 17 PRJNA953653  Yes      Yes      Yes               5          5 PAIRED Tissue
