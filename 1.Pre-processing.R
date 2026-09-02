# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 1 │ WORKSPACE & PACKAGES                                           #
# ─────────────────────────────────────────────────────────────────────────── #

setwd("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired")

library(limma);        library(edgeR)
library(dplyr);        library(tibble)
library(stringr);      library(matrixStats)
library(ggplot2);      library(ggrepel)
library(patchwork);    library(scales)
library(ComplexHeatmap); library(circlize); library(grid)
library(clusterProfiler); library(enrichplot)
library(rtracklayer)
library(vegan)

rm(list = ls()); gc()

dir.create("1.Final_Paired_results", showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 2 │ DATA LOADING                                                   #
# ─────────────────────────────────────────────────────────────────────────── #
#---------------------- Read the count data ---------------------------------------------
name_map <- c("S21" = "S01", "S23" = "S02", "S24" = "S03", "S25" = "S04", "S26" = "S05",
              "S30" = "V01", "S33" = "V02", "S34" = "V03", "S35" = "V04", "S37" = "V05")

count_data_all <- read.csv("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Raw_Files/ASVS.counts.genes.csv", check.names = FALSE) %>%
  rename_with(~ ifelse(.x %in% names(name_map), name_map[.x], .x))


#--------------------- Read the annotation file -------------------------------------------
annotation_file <- read.csv("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Raw_Files/VC_genes-Ng-file.csv", check.names = FALSE)%>% 
                   distinct(`N16961 locus`, .keep_all = T)


#-------------------Read the Uniprot proteome file------------------------------------------------------------
uniprotkb_proteome_file <- read.csv(
  "D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Raw_Files/uniprotkb_proteome_UP000000584_2026_06_17.csv", check.names = FALSE) %>%
  dplyr::select(`Protein names`, `Gene Names`) %>%
  mutate(
    GENEID    = ifelse(is.na(str_extract(`Gene Names`, "VC_[A-Z0-9]+")), "0",
                       str_extract(`Gene Names`, "VC_[A-Z0-9]+")),
    Gene_Name = ifelse(is.na(na_if(trimws(gsub("VC_[A-Z0-9]+", "", `Gene Names`)), "")), "0",
                       na_if(trimws(gsub("VC_[A-Z0-9]+", "", `Gene Names`)), ""))
  ) %>%dplyr::select(GENEID, Gene_Name)

uniprotkb_proteome_file[1:5,]

#-------------------------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 3 │ MERGE & ANNOTATE                                               #
# ─────────────────────────────────────────────────────────────────────────── #
Merged_Data2 <- count_data_all %>%
  inner_join(annotation_file, by = "Geneid") %>%
  mutate(GENEID = `N16961 locus`) %>%as.data.frame() %>%
  left_join(uniprotkb_proteome_file, by = "GENEID") %>%
  mutate(GeneName_Merged = case_when(
  !is.na(`C6706 symbol`) & `C6706 symbol` != "0" ~ `C6706 symbol`,
      TRUE ~ Gene_Name))

Merged_Data2[1:5,]

write.csv(Merged_Data2, "D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Merged_Data2.csv", row.names = FALSE)

#----------------------------------------------------------------------------
# Gene symbol lookup
gene_symbols <- Merged_Data2 %>%
  distinct(GENEID, .keep_all = TRUE) %>%
  mutate(
    Symbol2     = case_when(
      !is.na(GeneName_Merged) &
        GeneName_Merged != "0" &
        GeneName_Merged != ""   ~ GeneName_Merged,
      TRUE                      ~ GENEID
    ),
    Description = `N16961 annotation`,
    Symbol      = ifelse(Symbol2 == GENEID, GENEID, paste0(GENEID, " = ", Symbol2))
  ) %>%
  dplyr::select(GENEID, Symbol2, Symbol, Description,`C6706 annotation`)

gene_symbols[1:5,]

# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 4 │ PAIRED COUNT MATRIX                                            #
# ─────────────────────────────────────────────────────────────────────────── #

Merged_Data_paired_count <- Merged_Data2 %>%dplyr::select(GENEID,S01, S02, S03, S04, S05,V01, V02, V03, V04, V05) %>% distinct(GENEID, .keep_all = TRUE)

dim(Merged_Data_paired_count)

Merged_Data_paired_count[1:5,]

stopifnot("Column order mismatch!" = all(colnames(Merged_Data_paired_count)[-1] == Paired_Sample$Sample))

write.csv(Merged_Data_paired_count,"1.Final_Paired_results/Merged_Data_paired_count.csv", row.names = FALSE)
#-----------------------------------------------------------------------------------------------------------------#
