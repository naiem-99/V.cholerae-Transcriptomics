# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 5 │ SAMPLE METADATA — PAIRED DESIGN TABLE                         #
# ─────────────────────────────────────────────────────────────────────────── #
Paired_Sample <- data.frame(
  Sample      = c("S01", "S02", "S03", "S04", "S05",
                  "V01", "V02", "V03", "V04", "V05"),
  Sample_type = c("Stool", "Stool", "Stool", "Stool", "Stool",
                  "Vomit", "Vomit", "Vomit", "Vomit", "Vomit"),
  Patient     = factor(c("P01", "P02", "P03", "P04", "P05",
                         "P01", "P02", "P03", "P04", "P05")))
Paired_Sample
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 6 │ DGEList + TMM NORMALISATION                                    #
# ─────────────────────────────────────────────────────────────────────────── #
set.seed(123)
dge       <- DGEList(counts = Merged_Data_paired_count[, -1])
dge$genes <- data.frame(GENEID = Merged_Data_paired_count$GENEID)
dge       <- calcNormFactors(dge, method = "TMM")
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 7 │ GENE FILTERING (CPM >= 1 in at least one sample)              #
# ─────────────────────────────────────────────────────────────────────────── #

genes_keep <- rowMaxs(cpm(dge)) >= 1
d          <- dge[genes_keep, ]
dim(d)
After_filter <-as.data.frame(d$genes)
checkkk_gene <- as.data.frame(dge)%>%anti_join(After_filter, by="GENEID")

write.csv(checkkk_gene,"1.Final_Paired_results/checkkk_gene.csv",row.names = F)
#GTF72_10075
#VC_0846
# ───#──────────────────────────────────────────────────────────────────────── #
#  SECTION 8 │ PAIRED DESIGN MATRIX (~Patient + group)                       #
# ─────────────────────────────────────────────────────────────────────────── #
group   <- factor(Paired_Sample$Sample_type, levels = c("Stool", "Vomit"))
Patient <- Paired_Sample$Patient

design              <- model.matrix(~Patient + group)

rownames(design)    <- Paired_Sample$Sample

design

colnames(design)
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 9 │ VOOM TRANSFORMATION                                            #
# ─────────────────────────────────────────────────────────────────────────── #

v <- voom(d, design, plot = TRUE)
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 10 │ LINEAR MODEL + CONTRAST + EBAYES                             #
# ─────────────────────────────────────────────────────────────────────────── #

fit  <- lmFit(v, design)

contrast_matrix <- makeContrasts(Stool_vs_Vomit = -groupVomit,levels = colnames(design))

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 11 │ DE RESULTS + GENE SYMBOLS                                     #
# ─────────────────────────────────────────────────────────────────────────── #

DE_results <- topTable(fit2, coef = "Stool_vs_Vomit",adjust = "BH", sort.by = "logFC", n = Inf) %>% as.data.frame()
DE_results <- DE_results %>%
  mutate(Expression = case_when(
    logFC >=  0.5 & adj.P.Val < 0.054 ~ "Increased in Stool",
    logFC <= -0.5 & adj.P.Val < 0.054 ~ "Increased in Vomit",
    TRUE ~ "NS"))%>%
       mutate(Expression2 = case_when(
      logFC >=  0.5 & P.Value < 0.05 ~ "Up-regulated in Stool",
      logFC <= -0.5 & P.Value < 0.05 ~ "Up-regulated in Vomit",
      TRUE                            ~ "Stable"))%>%
        left_join(gene_symbols, by = "GENEID") %>%
        as.data.frame() %>%
        mutate(Expression = factor(
        Expression,
      levels = c("Increased in Stool", "NS", "Increased in Vomit")))

DE_results[1:5,]

#-----------------------------------------------------------------------------

library(dplyr)
library(dplyr)
library(readr)

# 1. Load the data using read_csv (better for tidyverse)
Table_S2 <- read_csv("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Revisions_to_RNAseq_V._cholerae_paper/table_s2-annotated-2026.07.21.csv")

# 2. OPTIONAL: Run this if you still get an error to inspect the names
# names(Table_S2)

# 3. Select and rename
Table_S2 <- Table_S2 %>% 
  dplyr::select(`Locus Tag`, `Gene Symbol`, `N16961 reference genome annotation`, 
         `NCBI Protein Search`, `Gene References`, `Pathogenicity Island`, 
         `Results from microarray analysis in LaRocque et al. 2005`) %>%
  dplyr::rename(GENEID = `Locus Tag`)

DE_results <- DE_results %>% inner_join(Table_S2, by="GENEID")%>%
              mutate(Heatmap_Gene_Name = paste0(`Gene Symbol`, " , ",`NCBI Protein Search`))

write.csv(DE_results, "1.Final_Paired_results/Stool_vs_Vomit_paired_DEgs.csv", row.names = FALSE)
#------------------------------------------------------------------------------------------------
