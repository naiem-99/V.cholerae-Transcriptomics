
# ---------------------------Differential Expression Analysis(limma) -----------------------------------------
#Step P02 : Preparing the paired sample list column
Paired_Sample <- data.frame(
  Sample_type = c("Stool", "Stool", "Stool", "Stool", "Stool",
                  "Vomit", "Vomit", "Vomit", "Vomit", "Vomit"),
  Sample = c("S21", "S23", "S24", "S25", "S26",
             "S30", "S33", "S34", "S35", "S37"))

Paired_Sample

set.seed(123) 

# Check if the data loaded correctly
head(Merged_Data_paired_count)
head(Paired_Sample)
#----------- Create DEGlist obj ----------------
count<- Merged_Data_paired_count
all(colnames(Merged_Data_paired_count)[-1] == Paired_Sample$Sample)
dge <- DGEList(counts = count)
dge
# Pre-processing (calculate normalization factors) 

dge <- calcNormFactors(dge,method = "TMM")
#https://f1000research.com/articles/5-1408
# ------- Filter out low-expressed genes --------
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
d <- dge[-drop,] 
dim(d)
#-----------------------------------------------
#--------------Make design matrix--------------------------------------------------------------
##https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
group <- factor(Paired_Sample$Sample_type, 
                levels = c("Stool","Vomit"))
design <- model.matrix(~0 + group)
design
colnames(design) <- levels(group)

#--------Transform counts using voom for linear modeling--------
v <- voom(d, design, plot = T)
v
#----------- Fiitting the count matrix to linear model-----------
fit <- lmFit(v, design)
# ------Define contrasts for each comparison--------------------
contrast_matrix <- makeContrasts( Stool_vs_Vomit= Stool - Vomit,levels = colnames(design))
contrast_matrix

# ------------Apply contrasts to the model----------------------
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 

# -----------Compute empirical Bayes statistics-----------------
fit2 <- eBayes(fit2)
# -----------Extract results for each contrast------------------
Stool_vs_Vomit_paired  <- topTable(fit2, coef="Stool_vs_Vomit", adjust="BH",sort.by = "logFC",n = Inf)

# -----------Create a new column indicating the expreesion type ------------------
Stool_vs_Vomit_paired  <-Stool_vs_Vomit_paired%>% mutate(Expression = 
                                         case_when(logFC >=0.5 & P.Value < 0.05 ~ "Up-regulated in Stool",
                                                   logFC <= -0.5 & P.Value < 0.05 ~ "Up-regulated in Vomit",TRUE ~ "Stable"))
dim(Stool_vs_Vomit_paired)
table(Stool_vs_Vomit_paired$Expression)
Stool_vs_Vomit_paired_0  <-Stool_vs_Vomit_paired%>% mutate(Expression = 
                                                           case_when(logFC > 0 & P.Value < 0.05 ~ "Up-regulated in Stool",
                                                                     logFC < 0 & P.Value < 0.05 ~ "Up-regulated in Vomit",
                                                                     TRUE ~ "Stable"))
table(Stool_vs_Vomit_paired_0$Expression)

# -----------  Save the result --------------------------------

write.csv(Stool_vs_Vomit_paired, file="1.Final_Paired_results/Stool_vs_Vomit_paired_DEgs.csv",row.names=F)


#---------------------------------------------------------------#

#***********************************----Volcano plot-------******************
#load the required libraries
library(ggplot2)
library(ggrepel)

#-------------------- Set the color for each DEGs------------------------------
cols_N2 <- c("Up-regulated in Stool" = "#FF1493", "Up-regulated in Vomit" = "#26b3ff", "Stable" = "grey50") 

sizes_N2 <- c("Up-regulated in Stool" = 3, "Up-regulated in Vomit" = 3, "Stable" = 3) 
alphas_N2 <- c("Up-regulated in Stool" = 1, "Up-regulated in Vomit" = 1, "Stable" = 0.47)
#-------------------------------------------------------------------------------

#------------------------Select top 5 DEGS for each group-----------------------
up2 = Stool_vs_Vomit_paired %>% filter(Expression == "Up-regulated in Stool") %>%arrange(P.Value,desc(abs(logFC))) %>%head(5)

down2 = Stool_vs_Vomit_paired %>%filter(Expression == "Up-regulated in Vomit") %>%arrange(P.Value,desc(abs(logFC))) %>%head(5)
top_5_genes_pvalue_N2 = bind_rows(up2,down2)
top_5_genes_pvalue_N2
#---------------------------------------------------------------------------------

#------------------Plot the Result------------------------------------------------

Stool_vs_Vomit_paired%>%ggplot(aes(x = logFC,
                            y = -log(P.Value, 10),
                            fill =   Expression,    
                            size =   Expression,
                            alpha =   Expression)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    colour = "black") +geom_hline(yintercept = -log10(0.05),size=0.4,linetype = "dashed") + 
  #geom_vline(xintercept = c(-1,1),linetype = "dashed") +
  scale_fill_manual(values = cols_N2) + 
  scale_size_manual(values = sizes_N2) + 
  scale_alpha_manual(values = alphas_N2) + 
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10))+ theme_bw()+theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.y  = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold") ,  # expression title size
    legend.position = "bottom",legend.text = element_text(size = 16))+xlim(-5.8,5.9)+
  geom_text_repel(data = top_5_genes_pvalue_N2, mapping = aes(logFC, -log(P.Value,10), label = GENEID), check_overlap = FALSE,size = 4.7)

Stool_vs_Vomit_paired
ggsave(filename = "1.Final_Paired_results/Stool_vs_Vomit_Paired_Volcano.png", width = 8.2, height = 7.5, units = "in",dpi= 1200) 



#---------------------------Heatmap for each Group ( Top 10 up regulated ; Functional Name)------------
#Load libraries
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

#------------Select the top 10 up and downregulated genes and filter them for Function gene name to plot-------
Merged_Data
Merged_Data_paired_annotated <- Merged_Data %>%dplyr::select(GENEID,S21,S23,S24,S25,S26,S30,S33,S34,S35,S37,`N16961 annotation`) %>%distinct(GENEID, .keep_all = TRUE)

#---------------Filter top 10 genes from DEGs result-----------------
top10_logFC <- bind_rows(Stool_vs_Vomit_paired %>%  filter(Expression == "Up-regulated in Stool") %>%  arrange(desc(logFC)) %>%  slice(1:10),
  Stool_vs_Vomit_paired %>%filter(Expression == "Up-regulated in Vomit") %>%arrange(logFC) %>%slice(1:10))
head(Merged_Data_paired_annotated)
head(top10_logFC)
#----------------------- Prepare Heatmap data-------------------------
Merged_Data_paired_annotated_top10 <- Merged_Data_paired_annotated %>%filter(GENEID %in% top10_logFC$GENEID) %>%left_join(  top10_logFC %>% select(GENEID, Expression),by = "GENEID")

Merged_Data_paired_annotated_top10 <- Merged_Data_paired_annotated_top10 %>%mutate(Function_short = sub(",.*$", "", `N16961 annotation`))

expr_mat <- Merged_Data_paired_annotated_top10 %>%column_to_rownames("GENEID") %>%
  select(S21, S23, S24, S25, S26, S30, S33, S34, S35, S37) %>%as.matrix()

expr_mat_scaled <- t(scale(t(log2(expr_mat + 1))))

#row_labels <- Merged_Data_paired_annotated_top10$`N16961 annotation`
row_ha <- rowAnnotation(Function = anno_text( Merged_Data_paired_annotated_top10$Function_short,
    gp = gpar(fontsize = 8),just = "left"),annotation_name_gp = gpar(fontsize = 10),width = unit(6, "cm"))

sample_group <- data.frame(Group = c(rep("Stool", 5), rep("Vomit", 5)))
rownames(sample_group) <- colnames(expr_mat_scaled)
ha_top <- HeatmapAnnotation( Group = sample_group$Group, col = list(Group = c("Stool" = "#EE1289", "Vomit" = "green")),annotation_legend_param = list(direction = "horizontal", nrow = 1))

Merged_Data_paired_annotated_top10$Expression <- factor(
  Merged_Data_paired_annotated_top10$Expression, levels = c("Up-regulated in Vomit", "Up-regulated in Stool"))

#---------------Plot the Data----------------

ht <- Heatmap(expr_mat_scaled,
  name = "Scaled log2(x+1)",
  top_annotation = ha_top,
  right_annotation = row_ha,
  row_split = Merged_Data_paired_annotated_top10$Expression,
  show_row_names = F,
  show_column_names = T,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,        # HIDE row dendrogram
  show_column_dend = FALSE,     # HIDE column dendrogram
 col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
 heatmap_legend_param = list(direction = "horizontal", nrow = 1))

ht
ht2<-draw(ht,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE)
ht2

#_-------------Save the result---------------------------
png("1.Final_Paired_results/Stool_vs_Vomit_Paired_Heatmap.png", width = 7,height = 7,units = "in", res = 400)

draw(ht,heatmap_legend_side = "bottom",annotation_legend_side = "bottom", merge_legends = TRUE)

dev.off()
#-----------------------------------------The End Differential expression analysis----------------------------

