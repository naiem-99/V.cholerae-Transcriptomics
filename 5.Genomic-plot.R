#*************************************************Genomic plot *********************************************
#------------load the required library --------------------------------
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggrepel)
Stool_vs_Vomit_paired<-read.csv("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/1.Final_Paired_results/Stool_vs_Vomit_paired_DEgs.csv",check.names = F)#%>%dplyr::rename(GeneID = GENEID)

Stool_vs_Vomit_paired
# Load GTF file
gtf <- rtracklayer::import("D:/ICDDRB_Feb/1.RNAseq/GCA_013085075.1_ASM1308507v1_genomic.gtf")
gtf
# Keep only genes
genes_gtf <- gtf[gtf$type == "gene"]
genes_gtf
# Build dataframe
gene_annot_01 <- data.frame( Geneid= mcols(genes_gtf)$gene_id,
  Chr = as.character(seqnames(genes_gtf)),Start = start(genes_gtf),
  Stop = end(genes_gtf),stringsAsFactors = FALSE)

gene_annot_01

dim(gene_annot_01)
write.csv(gene_annot_01,"gene_annot_01.csv")

Merged_Data_paired_count_all <- Merged_Data %>%dplyr::select(Geneid,GENEID,S21,S23,S24,S25,S26,S30,S33,S34,S35,S37) %>%
  distinct(GENEID, .keep_all = TRUE)
dim(Merged_Data_paired_count_all)

merged_genenno<- gene_annot_01 %>%inner_join(Merged_Data_paired_count_all,by="Geneid")
dim(merged_genenno)


# -------------------Add midpoint-----------------------------------------
gene_annot <- merged_genenno %>%mutate(Mid_gene = (Start + Stop) / 2)
table(gene_annot$Chr)

# -----------------Merge DE results with gene coordinates--------------------
res_df_all <- Stool_vs_Vomit_paired %>%left_join(gene_annot, by = "GENEID")
dim(res_df_all)
res_df_all_1 <- Stool_vs_Vomit_paired %>%left_join(gene_annot, by = "GENEID")
head(res_df_all)


# ------------3. Prepare chromosome offsets---------------------------
#-------------------------------
chr_lengths <- gene_annot %>%group_by(Chr) %>%
  summarise(Chr_len = max(Stop)) %>% arrange(Chr)

chr_offsets <- chr_lengths %>%mutate(offset = lag(cumsum(Chr_len), default = 0))
#if you want offset as not initial
#res_df_all <- res_df_all %>%left_join(chr_offsets %>% select(Chr, offset), by = "Chr") %%>%mutate(Mid_gene = Mid_gene + offset)

res_df_all <- res_df_all %>%left_join(chr_offsets %>% select(Chr, offset), by = "Chr") #%>%mutate(Mid_gene = Mid_gene + offset)
dev.off()
##################################################################

#********************************************************

# Identify top 10 up and down genes
top15_up <- res_df_all %>% filter(Expression == "Up-regulated in Stool") %>% arrange(desc(logFC)) %>% slice_head(n = 15)

top15_down <- res_df_all %>%filter(Expression == "Up-regulated in Vomit") %>%arrange(logFC) %>%slice_head(n = 15)
top30 <- bind_rows(top15_up, top15_down)
genome_position <- ggplot(res_df_all, aes(x = Mid_gene, y = logFC, color = Expression)) +geom_point(alpha = 0.6, size = 1.2) +
  
  # highlight top genes with blue/red fill
  geom_point(
    data = top30,aes(fill = Expression),shape = 21,color = "black",
    size = 3,stroke = 0.6
  ) +
  
  ggrepel::geom_text_repel(data = top30,aes(label = GENEID),size = 3,color = "black",max.overlaps = Inf) +
  
  scale_color_manual(values = c(
    "Stable" = "grey",
    "Up-regulated in Stool" = "red",
    "Up-regulated in Vomit" = "blue"
  )) +
  
  scale_fill_manual(values = c(
    "Up-regulated in Stool" = "red",
    "Up-regulated in Vomit" = "blue"
  ), guide = "none") +
  
  theme_light(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Genomic Position", y = "log2 Fold Change") +
  facet_wrap(~ Chr, scales = "free_x", ncol = 1)

genome_position
ggsave(filename = "1.Final_Paired_results/genome_position_paired.png", width = 10, height = 5.5, units = "in",dpi= 1000) 

#***********************************************************************************





