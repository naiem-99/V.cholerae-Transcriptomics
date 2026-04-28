
#----------------------------- KEGG Geneset Enrichment Analysis----------------------------------------
#load required library 
library(clusterProfiler)
library(ggplot2)
library(scales)  # for nice formatting of colors
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(enrichplot)
#------------------------------------------------
# 1. Create a named numeric vector of logFC
geneList <- setNames(Stool_vs_Vomit_paired$logFC, Stool_vs_Vomit_paired$GENEID)

# 2. Sort in decreasing order (highest logFC first)
geneList <- sort(geneList, decreasing = TRUE)
options(timeout = 1000)
#---------------Run GSEA for kegg pathway-----------------------------------

gseKEGG <- gseKEGG(geneList=geneList,organism = "vch",keyType = "kegg", minGSSize = 10,verbose = T)

saveRDS(gseKEGG,"1.Final_Paired_results/kk2_gseKEGG_stool_vomit_paired.rds")

#----------------------------------------------------------------------------
gseKEGG <-readRDS("1.Final_Paired_results/kk2_gseKEGG_stool_vomit_paired.rds")
GSEA_result <-as.data.frame(gseKEGG@result)
head(GSEA_result)
write.csv(GSEA_result,"1.Final_Paired_results/GSEA_result.csv")



gseaplot(gseKEGG,geneSetID = 3,title = gseKEGG@result$Description[3])
gseaplot( gseKEGG, geneSetID = "vch02020", title = "Two-component system")


#-------------------------------------------------------------------------------------

# Convert enrichment results to a data frame
df_gseKEGG_stool_vomit_paired <- as.data.frame(gseKEGG)

Processed_gseKEGG <- df_gseKEGG_stool_vomit_paired %>% mutate(Description = gsub(" - .*", "", Description)  # remove everything after " - ")
Processed_gseKEGG

#Processed_gseKEGG$qvalue_scaled <- Processed_gseKEGG$qvalue / max(abs(Processed_gseKEGG$qvalue))

y<-ggplot(Processed_gseKEGG, aes(x = reorder(Description, NES), y = NES, fill = log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#1E90FF" ,high = "#BA55D3",name = "log10(p.adjust)")+
  labs(
    title = "Enriched Pathways in Stool and Vomit ",
    x = "KEGG Pathway",
    y = "NES"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
x
y
#---------------------------2nd try--------------------------
fill_range <- range(log10(Processed_gseKEGG$p.adjust), na.rm = TRUE)

y <- ggplot(
  Processed_gseKEGG,
  aes(
    x = reorder(Description, NES),
    y = NES,
    fill = log10(p.adjust)
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(
    low = "#1E90FF",
    high = "#BA55D3",
    breaks = fill_range,     # 🔹 show only min & max
    labels = round(fill_range, 2),
    name = "log10(p.adjust)"
  ) +
  labs(
    title = "Enriched Pathways in Stool and Vomit ",
    x = "KEGG Pathway",
    y = "NES"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black", size = 1)
  )

y

ggsave(filename = "1.Final_Paired_results/Stool_vs_Vomit_GSEA_paired_2.png", width = 8, height = 6, units = "in",dpi= 1000) 


#-----------------------------End of KEGG pathways----------------------------

