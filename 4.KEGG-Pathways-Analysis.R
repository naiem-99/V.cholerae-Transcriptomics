# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 15 │ FIGURE 4 — KEGG GSEA BARPLOT                                 #
# ─────────────────────────────────────────────────────────────────────────── #

geneList_gsea <- DE_results %>%arrange(desc(logFC)) %>%dplyr::select(GENEID, logFC) %>%tibble::deframe()

options(timeout = 1000)

rds_path <- "1.Final_Paired_results/gseKEGG_result.rds"

gseKEGG_res <- gseKEGG(
  geneList     = geneList_gsea,
  organism     = "vch",
  keyType      = "kegg",
  minGSSize    = 10,
  pvalueCutoff = 1,
  verbose      = FALSE)

saveRDS(gseKEGG_res, rds_path)

GSEA_table <- as.data.frame(gseKEGG_res@result) %>% mutate(Group = if_else(NES > 0, "Stool", "Vomit"))

GSEA_table[1:1,]

write.csv(GSEA_table,"1.Final_Paired_results/GSEA_result.csv",row.names = FALSE)

#--------------------------------------------------------------------------------#

gsea_plot_df <- GSEA_table %>%
  mutate(Pathway = gsub(" - .*$", "", Description),
         logP    = -log10(pvalue)) %>%
  filter(pvalue < 0.05) %>%
  {
    bind_rows(
      filter(., NES > 0) %>% arrange(pvalue, desc(NES)),   # most significant first within Stool
      filter(., NES < 0) %>% arrange(pvalue,      NES)     # most significant first within Vomit
    )
  } %>%
  arrange(NES)

gsea_plot_df[1:1,]

dim(gsea_plot_df)
table(gsea_plot_df$Group)

write.csv(gsea_plot_df,"1.Final_Paired_results/gsea_plot_df.csv",row.names = FALSE)
gsea_plot_df <-read.csv("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/1.Final_Paired_results/gsea_plot_df.csv",check.names = F)

#-----------------------------------------------------------------------------------------
fig_gsea <- ggplot(gsea_plot_df,
                   aes(x = reorder(Pathway, NES), y = NES, fill = logP))+
                   #aes(x = reorder(Pathway, logP * sign(NES)), y = NES, fill = logP))+
  
  geom_col(width = 0.75, colour = "black", linewidth = 0.3) +
  
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "black") +
  
  coord_flip(clip = "off") +
  
  annotate("text",
           x = Inf, y = max(gsea_plot_df$NES) * 0.5,
           label = "Stool", size = 5, fontface = "bold",
           colour = "black", hjust = 0.5, vjust = -0.8) +
  
  annotate("text",
           x = Inf, y = min(gsea_plot_df$NES) * 0.5,
           label = "Vomit", size = 5, fontface = "bold",
           colour = "black", hjust = 0.5, vjust = -0.8) +
  
  scale_fill_gradient(low  = "#fff5eb", high = "#6a00a8",
                      name = expression(-log[10]~italic(p)~"-value")) +
  
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  
  labs(x = NULL, y = "Normalised Enrichment Score (NES)") +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title.x      = element_text(face = "bold", size = 14),
    axis.text         = element_text(size = 11, colour = "black"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(size = 10),
    legend.position   = "right",
    legend.key.width  = unit(0.35, "cm"),
    legend.key.height = unit(0.9,  "cm"),
    plot.caption      = element_text(size = 9, face = "italic"),
    plot.margin       = margin(t = 25, r = 10, b = 5, l = 5, unit = "pt"))

fig_gsea

ggsave("1.Final_Paired_results/Fig4_GSEA.png",fig_gsea, width = 10, height = 12, units = "in", dpi = 1200)
#---------------------------------------------------------------------------------------------------------#
