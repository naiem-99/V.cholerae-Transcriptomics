
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 16 │ FIGURE 5 — GENOMIC POSITION PLOT                             #
#  Pathogenicity islands shown as shaded regions with legend labels           #
# ─────────────────────────────────────────────────────────────────────────── #

gtf_all   <- rtracklayer::import("D:/ICDDRB_Feb/1.Sadia_Apu_RNAseq/1.Final_VC_RNAseq_10_11/5.Final_Analysis_paired/Raw_Files/GCA_013085075.1_ASM1308507v1_genomic.gtf")

gtf_genes <- gtf_all[gtf_all$type == "gene"]
gene_coords <- data.frame(
  Chr           = as.character(seqnames(gtf_genes)),
  Start         = start(gtf_genes),
  Stop          = end(gtf_genes),
  Width         = width(gtf_genes),
  Strand        = as.character(strand(gtf_genes)),
  Mid_gene      = (start(gtf_genes) + end(gtf_genes)) / 2,
  as.data.frame(mcols(gtf_genes)),
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(Geneid = gene_id)

gene_id_bridge <- Merged_Data2 %>%dplyr::select(Geneid, GENEID)

gene_coords_full <- gene_id_bridge %>%inner_join(gene_coords, by = "Geneid")

res_genome <- DE_results %>%left_join(gene_coords_full, by = "GENEID") %>%
                          filter(!is.na(Mid_gene)) %>%
                          mutate(Expression = factor(
                          as.character(Expression),
                         levels = c("Increased in Vomit", "NS", "Increased in Stool")))

write.csv(res_genome, "1.Final_Paired_results/res_genome.csv", row.names = FALSE)


res_genome <-read.csv("1.Final_Paired_results/res_genome.csv", check.names = F)
# ── Top gene labels ──────────────────────────────────────────────────────────

label_genome <- res_genome %>%filter(Expression %in% c("Increased in Stool", "Increased in Vomit")) 


# ── Pathogenicity island coordinates ─────────────────────────────────────────
# Computed from res_genome using Pathogenicity Island column from Table_S2
# Factor levels control legend order

island_colours <- c("VPI-I"  = "#fbb4d4","VPI-II" = "#e7298a",
                    "VSP-I"  = "#b3cde3","VSP-II" = "#1874cd")

#---------------------------------------------------------------------

#-----------------------------------------------------------------------
# Check which boundary genes exist in res_genome
boundary_check <- c( # For VPI-1
                     "VC_0817","VC_0818",
                    "VC_0846","VC_0847",
                    # For VPI-2 
                    "VC_1758","VC_1759",
                    "VC_1808","VC_1809",
                    # For VSP-1
                    "VC_0175","VC_0176",
                    "VC_0184","VC_0185",
                    # For VSP-2
                    "VC_0490","VC_0491",
                    "VC_0515","VC_0516")

res_genome_island<-res_genome %>%
  filter(GENEID %in% boundary_check) %>%
  dplyr::select(GENEID, Chr, Start, Stop, `Pathogenicity Island`) %>%
  arrange(GENEID)%>%as.data.frame()


res_genome_island
write.csv(res_genome_island, "1.Final_Paired_results/res_genome_island.csv", row.names = FALSE)

islands <- data.frame(Chr    = "CP047295.1",Island = factor(
          c("VPI-I",   "VPI-II",  "VSP-I",   "VSP-II"),levels = c("VPI-I","VPI-II","VSP-I","VSP-II")),
          xmin = c(2219669,  1129597,  2948154,  2580596),   
          xmax = c(2257361,  1186335,  2960979,  2607114) ) %>%mutate(x = (xmin + xmax) / 2)

print(islands)
write.csv(islands, "1.Final_Paired_results/islands.csv", row.names = FALSE)
#-------------------------------------------------------------------------------
#---------------------------------------------------------------------------

#---------------------------------------------------------------
#PLOT COLOURS & THEME                                           #

COL_STOOL  <- "#FF1493"
COL_VOMIT  <- "#26b3ff"
COL_STABLE <- "#9EB3C2"
plot_theme <- theme_bw() +
  theme(
    axis.title       = element_text(size = 15, face = "bold"),
    axis.text        = element_text(size = 14, colour = "black"),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 16),
    legend.position  = "bottom",
    panel.grid.minor = element_blank())
# ── Plot ──────────────────────────────────────────────────────────────────────
fig_genome <- ggplot(res_genome,
                     aes(x = Mid_gene, y = logFC, colour = Expression)) +
  
  # Island shaded rectangles — fill mapped to Island for legend
  geom_rect(data        = islands,
            aes(xmin  = xmin,
                xmax  = xmax,
                ymin  = -Inf,
                ymax  = Inf,
                fill  = Island),
            inherit.aes = FALSE,
            alpha       = 0.30,
            colour      = NA) +
  
  # Island border lines
  geom_vline(data        = islands,
             aes(xintercept = xmin),
             inherit.aes = FALSE,
             linewidth   = 0.4,
             linetype    = "solid",
             colour      = "grey40",
             alpha       = 0.65) +
  
  geom_vline(data        = islands,
             aes(xintercept = xmax),
             inherit.aes = FALSE,
             linewidth   = 0.4,
             linetype    = "solid",
             colour      = "grey40",
             alpha       = 0.65) +
  
  # All genes — small semi-transparent points
  geom_point(size = 1.2, alpha = 0.5) +
  
  # Top-ranked highlighted genes — outlined circles
  geom_point(data   =label_genome,
             aes(fill = Expression),
             shape  = 21,
             colour = "black",
             size   = 3.2,
             stroke = 0.5,
             show.legend = FALSE) +
  
  # Gene name labels
  geom_text_repel(data           = label_genome,
                  mapping        = aes(x = Mid_gene, y = logFC,
                                       label = `Gene Symbol`),
                  inherit.aes    = FALSE,
                  size           = 3,
                  fontface       = "italic",
                  colour         = "black",
                  segment.size   = 0.3,
                  segment.colour = "grey40",
                  box.padding    = 0.4,
                  seed           = 42,
                  max.overlaps   = Inf) +
  
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey50") +
  
  # ── Colour scale: Expression (point outline colour) ────────────────────────
  scale_colour_manual(
    name   = "Expression",
    limits = c("Increased in Vomit", "NS", "Increased in Stool"),
    breaks = c("Increased in Vomit", "NS", "Increased in Stool"),
    values = c(
      "Increased in Stool" = COL_STOOL,
      "NS"                 = "grey60",
      "Increased in Vomit" = COL_VOMIT
    ),
    guide  = guide_legend(
      order        = 1,
      override.aes = list(size = 3, alpha = 1, shape = 16)) ) +
  
  # ── Fill scale: island shading (appears in legend) + point fills ──────────
  scale_fill_manual(name   = "Pathogenicity Island",
        values = c(
      # Island shading colours
      "VPI-I"  = "#fbb4d4",
      "VPI-II" = "#e7298a",
      "VSP-I"  = "#b3cde3",
      "VSP-II" = "#1874cd",
      # Significant gene point fills (same as colour — not shown in legend)
      "Increased in Stool" = COL_STOOL,
      "Increased in Vomit" = COL_VOMIT),
    breaks = c("VPI-I", "VPI-II", "VSP-I", "VSP-II"),   # only islands in legend
    guide  = guide_legend(
      order        = 2,
      override.aes = list(
        alpha  = 0.5,
        colour = NA,
        size   = 5,
        shape  = 22    # square swatch — matches shaded rectangle
      ) )) +
  
  scale_x_continuous(
    breaks = function(x) {
      if (max(x, na.rm = TRUE) < 1.5e6)
        c(0, 0.3e6, 0.6e6, 0.9e6, 1.0e6)
      else
        c(0, 1e6, 2e6, 3e6)
    },
    labels = function(x) paste0(x / 1e6, " Mb"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  facet_wrap(~ Chr, scales = "free_x", ncol = 1) +
  
  labs(
    x      = "Genomic position",
    y      = expression(log[2]~"Fold Change"),
    colour = "Expression"
  ) +
  
  plot_theme +
  theme(
    panel.spacing    = unit(0.6, "lines"),
    strip.text       = element_text(size = 11),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(0.5, "cm"),
    legend.spacing.x = unit(0.8, "cm")
  )

fig_genome

ggsave("1.Final_Paired_results/Fig5_GenomicPosition.png",fig_genome, width = 12, height = 6.5, units = "in", dpi = 1200)
#-------------------------------------------------------------------------------
