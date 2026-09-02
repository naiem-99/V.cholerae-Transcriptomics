 ═══════════════════════════════════════════════════════════════════════════ #
#  SECTION 17 │ SUPPLEMENTAL FIGURE — VIRULENCE GENE EXPRESSION              #
# ═══════════════════════════════════════════════════════════════════════════ #
vir_genes  <- c("VC_1457", "VC_1456", "VC_0828", "VC_0838")
vir_labels <- c("VC_1457" = "ctxA", "VC_1456" = "ctxB",
                "VC_0828" = "tcpA", "VC_0838" = "toxT")

patient_cols <- c("P01" = "#2ca02c",
                  "P02" = "#1f77b4",
                  "P03" = "#ff7f0e",
                  "P04" = "#9467bd",
                  "P05" = "#17becf")

# Pull raw counts → log2(x+1) → long format → attach metadata
vir_df <- Merged_Data2 %>% dplyr::select(GENEID,
                S01, S02, S03, S04, S05,
                V01, V02, V03, V04, V05,GeneName_Merged)%>%
         filter(GENEID %in% vir_genes) %>%
         mutate(across(S01:V05, ~ log2(. + 1))) %>%
         tidyr::pivot_longer(cols = -c(GENEID,GeneName_Merged), names_to = "Sample", values_to = "log2count") %>%
         left_join(Paired_Sample, by = "Sample") %>%
         mutate(Gene = factor(vir_labels[GENEID], levels = c("ctxA","ctxB","tcpA","toxT")),
         Condition = factor(Sample_type, levels = c("Vomit","Stool")))



vir_df[1:10,]

fig_virulence <- ggplot(vir_df,aes(x     = Condition,
                            y     = log2count,
                            colour = Patient,
                            group  = Patient)) +
  
  # Black line connecting each patient's Vomit → Stool
  geom_line(colour = "black", linewidth = 0.5, alpha = 0.7) +
  # Triangle (17) for Vomit, circle (16) for Stool; coloured by patient
  geom_point(aes(shape = Condition), size = 4, alpha = 0.95) +
  scale_shape_manual(values = c("Vomit" = 17, "Stool" = 16), guide = "none") +
  scale_colour_manual(values = patient_cols, guide = "none") +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  labs(x       = NULL,y       = "Expression\n(log2[count + 1])")+
  theme_classic(base_size = 12) +
  theme(
    axis.title.y     = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 11, colour = "black"),
    axis.text.y      = element_text(size = 10, colour = "black"),
    strip.text       = element_text(size = 12, face = "italic"),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position  = "none",
    plot.caption     = element_text(size = 9, face = "italic",
                                    hjust = 0, colour = "grey40"))

fig_virulence
ggsave("1.Final_Paired_results/FigS1_Virulence_Genes01.png",fig_virulence, width = 10, height = 5, units = "in", dpi = 1200)

#-----Run Wilcoxon signed-rank test for each gene (paired = TRUE)----------

vir_stats <- vir_df %>%group_by(Gene)%>%summarise(p_value = wilcox.test(
            log2count[Condition == "Vomit"],
            log2count[Condition == "Stool"],
            paired    = TRUE,
            exact     = FALSE)$p.value,
    mean_Vomit = mean(log2count[Condition == "Vomit"]),
    mean_Stool = mean(log2count[Condition == "Stool"]),
    .groups = "drop") %>% mutate(p_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ paste0("p = ", round(p_value, 2))))

print(vir_stats)
#-------------------------------------------------------------#
fig_virulence <- fig_virulence +geom_text(
                  data  = vir_stats,
                  mapping     = aes(x = 1.5,
                                    y = Inf, 
                                label = p_label),
                  inherit.aes = FALSE,
                  size        = 3.8,
                  vjust       = 1.5,
                  fontface    = "bold",
                  colour      = "black")
#---------------------------------------------------------------#
fig_virulence
#-------------------------------------------------------------#
ggsave("1.Final_Paired_results/FigS1_Virulence_Genes.png",fig_virulence, width = 10, height = 5, units = "in", dpi = 1200)
