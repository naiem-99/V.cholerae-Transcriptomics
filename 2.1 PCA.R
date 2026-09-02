
# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 12 │ FIGURE 1 — PCA                                                #
# ─────────────────────────────────────────────────────────────────────────── #

logCPM   <- v$E
pca_out  <- prcomp(t(logCPM), center = TRUE, scale. = FALSE)
var_pct  <- round(100 * summary(pca_out)$importance[2, ], 1)

pca_df <- data.frame(
  PC1       = pca_out$x[, 1],
  PC2       = pca_out$x[, 2],
  Sample    = Paired_Sample$Sample,
  Condition = factor(Paired_Sample$Sample_type, levels = c("Vomit", "Stool")),
  Patient   = Paired_Sample$Patient)

CI_LEVEL <- 0.95

pca_df[1:5,]

fig_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition)) +stat_ellipse(aes(group = Condition, fill = Condition),
               geom = "polygon", type = "norm", level = CI_LEVEL,
               alpha = 0.20, colour = NA, show.legend = FALSE) +
  
  stat_ellipse(aes(group = Condition),
               geom = "path", type = "norm", level = CI_LEVEL,
               linewidth = 0.5, colour = "black", show.legend = FALSE) +
  
  geom_point(size = 4, alpha = 0.9) +
  
  geom_text_repel(aes(label = Sample),
                  size = 3, colour = "black",
                  box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = Inf, seed = 42, show.legend = FALSE) +
  
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  
  scale_colour_manual(values = c("Stool" = "#FF1493", "Vomit" = "#26b3ff"),
                      name = "Group") +
  scale_fill_manual(values   = c("Stool" = "#FF1493", "Vomit" = "#26b3ff"),
                    guide = "none") +
  
  labs(x = paste0("PC1 (", var_pct[1], "% variance)"),
       y = paste0("PC2 (", var_pct[2], "% variance)")) +
  
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 10, colour = "black"),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))


fig_pca
ggsave("1.Final_Paired_results/Fig1_PCA.png",fig_pca, width = 5, height = 6, dpi = 800)

#----------------------------------------------------------------------

#---------------------- Perform permanova Test ------------------------
logCPM_t <- t(v$E)   # transpose: samples = rows, genes = columns

# With patient-restricted permutations (correct for paired design)
set.seed(123)
perm <- how(blocks = Paired_Sample$Patient, nperm = 999)
permanova_result <- adonis2(logCPM_t ~ Sample_type,data = Paired_Sample,permutations = perm,method = "euclidean")
print(permanova_result)

# Extract values from PERMANOVA result
perm_R2 <- round(permanova_result$R2[1], 3)
perm_p  <- round(permanova_result$`Pr(>F)`[1], 3)


#---------------------2nd try
set.seed(123)
#permanova_result <- adonis2(logCPM_t ~ Sample_type,data = Paired_Sample,permutations = 31,strata = Paired_Sample$Patient, method = "euclidean")
#print(permanova_result)
#-------------------------------------------------------------------------

fig_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition)) +
  
  stat_ellipse(aes(group = Condition, fill = Condition),
               geom = "polygon", type = "norm", level = CI_LEVEL,
               alpha = 0.20, colour = NA, show.legend = FALSE) +
  
  stat_ellipse(aes(group = Condition),
               geom = "path", type = "norm", level = CI_LEVEL,
               linewidth = 0.5, colour = "black", show.legend = FALSE) +
  
  geom_point(size = 4, alpha = 0.9) +
  
  geom_text_repel(aes(label = Sample),
                  size = 3, colour = "black",
                  box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = Inf, seed = 42, show.legend = FALSE) +
  
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  
  annotate("text",
           x        = Inf,
           y        = Inf,
           label    = paste0("R² = ", perm_R2,
                             "\np = ", perm_p),
           hjust    = 1.1,
           vjust    = 1.4,
           size     = 3.2,
           colour   = "black",
           fontface = "italic")+
  
  scale_colour_manual(values = c("Stool" = "#FF1493", "Vomit" = "#26b3ff"),
                      name = "Group") +
  scale_fill_manual(values = c("Stool" = "#FF1493", "Vomit" = "#26b3ff"),
                    guide = "none") +
  
  labs(x = paste0("PC1 (", var_pct[1], "% variance)"),
       y = paste0("PC2 (", var_pct[2], "% variance)")) +
  
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10, colour = "black"),
    legend.position = "top",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    panel.border    = element_rect(colour = "black", fill = NA,
                                   linewidth = 0.6))
fig_pca
ggsave("1.Final_Paired_results/Fig1_PCA_p_val.png",fig_pca, width = 5, height = 6, dpi = 800)
#-------------------------------------------------------------------------------------------#
