#*#***************************PCA pot************************************************
#Run PCA
#devtools::install_github("kassambara/factoextra")
library(factoextra)
library(ggrepel)

#------------------------------------------------------------------------------------
logCPM <- v$E   # log2 counts-per-million matrix
logCPM
pca <- prcomp(t(logCPM), scale. = F,center = TRUE)
pca_df <- data.frame(pca$x, Paired_Sample)
library(dplyr)
pca_df <-pca_df%>%mutate(Condition=as.factor(Sample_type))
#https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#google_vignette
groups <- as.factor(pca_df$Condition)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  ## Points
  geom_point(size = 4, alpha = 0.9) + geom_text_repel(aes(label = Sample), size = 3, show.legend = FALSE) +
  ## FILLED ellipse (surface)
  stat_ellipse(
    aes(group = Condition, fill = Condition),
    geom = "polygon",
    type = "norm",
    level = 0.8,
    alpha = 0.25,
    color = NA,
    show.legend = FALSE) +
  ## OUTLINE ellipse
  stat_ellipse(aes(group = Condition),
    geom = "path",type = "norm",level = 0.8,
    linewidth = 0.4,color = "black",show.legend = FALSE) +
  ## Axes
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  ## Labels
  xlab(paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "% variance)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "% variance)")) +
  ## Colors
  scale_color_manual(values = c("Stool" = "#1874CD", "Vomit" = "#FF69B4")) +
  scale_fill_manual(values = c("Stool" = "#1874CD", "Vomit" = "#FF69B4")) +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.line = element_blank()) +
  labs(color = "Groups", fill = "Groups")
ggsave(filename = "1.Final_Paired_results/pca_paired.png", width = 4, height = 4, units = "in",dpi= 1000) 

