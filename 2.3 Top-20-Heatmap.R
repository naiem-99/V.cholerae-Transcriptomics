# ─────────────────────────────────────────────────────────────────────────── #
#  SECTION 14 │ FIGURE 3 — HEATMAP (Top 20 by adj.P.Val, Vomit ----> Stool order) #
# ─────────────────────────────────────────────────────────────────────────── #

# ── Top 20 genes by adj.P.Val only ────────────────────────────────────────────
hm_genes <- DE_results %>%arrange(adj.P.Val) %>%head(20)

write.csv(hm_genes, "1.Final_Paired_results/hm_genes_20_overall.csv",row.names = FALSE)

# ── Build expression matrix ───────────────────────────────────────────────────
hm_data <- Merged_Data2 %>%
  dplyr::select(GENEID,
                S01, S02, S03, S04, S05,
                V01, V02, V03, V04, V05,
                `N16961 annotation`) %>%
  distinct(GENEID, .keep_all = TRUE) %>%
  filter(GENEID %in% hm_genes$GENEID) %>%
  left_join(hm_genes %>% dplyr::select(GENEID, Expression, Symbol, Heatmap_Gene_Name),
            by = "GENEID") %>%
  left_join(gene_symbols %>% dplyr::select(GENEID, Description),
            by = "GENEID") %>%
  mutate(Function_label = Heatmap_Gene_Name)

# ── Reordered columns: Vomit first, Stool second ─────────────────────────────
expr_raw <- hm_data %>%
  column_to_rownames("GENEID") %>%
  dplyr::select(V01, V02, V03, V04, V05, S01, S02, S03, S04, S05) %>%
  as.matrix()

# Raw z-scores — no capping
expr_scaled <- t(scale(t(log2(expr_raw + 1))))

# ── Colour scale spans actual z-score range ──────────────────────────────────
zmax <- max(abs(expr_scaled), na.rm = TRUE)

col_fun <- colorRamp2(
  c(-zmax, -zmax/2, 0, zmax/2, zmax),
  c("#053061", "#2c7bb6", "#ffffff", "#ff8080", "#ff0000"))

# ── Top annotation: Vomit block first, Stool block second ────────────────────
col_annotation <- HeatmapAnnotation(
  Group = factor(
    c("Vomit","Vomit","Vomit","Vomit","Vomit",
      "Stool","Stool","Stool","Stool","Stool"),
    levels = c("Vomit","Stool")
  ),
  col   = list(Group = c("Stool" = "#FF1493", "Vomit" = "#26b3ff")),
  border                  = TRUE,
  simple_anno_size        = unit(0.35, "cm"),
  show_annotation_name    = FALSE,
  annotation_legend_param = list(direction = "horizontal", nrow = 1)
)

# ── Right annotation: gene function text only ────────────────────────────────
row_annotation <- rowAnnotation(
  Function = anno_text(hm_data$Function_label,
                       gp = gpar(fontsize = 10), just = "left", location = 0),
  width = unit(10.8, "cm"))

# ── Build heatmap ─────────────────────────────────────────────────────────────
ht <- Heatmap(
  expr_scaled,
  name                 = "Gene expression\n(Scaled log2[count+1])",
  col                  = col_fun,
  top_annotation       = col_annotation,
  right_annotation     = row_annotation,
  row_split            = hm_data$Expression,
  row_title_gp         = gpar(fontsize = 0, fontface = "bold"),
  row_title_rot        = 0,
  cluster_rows         = TRUE,
  cluster_columns      = FALSE,
  show_row_dend        = FALSE,
  show_column_dend     = FALSE,
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 9),
  border               = TRUE,
  rect_gp              = gpar(col = "white", lwd = 0.6),
  heatmap_legend_param = list(
    title_gp       = gpar(fontsize = 9, fontface = "bold"),
    labels_gp      = gpar(fontsize = 8),
    at             = c(-zmax, 0, zmax),
    labels         = c("Low", "Mid", "High"),
    direction      = "horizontal",
    legend_width   = unit(4, "cm"),
    title_position = "topcenter"))

# ── Save ──────────────────────────────────────────────────────────────────────
png("1.Final_Paired_results/Fig3_Heatmap.png",width = 10.5, height = 8.5, units = "in", res = 1200)
P <- draw(ht,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  merge_legends          = TRUE,
  padding                = unit(c(2, 2, 2, 2), "mm"))
P
dev.off()
