library(Seurat)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)

monocytes <- readRDS("Monocytes_without_48h.rds")

gene_list <- c("Ly6c2", "Sell", "Mmp8", "Dusp16", "G0s2", "Atrnl1", "F13a1", "Cd83", 
               "Slc12a2", "Fabp4", "Ikzf3", "Tgfbr3", "Kit", "CD177", "Cd36", "Eno3", 
               "Tgm2", "Vegfa", "Fn1", "Spn", "Ccr2", "Treml4", "Cx3cr1", "Siglec1", 
               "Itgax", "Itgam")

marker_genes <- c("Sell", "H2-Ab1", "Ly6c2", "Nos2")
VlnPlot(monocytes, features = marker_genes, pt.size = 0)
VlnPlot(monocytes, features = marker_genes, pt.size = 0, slot = "counts", log = TRUE)
FeaturePlot(monocytes, features = marker_genes)

gene_list <- c("Ly6c2", "Sell", "Mmp8", "Dusp16", "G0s2", "Atrnl1", "F13a1", "Cd83", 
               "Slc12a2", "Fabp4", "Ikzf3", "Tgfbr3", "Kit", "CD177", "Cd36", "Eno3", 
               "Tgm2", "Vegfa", "Fn1", "Spn", "Ccr2", "Treml4", "Cx3cr1", "Siglec1", 
               "Itgax", "Itgam")

marker_genes <- c("Sell", "H2-Ab1", "Ly6c2", "Nos2")
VlnPlot(monocytes, features = marker_genes, pt.size = 0)
VlnPlot(monocytes, features = marker_genes, pt.size = 0, slot = "counts", log = TRUE)
FeaturePlot(monocytes, features = marker_genes)

RidgePlot(monocytes, features = marker_genes)
DotPlot(monocytes, features = marker_genes) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

gene_list <- unique(gene_list)
gene_list <- intersect(gene_list, rownames(monocytes))

dir.create("gene_plots", showWarnings = FALSE)
dir.create("gene_stats", showWarnings = FALSE)

plot_gene_expression <- function(gene) {
  if (!gene %in% rownames(monocytes)) {
    message(paste("Gene", gene, "not found in monocytes. Skipping..."))
    return(NULL)
  }
  
  p_umap <- FeaturePlot(monocytes, features = gene, split.by = "time_point",
                        order = TRUE, cols = c("lightgrey", "red"), ncol = 4) + 
    ggtitle(paste(gene, "Expression (UMAP)"))
  
  p_vln <- VlnPlot(monocytes, features = gene, split.by = "time_point",
                   group.by = "seurat_clusters", pt.size = 0.1, ncol = 1) +
    ggtitle(paste(gene, "Expression by Subtype and Time"))
  
  p_combined <- p_umap / p_vln + plot_layout(heights = c(1, 2))
  
  ggsave(paste0("gene_plots/", gene, "_expression.png"), plot = p_combined,
         width = 10, height = 8, dpi = 300)
  
  return(p_combined)
}

lapply(gene_list, plot_gene_expression)

calculate_gene_stats <- function(gene) {
  if (!gene %in% rownames(monocytes)) {
    return(data.frame(Gene = gene, Error = "Not found"))
  }
  
  stats <- monocytes@meta.data %>%
    mutate(expression = FetchData(monocytes, vars = gene)[, 1]) %>%
    group_by(seurat_clusters, time_point) %>%
    summarise(
      gene = gene,
      median_expression = median(expression),
      mean_expression = mean(expression),
      expr_cell_ratio = sum(expression > 0) / n(),
      .groups = "drop"
    )
  
  write.csv(stats, paste0("gene_stats/", gene, "_stats.csv"), row.names = FALSE)
  
  return(stats)
}

all_stats <- lapply(gene_list, calculate_gene_stats) %>% bind_rows()
write.csv(all_stats, "gene_stats/all_genes_summary.csv", row.names = FALSE)

expr_data <- FetchData(monocytes, vars = c(gene_list, "time_point"))
time_stats <- expr_data %>%
  pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression") %>%
  group_by(gene, time_point) %>%
  summarise(
    mean_exp = median(expm1(expression)),
    pct_exp = median(expression > 0) * 100,
    cells = n(),
    .groups = "drop"
  )

heatmap_data <- time_stats %>%
  select(gene, time_point, mean_exp) %>%
  pivot_wider(names_from = time_point, values_from = mean_exp) %>%
  column_to_rownames("gene")

p_heatmap <- pheatmap::pheatmap(
  heatmap_data,
  scale = "row",
  color = viridis::viridis(100),
  main = "Median Gene Expression Across Time Points",
  fontsize_row = 8,
  fontsize_col = 10,
  cluster_cols = FALSE
)

p_curves <- ggplot(time_stats, aes(x = time_point, y = mean_exp, group = gene, color = gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Time Point", y = "Median Expression (log10)") +
  theme_minimal() +
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(
    data = filter(time_stats, time_point == last(time_point)),
    aes(label = gene), direction = "y", hjust = 0, segment.size = 0.2
  )

p_bubble <- ggplot(time_stats, aes(x = time_point, y = gene, size = pct_exp, color = mean_exp)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(2, 10)) +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Time Point", y = "Gene", size = "% Cells Expressing", color = "Median Expression") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("expression_heatmap_median.png", p_heatmap, width = 8, height = 10, dpi = 300)
ggsave("expression_curves_median.png", p_curves, width = 8, height = 6, dpi = 300)
ggsave("expression_bubble_median.png", p_bubble, width = 8, height = 10, dpi = 300)

Idents(monocytes) <- "time_point"
time_markers <- FindAllMarkers(monocytes, features = gene_list, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(time_stats, "gene_expression_time_stats.csv", row.names = FALSE)

time_levels <- sort(unique(monocytes$time_point))
monocytes$time_point <- factor(monocytes$time_point, levels = time_levels)

expr_data <- FetchData(monocytes, vars = c(gene_list, "time_point"))
baseline <- expr_data[monocytes$time_point == time_levels[1], ]

dir.create("timepoint_foldchange_plots", showWarnings = FALSE)

for (tp in time_levels[-1]) {
  fc_data <- sapply(gene_list, function(gene) {
    mean(expr_data[expr_data$time_point == tp, gene]) / mean(baseline[, gene])
  })
  
  plot_data <- data.frame(
    gene = factor(gene_list, levels = gene_list),
    log2_fold_change = log2(fc_data)
  )
  
  p <- ggplot(plot_data, aes(x = 1, y = gene, fill = log2_fold_change)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-2, 2), na.value = "grey50", name = "Log2 Fold Change") +
    geom_text(aes(label = round(log2_fold_change, 2)), color = "black", size = 3) +
    labs(title = paste("Expression Change:", time_levels[1], "vs", tp), x = NULL, y = "Gene") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(paste0("timepoint_foldchange_plots/foldchange_", tp, ".png"), plot = p,
         width = 6, height = max(6, length(gene_list) * 0.3), dpi = 300)
}

all_fc <- sapply(time_levels[-1], function(tp) {
  sapply(gene_list, function(gene) {
    mean(expr_data[expr_data$time_point == tp, gene]) / mean(baseline[, gene])
  })
})

all_plot_data <- as.data.frame(all_fc) %>%
  mutate(gene = gene_list) %>%
  pivot_longer(cols = -gene, names_to = "time_point", values_to = "fold_change") %>%
  mutate(log2_fc = log2(fold_change))

p_summary <- ggplot(all_plot_data, aes(x = time_point, y = gene, fill = log2_fc)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-2, 2), name = "Log2 FC") +
  geom_text(aes(label = round(log2_fc, 1)), size = 3, color = "black") +
  labs(title = paste("All Time Points vs", time_levels[1]), x = "Time Point", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("timepoint_foldchange_plots/all_timepoints_summary.png", p_summary,
       width = 8, height = max(8, length(gene_list) * 0.3), dpi = 300)

expr_data <- AverageExpression(monocytes, assays = "RNA", features = gene_list, group.by = "time_point")$RNA
cv_data <- apply(expr_data, 1, function(x) sd(x)/mean(x))
stable_genes <- names(sort(cv_data)[1:10])

expr_data <- FetchData(monocytes, vars = gene_list)
zscore_data <- t(scale(t(expr_data)))
zscore_data <- cbind(zscore_data, time_point = monocytes$time_point)

zscore_long <- zscore_data %>%
  as.data.frame() %>%
  pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "zscore") %>%
  group_by(gene, time_point) %>%
  summarise(mean_zscore = mean(zscore), .groups = "drop")

ggplot(zscore_long, aes(x = time_point, y = mean_zscore, group = gene, color = gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "Time Point", y = "Z-score Normalized Expression", 
       title = "Gene Expression Trends (Z-score Normalized)") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggrepel::geom_text_repel(
    data = filter(zscore_long, time_point == last(time_point)),
    aes(label = gene), direction = "y", hjust = 0
  )

time_levels <- sort(unique(monocytes$time_point))
monocytes$time_point <- factor(monocytes$time_point, levels = time_levels)

expr_data <- FetchData(monocytes, vars = gene_list)
baseline <- expr_data[monocytes$time_point == time_levels[1], ]
fold_change <- sapply(gene_list, function(gene) {
  sapply(time_levels, function(tp) {
    mean(expr_data[monocytes$time_point == tp, gene]) / mean(baseline[, gene])
  })
})

fold_change_long <- as.data.frame(fold_change) %>%
  mutate(time_point = time_levels) %>%
  pivot_longer(cols = -time_point, names_to = "gene", values_to = "fold_change")

ggplot(fold_change_long, aes(x = time_point, y = gene, fill = log2(fold_change))) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2, 2)) +
  labs(x = "Time Point", y = "Gene", fill = "Log2 Fold Change",
       title = "Expression Change Relative to Baseline Time") +
  theme_minimal()

time_stats <- FetchData(monocytes, vars = c(gene_list, "time_point")) %>%
  pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression") %>%
  group_by(gene, time_point) %>%
  summarise(mean_exp = mean(expression), .groups = "drop")

ggplot(time_stats, aes(x = time_point, y = mean_exp)) +
  geom_line(aes(group = 1), color = "steelblue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  labs(x = "Time Point", y = "Mean Expression (log-normalized)",
       title = "Individual Gene Expression Trends Over Time") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Individual Gene Expression Trends Over Time.png")

expr_data <- AverageExpression(monocytes, assays = "RNA", features = gene_list, group.by = "time_point")$RNA

calculate_stability <- function(expr_matrix) {
  apply(expr_matrix, 1, function(x) {
    cv <- sd(x) / mean(x)
    return(1 - cv)
  })
}

stability_all <- calculate_stability(expr_data)

stability_df <- data.frame(
  Gene = names(stability_all),
  Stability = stability_all,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Stability))

head(stability_df)

top_genes <- head(stability_df$Gene, 25)
scaled_data <- t(scale(t(expr_data[top_genes, ])))

hh <- Heatmap(
  scaled_data,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 12),
  column_title = "Stable Gene Expression Across Time",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title_position = "leftcenter-rot", legend_height = unit(4, "cm")))

dev.off()
draw(hh)
ggsave("gene_stability_heatmap.png", width = 12, height = 8)

stability_df <- stability_df[order(-stability_df$Stability), ]
gene_order <- stability_df$Gene
scaled_data_ordered <- scaled_data[gene_order, ]

time_order <- c("0h", "4h", "24h")
scaled_data_ordered <- scaled_data_ordered[, time_order]

ht <- Heatmap(
  scaled_data_ordered,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = gene_order,
  column_order = time_order,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 12),
  column_title = "Gene Expression by Time (Ordered by Stability)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title_position = "leftcenter-rot", legend_width = unit(6, "cm")))

pdf("Stability_Ordered_Heatmap.pdf", width = 8 + ncol(scaled_data)*0.3, height = 6 + nrow(scaled_data)*0.4)
draw(ht)
dev.off()

png("Stability_Ordered_Heatmap.png", width = 2000, height = 1200 + nrow(scaled_data)*40, res = 300)
draw(ht)
dev.off()

sessionInfo()
