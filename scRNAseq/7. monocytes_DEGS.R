library(Seurat)
library(dplyr)
if (!require("tibble")) install.packages("tibble")
library(tibble)  


exclude_genes <- c("Itgax", "Itgam", "Siglec1") 

new_meta <- monocytes@meta.data %>%
  mutate(
    new_group = case_when(
      coexpressed == "Co-expressed" & 
        subtype %in% c("Nonclassical Monocytes") ~ "Coexp_Monocytes",
      
      coexpressed == "Co-expressed" ~ "Coexp_NonMonocytes",
      
      coexpressed == "Other" & 
        subtype %in% c("Classical Monocytes", "Intermediate Monocytes", "Nonclassical Monocytess") ~ paste0("Other_", subtype),

      TRUE ~ "Remove"
    )
  )


keep_cells <- rownames(new_meta)[new_meta$new_group != "Remove"]
new_counts <- monocytes[["RNA"]]$counts[, keep_cells]
new_data <- monocytes[["RNA"]]$data[, keep_cells]

new_counts <- new_counts[!rownames(new_counts) %in% exclude_genes, ]
new_data <- new_data[!rownames(new_data) %in% exclude_genes, ]


shared_genes <- intersect(rownames(new_counts), rownames(new_data))
new_counts <- new_counts[shared_genes, ]
new_data <- new_data[shared_genes, ]


new_seurat <- CreateSeuratObject(
  counts = new_counts,
  meta.data = new_meta[keep_cells, ],
  assay = "RNA"
)


new_seurat[["RNA"]]$data <- new_data


table(new_seurat$new_group)


print(paste("successfully remove:", 
            all(!exclude_genes %in% rownames(new_seurat))))


dim(new_seurat) 


new_seurat <- NormalizeData(new_seurat)


new_seurat <- FindVariableFeatures(new_seurat, 
                                   nfeatures = 2000,
                                   selection.method = "vst")


new_seurat <- ScaleData(new_seurat)
new_seurat <- RunPCA(new_seurat, npcs = 30)
new_seurat <- RunUMAP(new_seurat, dims = 1:20)
new_seurat <- FindNeighbors(new_seurat, dims = 1:20)
new_seurat <- FindClusters(new_seurat, resolution = 0.8)


group_colors <- c(
  "Coexp_Monocytes" = "#E15759",
  "Coexp_NonMonocytes" = "#FF9DA7",
  "Other_Class" = "#59A14F",
  "Other_Nonclass&Macro" = "#B383B9",
  "Other_Inter" = "#B6992D",
  "Other_NonClass" = "#C6307C"
)

DimPlot(new_seurat, 
        group.by = "new_group",
        cols = group_colors,
        pt.size = 0.05) +
  ggtitle("New Object Composition")
ggsave("new_seurat_for_gene_analysis.png", width = 7, height = 5)


saveRDS(new_seurat, "filtered_seurat_object.rds")


comparison_groups <- list(
  c("Coexp_Monocytes", "Other_Class"),
  c("Coexp_Monocytes", "Other_NonClass"),
  c("Coexp_Monocytes", "Coexp_NonMonocytes")
)


diff_results <- lapply(comparison_groups, function(pair) {
  FindMarkers(
    object = new_seurat,
    ident.1 = pair[1],
    ident.2 = pair[2],
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox",
    only.pos = FALSE
  ) %>%
    tibble::rownames_to_column("gene") %>%  
    mutate(
      comparison = paste(pair[1], "vs", pair[2]),
      direction = ifelse(avg_log2FC > 0, "Up", "Down")
    )
})


head(diff_results[[1]]) 


sig_genes <- do.call(rbind, diff_results) %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.58)


write.csv(sig_genes, "significant_DE_genes.csv")


library(EnhancedVolcano)
library(ggplot2)
library(patchwork)


sig_genes <- sig_genes %>%
  mutate(
    status = case_when(
      p_val_adj >= 0.05 ~ "NS",
      p_val_adj < 0.05 & avg_log2FC > 0.58 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.58 ~ "Down",
      TRUE ~ "NS"
    )
  )


status_colors <- c("NS" = "grey30", 
                   "Up" = "firebrick", 
                   "Down" = "royalblue")


plot_volcano <- function(data, title) {
  # p=0
  if(min(data$p_val_adj) == 0) {
    min_p <- min(data$p_val_adj[data$p_val_adj > 0], na.rm = TRUE)
    data$p_val_adj[data$p_val_adj == 0] <- min_p * 0.1
  }
  
  EnhancedVolcano(data,
                  lab = data$gene,
                  x = "avg_log2FC",
                  y = "p_val_adj",
                  title = title,
                  subtitle = "",
                  pCutoff = 0.05,
                  FCcutoff = 0.58,
                  cutoffLineType = "dashed",
                  cutoffLineCol = "grey40",
                  pointSize = 1.5,
                  labSize = 3.0,
                  colAlpha = 0.8,
                  legendPosition = "right",
                  colCustom = status_colors[data$status],  
                  selectLab = data %>%
                    arrange(p_val_adj) %>%
                    head(10) %>%
                    pull(gene)
  ) + 
    scale_color_manual(values = status_colors) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}


plot_list <- lapply(unique(sig_genes$comparison), function(comp) {
  df <- sig_genes %>% 
    filter(comparison == comp) %>%
    mutate(status = factor(status, levels = names(status_colors)))
  
  plot_volcano(df, comp)
})


print(plot_list[[1]]) 



test_df <- sig_genes %>% 
  filter(comparison == unique(sig_genes$comparison)[1])
p <- plot_volcano(test_df, "Test Volcano")
print(p)



length(status_colors[test_df$status]) == nrow(test_df)



combined_plot <- wrap_plots(plot_list, ncol = 2)

ggsave("combined_volcano.pdf", combined_plot, width = 16, height =10)
ggsave("combined_volcano.png", combined_plot, width = 16, height = 10, dpi = 300)


walk2(plot_list, unique(sig_genes$comparison), ~{
  filename <- paste0("volcano_", gsub(" vs ", "_", .y), ".pdf")
  ggsave(filename, .x, width = 6, height = 5)
})



library(dplyr)
library(tidyr)
library(openxlsx)
library(ComplexHeatmap)


group_stats <- sig_genes %>%
  group_by(comparison, direction) %>%
  summarise(
    gene_count = n(),
    avg_log2FC = mean(avg_log2FC),
    max_log2FC = max(abs(avg_log2FC)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = c(gene_count, avg_log2FC, max_log2FC),
    values_fill = 0
  )

write.csv(group_stats, "significant_DE_genes_number_stats.csv")


group_stats <- group_stats %>%
  mutate(
    total_genes = gene_count_Up + gene_count_Down,
    Up_percent = gene_count_Up / total_genes * 100,
    Down_percent = gene_count_Down / total_genes * 100
  ) %>%
  select(
    Comparison = comparison,
    `Total Genes` = total_genes,
    `Upregulated` = gene_count_Up,
    `Downregulated` = gene_count_Down,
    `Up %` = Up_percent,
    `Down %` = Down_percent,
    `Avg Up FC` = avg_log2FC_Up,
    `Avg Down FC` = avg_log2FC_Down,
    `Max FC` = max_log2FC_Up
  )


gene_freq <- sig_genes %>%
  group_by(gene) %>%
  summarise(
    comparison_count = n(),
    comparisons = paste(unique(comparison), collapse = ";"),
    mean_log2FC = mean(avg_log2FC),
    .groups = "drop"
  ) %>%
  arrange(desc(comparison_count), desc(abs(mean_log2FC)))


core_genes <- gene_freq %>%
  filter(comparison_count >= 2) %>%
  mutate(gene_type = case_when(
    mean_log2FC > 1 ~ "Consistent Up",
    mean_log2FC < -1 ~ "Consistent Down",
    TRUE ~ "Variable"
  ))


expr_matrix <- AverageExpression(
  new_seurat,
  assays = "RNA",
  features = core_genes$gene,
  group.by = "new_group"
)$RNA


scaled_matrix <- t(scale(t(expr_matrix)))


heatmap <- Heatmap(
  scaled_matrix,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_km = 3, 
  column_km = 2,
  show_row_names = FALSE,
  row_title = "Gene Clusters",
  column_title = "Expression Trends",
  heatmap_legend_param = list(
    title = "Expression\nZ-score",
    title_position = "leftcenter-rot"
  )
)


pdf("core_genes_heatmap.pdf", width = 10, height = 8)
draw(heatmap)
dev.off()


row_order <- row_order(heatmap)
gene_clusters <- lapply(seq_along(row_order), function(i){
  data.frame(
    Gene = rownames(scaled_matrix)[row_order[[i]]],
    Cluster = paste("Cluster", i),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()


cluster_enrich <- lapply(split(gene_clusters, gene_clusters$Cluster), function(df){
  genes <- df$Gene
  
  # GO
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # KEGG
  ekegg <- enrichKEGG(
    gene = genes,
    organism = "mmu",
    keyType = "kegg"
  )
  
  list(GO = ego, KEGG = ekegg)
})


wb <- createWorkbook()


addWorksheet(wb, "Group Stats")
writeData(wb, "Group Stats", group_stats)


addWorksheet(wb, "Core Genes")
writeData(wb, "Core Genes", core_genes)


for(cl in names(cluster_enrich)){
  sheet_name <- paste("Cluster", cl)
  addWorksheet(wb, sheet_name)
  
  # GO
  go_res <- cluster_enrich[[cl]]$GO@result
  writeData(wb, sheet_name, x = "GO Enrichment", startRow = 1)
  writeData(wb, sheet_name, go_res, startRow = 3)
  
  # KEGG
  kegg_res <- cluster_enrich[[cl]]$KEGG@result
  writeData(wb, sheet_name, x = "KEGG Pathways", startRow = nrow(go_res)+5)
  writeData(wb, sheet_name, kegg_res, startRow = nrow(go_res)+7)
}

# save
saveWorkbook(wb, "DE_genes_statistics.xlsx", overwrite = TRUE)
