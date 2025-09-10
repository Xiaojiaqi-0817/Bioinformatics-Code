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
library(GPTCelltypeSXY)
library(apiSXY)
library(openxlsx)

monocytes <- readRDS("Monocytes_without_48h.rds")


cluster_markers <- FindAllMarkers(
  monocytes,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox")

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  filter(p_val_adj < 0.01)

write.csv(top_markers, file = "monocyte_cluster_markers.csv", row.names = FALSE)

DoHeatmap(monocytes, features = top_markers$gene) + NoLegend()
ggsave("heatmap-topgenes10.png", width = 8, height = 15, dpi = 600)


Sys.setenv(OPENAI_API_KEY = 'YOUR_API_KEY_HERE')

res_10 <- gptcelltype(cluster_markers, 
                      tissuename = "monocytes in inflammated skin tissue", 
                      model = 'gpt-4',
                      topgenenumber = 10) 
res_20 <- gptcelltype(cluster_markers, 
                      tissuename = "monocytes in inflammated skin tissue", 
                      model = 'gpt-4',
                      topgenenumber = 20) 
res_50 <- gptcelltype(cluster_markers, 
                      tissuename = "monocytes in inflammated skin tissue", 
                      model = 'gpt-4',
                      topgenenumber = 50) 
res_all <- gptcelltype(cluster_markers, 
                       tissuename = "monocytes in inflammated skin tissue", 
                       model = 'gpt-4') 
res_10_1 <- gptcelltype(cluster_markers, 
                        tissuename = "inflamted skin immune cells", 
                        model = 'gpt-4',
                        topgenenumber = 10) 
res_20_1 <- gptcelltype(cluster_markers, 
                        tissuename = "inflamted skin immune cells", 
                        model = 'gpt-4',
                        topgenenumber = 20) 
res_50_1 <- gptcelltype(cluster_markers, 
                        tissuename = "inflamted skin immune cells", 
                        model = 'gpt-4',
                        topgenenumber = 50) 
res_all_1 <- gptcelltype(cluster_markers, 
                         tissuename = "inflamted skin immune cells", 
                         model = 'gpt-4') 

res_10_df <- enframe(res_10, name = "Cluster", value = "Cskin_10")
res_20_df <- enframe(res_20, name = "Cluster", value = "_skin_20")
res_50_df <- enframe(res_50, name = "Cluster", value = "skin_50")
res_all_df <- enframe(res_all, name = "Cluster", value = "skin_all")

res_10_1_df <- enframe(res_10_1, name = "Cluster", value = "Cell_10")
res_20_1_df <- enframe(res_20_1, name = "Cluster", value = "Cell_20")
res_50_1_df <- enframe(res_50_1, name = "Cluster", value = "Cell_50")
res_all_1_df <- enframe(res_all_1, name = "Cluster", value = "Cell_all")

final_table <- reduce(
  list(
    res_10_df, res_20_df, res_50_df, res_all_df,
    res_10_1_df, res_20_1_df, res_50_1_df, res_all_1_df
  ),
  full_join,
  by = "Cluster"
)

final_table <- final_table %>% arrange(Cluster)

write.xlsx(final_table, "celltype_comparison_cluster_cell_1.xlsx")


monocytes <- RenameIdents(monocytes, res_all)
table(Idents(monocytes))

DimPlot(monocytes, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("cluster_cell_1.png", width = 5, height = 5, dpi = 600)


DimPlot(monocytes, group.by = "subtype", split.by = "coexpressed", 
        pt.size = 0.5, reduction = "umap", label = TRUE, ncol = 2) +
  ggtitle("Distribution of Coexpressed Cells by Subtype")

ggsave("Monocytes_subtype_umap_subset2_15_0.2.png", width = 8, height = 4, dpi = 600)

monocytes$time_point <- factor(monocytes$time_point, levels = c("0h", "4h", "24h"))

DimPlot(monocytes, group.by = "coexpressed", split.by = "time_point", pt.size = 0.5,
        cols = c("Other" = "gray90", "Co-expressed" = "red"), order = "Co-expressed") +
  ggtitle("Coexpressed Cells in Monocytes Subtypes (Time Series)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("triple_cells_by_custom_time_order_15_0.2_subset2.png", width = 10, height = 4, dpi = 600)


DimPlot(monocytes, split.by = "coexpressed", pt.size = 0.5, reduction = "tsne", label = TRUE) +
  ggtitle("t-SNE")

ggsave("Monocytes_subtype_tsne_15_0.2_subset2.png", width = 8, height = 4, dpi = 600)

saveRDS(monocytes, file = "Monocytes_without_48h.rds")

sessionInfo()
