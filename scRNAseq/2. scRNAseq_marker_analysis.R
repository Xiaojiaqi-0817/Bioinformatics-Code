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

options(stringsAsFactors = FALSE)
set.seed(123)
PlotColor <- c(  "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#F1CE63", 
                          "#B383B9", "#8CD17D",  "#E15759", "#499894", "#B6992D", 
                          "#86BCB6", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                          "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
)


sce1 <- readRDS("sce1_with_singleR_annotation_without_48h_reunited.rds")


sce1 <- FindVariableFeatures(sce1, selection.method = "vst", nfeatures = 2000,
                             mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

top20 <- head(VariableFeatures(sce1), 20)

plot <- VariableFeaturePlot(sce1)

LabelPoints(plot = plot, points = top20, repel = TRUE)
ggsave("Variable_genes.png", width = 10, height = 4)

sce1 <- ScaleData(sce1, vars.to.regress = "percent.mt")

sce1 <- RunPCA(sce1, npcs = 50, verbose = FALSE)

elbow <- ElbowPlot(sce1, ndims = 30) + geom_vline(xintercept = 20, linetype = 2, color = "red")
ggsave("ElbowPlot.png")

DimHeatmap(sce1, dims = 1:17, cells = 500, balanced = TRUE)
ggsave("heatmap.png", width = 8, height = 8, dpi = 600)

sce1 <- FindNeighbors(sce1, dims = 1:17)
sce1 <- FindClusters(sce1, resolution = 0.2, algorithm = "Louvain")
sce1 <- RunUMAP(sce1, dims = 1:17)

DimPlot(sce1, label = TRUE) + ggtitle("Cell clusters")
ggsave("UMAP_clusters.png", width = 8, height = 6)

sce1 <- RunTSNE(sce1, dims = 1:17, perplexity = 30)
DimPlot(sce1, split.by = "coexpressed", pt.size = 0.5, reduction = "tsne", label = TRUE) +
  ggtitle("t-SNE")

ggsave("sce1_subtype_tsne_17_subset2.png", width = 8, height = 4, dpi = 600)

time_mapping <- list(
  "0h" = c("28", "32", "36", "40", "44", "48"),
  "4h" = c("29", "33", "37", "41", "45", "49"),
  "24h" = c("30", "34", "38", "42", "46", "50")
)

sce1@meta.data$time_mapping <- NA

for (time in names(time_mapping)) {
  sce1@meta.data$time_point[sce1@meta.data[["orig.ident"]] %in% time_mapping[[time]]] <- time
}

table(sce1@meta.data$time_point, useNA = "ifany")
head(sce1@meta.data[c("orig.ident", "time_point")])
table(sce1@meta.data$time_point)

sce1@meta.data$time_point <- factor(sce1@meta.data$time_point, levels = c("0h", "4h", "24h"))
levels(sce1@meta.data$time_point)

expressed_cells <- colnames(sce1)[(
  GetAssayData(sce1, slot = "data")["Siglec1", ] > 0 &
    GetAssayData(sce1, slot = "data")["Itgax", ] > 0 &
    GetAssayData(sce1, slot = "data")["Itgam", ] > 1.5
)]

sce1@meta.data$coexpressed <- ifelse(colnames(sce1) %in% expressed_cells, "Co-expressed", "Other")

table(sce1@meta.data$coexpressed)
table(sce1@meta.data$time_point)
sce1_triple <- subset(sce1, cells = expressed_cells)

celltype_coexp_counts <- sce1@meta.data %>%
  filter(coexpressed == "Co-expressed") %>%
  group_by(seurat_clusters, time_point) %>%
  summarise(coexp_count = n(), .groups = 'drop') %>%
  complete(seurat_clusters, time_point, fill = list(coexp_count = 0))

counts_wide <- celltype_coexp_counts %>%
  pivot_wider(names_from = time_point, values_from = coexp_count, names_sort = TRUE)

triple_sce1_cells <- intersect(colnames(sce1_triple), colnames(sce1))

cell_types <- na.omit(unique(sce1@meta.data[["singleR"]]))
n_types <- length(cell_types)

celltype_colors <- hue_pal(l = 70, c = 100)(n_types)
names(celltype_colors) <- cell_types
background_color <- "#F0F0F0"
  
DimPlot(sce1, group.by = "singleR", cols = celltype_colors, pt.size = 0.5,
        label = TRUE, label.size = 5, repel = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black", family = "Arial"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 9, family = "Arial"),
        legend.title = element_blank(), legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

ggsave("Figure_UMAP_AllTime.pdf", width = 8, height = 6, device = cairo_pdf)

DimPlot(sce1, group.by = "cell_type", split.by = "coexpressed", 
        pt.size = 0.5, reduction = "umap", label = TRUE, ncol = 2) +
  ggtitle("Distribution of Coexpressed Cells by Subtype")

ggsave("sce1_subtype_umap_17_0.2.png", width = 8, height = 4, dpi = 600)

sce1$time_point <- factor(sce1$time_point, levels = c("0h", "4h", "24h"))

DimPlot(sce1, group.by = "coexpressed", split.by = "time_point", pt.size = 0.5,
        cols = c("Other" = "gray90", "Co-expressed" = "red"), order = "Co-expressed") +
  ggtitle("Triple+ Cells in sce1 Subtypes (Time Series)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("triple_cells_by_custom_time_order_15_0.2.png", width = 10, height = 4, dpi = 600)



cluster_markers <- FindAllMarkers(sce1, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25, test.use = "wilcox")

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  filter(p_val_adj < 0.01)

write.csv(top_markers, file = "monocyte_top_markers.csv", row.names = FALSE)
write.csv(cluster_markers, file = "monocyte_cluster_markers.csv", row.names = FALSE)

p1 <- DoHeatmap(sce1, features = top_markers$gene) + NoLegend()
ggsave("heatmap-topgenes10.png", p1, width = 8, height = 15, dpi = 600)

library(GPTCelltype)

markers_list <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10) %>%
  summarise(genes = list(gene)) %>%
  deframe()

Sys.setenv(OPENAI_API_KEY = 'YOUR_API_KEY_HERE')

res_10 <- gptcelltype(cluster_markers, tissuename = "inflamted skin immune cells", 
                      model = 'gpt-4', topgenenumber = 10)
res_20 <- gptcelltype(cluster_markers, tissuename = "inflamted skin immune cells", 
                      model = 'gpt-4', topgenenumber = 20)
res_50 <- gptcelltype(cluster_markers, tissuename = "inflamted skin immune cells", 
                      model = 'gpt-4', topgenenumber = 50)
res_all <- gptcelltype(cluster_markers, tissuename = "inflamted skin immune cells", 
                       model = 'gpt-4')

res_10_2 <- gptcelltype(cluster_markers, tissuename = "inflamted skin tissue", 
                        model = 'gpt-4', topgenenumber = 10)
res_20_2 <- gptcelltype(cluster_markers, tissuename = "inflamted skin tissue", 
                        model = 'gpt-4', topgenenumber = 20)
res_50_2 <- gptcelltype(cluster_markers, tissuename = "inflamted skin tissue", 
                        model = 'gpt-4', topgenenumber = 50)
res_all_2 <- gptcelltype(cluster_markers, tissuename = "inflamted skin tissue", 
                         model = 'gpt-4')

library(tibble)
res_10_df <- enframe(res_10, name = "Cluster", value = "InfSkinImmuneCells_10")
res_20_df <- enframe(res_20, name = "Cluster", value = "InfSkinImmuneCells_20")
res_50_df <- enframe(res_50, name = "Cluster", value = "InfSkinImmuneCells_50")
res_all_df <- enframe(res_all, name = "Cluster", value = "InfSkinImmuneCells_all")

res_10_blood_df <- enframe(res_10_2, name = "Cluster", value = "InfSkin_10")
res_20_blood_df <- enframe(res_20_2, name = "Cluster", value = "InfSkin_20")
res_50_blood_df <- enframe(res_50_2, name = "Cluster", value = "InfSkin_50")
res_all_blood_df <- enframe(res_all_2, name = "Cluster", value = "InfSkin_all")

final_table <- reduce(list(res_10_df, res_20_df, res_50_df, res_all_df,
                           res_10_blood_df, res_20_blood_df, res_50_blood_df, res_all_blood_df),
                      full_join, by = "Cluster")

final_table <- final_table %>% arrange(Cluster)

library(openxlsx)
write.xlsx(final_table, "celltype_comparison_cluster_cell_1.xlsx")


sce1 <- RenameIdents(sce1, res_10)
sce1$Cluster_17_0.5_anno <- Idents(sce1)
sce1$seurat_clusters <- Idents(sce1)

num_tab <- table(Idents(sce1))
freq_tab <- prop.table(num_tab)

bar1 <- barplot(height = freq_tab, width = 1, xlim = c(0, 5), col = c(1:10), 
                legend = rownames(freq_tab), xlab = "")
ggsave("barplot-subtype_freuqtab.png", bar1, width = 8, height = 15, dpi = 600)

DimPlot(sce1, split.by = "coexpressed", pt.size = 0.5, label = TRUE, repel = TRUE,
        order = "Co-expressed") +
  ggtitle("Coexpressed Cells in sce1 Subtypes (Time Series)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("triple_cells_by_annotation_17_0.2.png", width = 10, height = 4, dpi = 600)


nature_umap <- DimPlot(
  sce1,
  # split.by = "coexpressed", 
  pt.size = 0.3,  
  alpha = 0.6,    
  label = TRUE,
  label.size = 4,  
  repel = TRUE,
  ncol = 2,
  cols = PlotColor,  
  raster = FALSE   
) +
  ggtitle("") +    
  theme_cowplot(font_size = 10) +  
  labs(x = "umap 1", y = "umap 2") +  
  theme(
    legend.position = "right",      
    legend.title = element_blank(),  
    legend.key.size = unit(5, "mm"), 
    strip.background = element_blank(),  
    strip.text = element_text(size = 12, face = "bold"),  
    axis.line = element_line(linewidth = 0.5),  
    axis.ticks = element_line(linewidth = 0.5)   
  )

nature_umap


ggsave("Figure4c_CellTypeDistribution.pdf",
       plot = nature_umap,
       width = 13,     
       height = 10,
       units = "cm",
       device = cairo_pdf)  


ggsave("Figure4c_CellTypeDistribution.png",
       plot = nature_umap,
       width = 13,     
       height = 10,
       units = "cm") 




celltype_coexp_counts <- sce1@meta.data %>%
  filter(coexpressed == "Co-expressed") %>%
  group_by(cell_type, time_point) %>%
  summarise(coexp_count = n(), .groups = 'drop') %>%
  complete(cell_type, time_point, fill = list(coexp_count = 0))

counts_wide <- celltype_coexp_counts %>%
  pivot_wider(names_from = time_point, values_from = coexp_count, names_sort = TRUE)

scRNA_triple <- subset(sce1, cells = expressed_cells)


palette <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#F1CE63",
                      "#B383B9","#8CD17D",  "#E15759", "#499894", "#B6992D")
                      

all_cell_types <- levels(factor(celltype_coexp_counts$cell_type))


color_mapping <- setNames(palette[1:length(all_cell_types)], all_cell_types)


plot_data <- celltype_coexp_counts %>%
  mutate(
    time_point = factor(time_point, levels = c("0h", "4h", "24h")),
    cell_type = factor(cell_type, levels = all_cell_types) 
  ) %>%
  filter(coexp_count > 0) %>% 

  mutate(percentage = round(coexp_count/as.numeric(table(sce1@meta.data$time_point)[time_point])*100, 2))


p <- ggplot(plot_data, 
            aes(x = time_point, y = coexp_count, 
                group = cell_type, color = cell_type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Cell Type",
    values = color_mapping,
    breaks = all_cell_types, 
    drop = FALSE 
  ) +

  geom_text(
    aes(label = paste0(coexp_count)),#, "\n(", percentage, "%)")),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  labs(
    title = "Triple Positive Cells Over Time",
    x = "Time Point",
    y = "Cell Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(0, max(plot_data$coexp_count) * 1.3))


ggsave("Fig4d_triple_positive_cells_color_corrected.pdf", p, width = 5, height = 5)
ggsave("Fig4d_triple_positive_cells_color_corrected.png", p, width = 5, height = 5, dpi = 300)


saveRDS(sce1, file = "scRNA-seq_data/sce1_with_singleR_annotation_without_48h.rds")

sessionInfo()
