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
library(igraph)
library(SingleR)
library(clustree)

options(stringsAsFactors = FALSE)
set.seed(123)

PlotColor <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#F1CE63", 
                        "#B383B9", "#8CD17D", "#E15759", "#499894", "#B6992D", 
                        "#86BCB6", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                        "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")
                        
PlotColor_2 <- c("#479D88", "#3C77AF", "#6CB8D2", "#415284", 
                          "#AECDE1", "#EE934E", "#9B5B33", "#E15759",
                          "#B383B9", "#8FA4AE", "#F5D2A8", "#C6307C",
                          "#59A14F", "#B6992D", "#D55640")
                          
sce1 <- readRDS("sce1_with_singleR_annotation_without_48h.rds")
neuron <- readRDS("Neurons_sce_annotated.rds")

Idents(sce1) <- "cell_type"
Idents(neuron) <- "cell_type"

genes_to_check <- c("Dnajb1", "Hsp90ab1", "Hspa5", "Tlr2", "Bag1", "Tlr4",
                    "Hsf1", "Ybx1", "Hsp90aa1", "Gps1", "Timp2", "Alb",
                    "Flvcr1", "Cd163", "Hp", "Lrp1", "Cd44", "Egfr",
                    "Areg", "Erbb2", "Ereg", "Nrd1", "Erbb3", "Egf", "Erbb4",
                    "Il1r1", "Il1r2", "Il1b", "Il1a", "Gnrhr", "Gpr101",
                    "Gapdh", "Ret", "Gfra4", "Gfral")

existing_genes <- genes_to_check[genes_to_check %in% rownames(sce1)]
missing_genes <- genes_to_check[!genes_to_check %in% rownames(sce1)]

desired_order_macro <- unique(sce1@active.ident)
sce1_subset <- subset(sce1, idents = desired_order_macro)
sce1_subset@active.ident <- factor(sce1_subset@active.ident, levels = desired_order_macro)

DotPlot(sce1_subset, features = existing_genes, cols = c("lightgrey", "#E64B35FF"), cluster.idents = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Gene Expression in Macrophages") +
  coord_flip()

ggsave("Multiple_gene_dotplot_macrophage.png", width = 6, height = 8)

existing_genes_neuron <- genes_to_check[genes_to_check %in% rownames(neuron)]
desired_order <- unique(neuron@active.ident)
neuron_subset <- subset(neuron, idents = desired_order)
existing_genes_unique <- make.unique(as.character(existing_genes_neuron))
neuron_subset@active.ident <- factor(neuron_subset@active.ident, levels = desired_order)

DotPlot(neuron_subset, features = existing_genes_unique, cols = c("lightgrey", "#C6307C"),
        cluster.idents = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Gene Expression in Neurons") +
  coord_flip()

ggsave("Multiple_gene_dotplot_neuron.png", width = 5, height = 8)

sce1$original_annotation <- Idents(sce1)
neuron$original_annotation <- Idents(neuron)

check_gene_presence <- function(genes, dataset1, dataset2) {
  results <- data.frame(
    Gene = genes,
    In_sce1 = genes %in% rownames(dataset1),
    In_neuron = genes %in% rownames(dataset2)
  )
  results$In_both <- results$In_sce1 & results$In_neuron
  return(results)
}

gene_report <- check_gene_presence(genes_to_check, sce1, neuron)

common_genes <- intersect(rownames(sce1), rownames(neuron))
sce1 <- sce1[common_genes, ]
neuron <- neuron[common_genes, ]

sce1$batch <- "dataset1"
neuron$batch <- "dataset2"

combined <- merge(sce1, neuron, add.cell.ids = c("ds1", "ds2"))
Idents(combined) <- combined$original_annotation

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, features = rownames(combined))

library(harmony)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "batch", reduction.use = "pca", 
                       dims.use = 1:20, project.dim = FALSE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

Idents(combined) <- "original_annotation"

p1 <- DimPlot(combined, reduction = "umap", group.by = "batch", pt.size = 0.1) +
  ggtitle("Batch correction with Harmony")

p2 <- DimPlot(combined, reduction = "umap", group.by = "original_annotation", 
              pt.size = 0.1, label = TRUE, repel = TRUE) +
  ggtitle("Cell types after integration")

ggsave("batch_correction.png", plot = p1, width = 8, height = 6)
ggsave("cell_types_after_integration.png", plot = p2, width = 10, height = 8)

genes_to_check <- c("Tgfbr2", "Gnrhr", "Il1r2", "Il1r1", "Egfr", 
                    "Cd163", "Lrp1", "Cd36", "Cd40", "Tlr4", "Tlr2")
genes_to_check <- unique(genes_to_check)

existing_genes <- genes_to_check[genes_to_check %in% rownames(combined)]
missing_genes <- genes_to_check[!genes_to_check %in% rownames(combined)]

desired_order_combi <- unique(combined@active.ident)
combi_subset <- subset(combined, idents = desired_order_combi)
combi_subset@active.ident <- factor(combi_subset@active.ident, levels = desired_order_combi)

exp_data <- LayerData(combi_subset, layer = "data")
existing_genes <- intersect(genes_to_check, rownames(exp_data))
gene_subset <- exp_data[existing_genes, ]
avg_exp <- Matrix::rowMeans(gene_subset)

expression_df <- data.frame(
  Gene = names(avg_exp),
  AvgExpression = as.vector(avg_exp)
) %>% 
  arrange(desc(AvgExpression)) %>%
  mutate(Rank = 1:n())

dot_plot <- DotPlot(combi_subset, features = expression_df$Gene, cols = c("lightgrey", "#E64B35FF"),
                    dot.scale = 6, scale = TRUE, cluster.idents = FALSE) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "italic"),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "Gene Expression in Integrated Dataset", y = "Cell Types", x = "Genes") +
  coord_flip()

ggsave("combined_gene_dotplot.png", plot = dot_plot, width = 9, height = 12)

avg_exp <- AverageExpression(combi_subset, assays = "RNA", layer = "data",
                             features = existing_genes, group.by = "cell_type")$RNA
scaled_avg <- t(scale(t(avg_exp)))

gene_cluster <- hclust(dist(scaled_avg), method = "complete")
cell_cluster <- hclust(dist(t(scaled_avg)), method = "complete")

gene_order <- rownames(scaled_avg)[gene_cluster$order]
cell_order <- colnames(scaled_avg)[cell_cluster$order]

combi_subset$celltype_clustered <- factor(as.character(combi_subset$cell_type), levels = cell_order)

dot_plot_clustered <- DotPlot(combi_subset, features = gene_order, cols = c("lightgrey", "#E64B35FF"),
                              dot.scale = 6, scale = TRUE, cluster.idents = FALSE) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, face = "italic"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "Clustered Gene Expression by Cell Type", y = "Cell Types", x = "Genes") +
  coord_flip()

ggsave("clustered_dotplot.png", plot = dot_plot_clustered, width = 10, height = 12)

transposed_avg <- t(scaled_avg)
heatmap_plot <- pheatmap(transposed_avg, cluster_rows = cell_cluster, cluster_cols = gene_cluster,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 7,
                         fontsize_col = 7, treeheight_row = 15, treeheight_col = 20,
                         main = "Gene Expression Clustering", angle_col = 90, silent = FALSE)

ggsave("clustered_heatmap.png", plot = heatmap_plot, width = 7, height = 4)

genes_to_check <- c("Tgfbr2", "Gnrhr", "Il1r2", "Il1r1", "Egfr", "Cd163",
                    "Lrp1", "Cd36", "Cd40", "Tlr4", "Tlr2", "Ybx1", "Dnaja1", "Bag1", 
                    "Timp2", "Gpc1", "Akt1", "Mapk3", "Mapk1", "Ret")
genes_to_check <- unique(genes_to_check)

existing_genes <- genes_to_check[genes_to_check %in% rownames(combined)]
missing_genes <- genes_to_check[!genes_to_check %in% rownames(combined)]

dot_plot <- DotPlot(combined, features = existing_genes, cols = c("lightgrey", "#E64B35FF"),
                    cluster.idents = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  labs(title = "Gene Expression in Integrated Dataset", y = "Cell Types", x = "Genes") +
  coord_flip()

ggsave("combined_gene_dotplot_all.png", plot = dot_plot, width = 10, height = 10)

sessionInfo()