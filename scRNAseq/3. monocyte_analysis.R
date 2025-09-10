library(Seurat)
library(ggplot2)

sce1 <- readRDS("sce1_with_singleR_annotation_without_48h.rds")

monocytes <- subset(sce1, subset = singleR == "Monocytes")

table(monocytes@meta.data[["singleR"]])
table(monocytes@meta.data[["time_point"]])
head(monocytes@meta.data[, c("orig.ident", "singleR")], 5)

table(Idents(monocytes))

monocytes[["percent.mt"]] <- PercentageFeatureSet(monocytes, pattern = "^mt-")
head(monocytes@meta.data)

monocytes <- NormalizeData(monocytes, normalization.method = "LogNormalize", scale.factor = 10000)

monocytes <- FindVariableFeatures(monocytes, selection.method = "vst", nfeatures = 2000,
                                  mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

top20 <- head(VariableFeatures(monocytes), 20)
plot <- VariableFeaturePlot(monocytes)
LabelPoints(plot = plot, points = top20, repel = TRUE)
ggsave("Variable_genes.png", width = 10, height = 4)

monocytes <- ScaleData(monocytes, vars.to.regress = "percent.mt")

monocytes <- RunPCA(monocytes, npcs = 50, verbose = FALSE)
elbow <- ElbowPlot(monocytes, ndims = 20)
ggsave("ElbowPlot_monocytes.png", elbow)

DimHeatmap(monocytes, dims = 1:20, cells = 500, balanced = TRUE)
ggsave("heatmap.png", width = 8, height = 8, dpi = 600)

monocytes <- FindNeighbors(monocytes, dims = 1:15)
monocytes <- FindClusters(monocytes, resolution = 0.2)
monocytes <- RunUMAP(monocytes, dims = 1:15)

DimPlot(monocytes, group.by = "seurat_clusters", split.by = "coexpressed", 
        pt.size = 0.5, reduction = "umap", label = TRUE, ncol = 2) +
  ggtitle("Distribution of Triple+ Cells by Subtype")

ggsave("Monocytes_subtype_umap_15_0.2.png", width = 8, height = 4, dpi = 600)

monocytes$time_point <- factor(monocytes$time_point, levels = c("0h", "4h", "24h"))

DimPlot(monocytes, group.by = "coexpressed", split.by = "time_point", pt.size = 0.5,
        cols = c("Other" = "gray90", "Co-expressed" = "red"), order = "Co-expressed") +
  ggtitle("Triple+ Cells in Monocytes Subtypes (Time Series)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("triple_cells_by_custom_time_order_15_0.2.png", width = 10, height = 4, dpi = 600)

monocytes <- RunTSNE(monocytes, dims = 1:15, perplexity = 30)
DimPlot(monocytes, split.by = "coexpressed", pt.size = 0.5, reduction = "tsne", label = TRUE) +
  ggtitle("t-SNE")

ggsave("Monocytes_subtype_tsne_15_0.2.png", width = 8, height = 4, dpi = 600)

saveRDS(monocytes, file = "Monocytes_without_48h.rds")

sessionInfo()
