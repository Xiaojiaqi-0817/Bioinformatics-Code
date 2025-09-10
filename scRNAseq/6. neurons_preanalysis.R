
library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
##install_github("immunogenomics/harmony")
library(gtable)

Nuerons <- readRDS("D:/data/C57_Raw_counts.RDS")
seurat_obj <- CreateSeuratObject(counts = Nuerons)
print(seurat_obj)
rm(seurat_obj)


patterns <- c("male_C57_Naive_0_rep1", "male_C57_Naive_0_rep2",
              "male_C57_Naive_0_rep3", "male_C57_Naive_0_rep4",
              "male_C57_Naive_0_rep5","female_C57_Naive_0_rep1",
             "female_C57_Naive_0_rep2")

matching_cells <- unlist(lapply(patterns, function(p) grep(p, colnames(seurat_obj), value = TRUE)))
seurat_subset <- subset(seurat_obj, cells = matching_cells)
print(seurat_subset)
head(seurat_subset@meta.data)


gene_names <- rownames(seurat_subset)
mt_genes <- grep("^mt-", gene_names, value = TRUE)

print(mt_genes)


seurat_subset [['percent.mt']]<-PercentageFeatureSet(seurat_subset,pattern = "^mt-")

seurat_subset <- NormalizeData(seurat_subset)

VlnPlot(seurat_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1<-FeatureScatter(seurat_subset,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(seurat_subset,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1
plot2


seurat_subset<-subset(seurat_subset,subset=nFeature_RNA>500 & nFeature_RNA<15000 & percent.mt<10)


seurat_subset<-NormalizeData(seurat_subset,normalization.method = "LogNormalize",scale.factor = 10000)

seurat_subset <- FindVariableFeatures(seurat_subset, selection.method = "vst", nfeatures = 2000)

top10<-head(VariableFeatures(seurat_subset),10)
top10

plot3<-VariableFeaturePlot(seurat_subset)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot3
plot4


all.genes <- rownames(seurat_subset)
seurat_subset <- ScaleData(seurat_subset,features = all.genes)

seurat_subset<-RunPCA(seurat_subset,features = VariableFeatures(object=seurat_subset))
print(seurat_subset[["pca"]],dims = 1:5,nfeatures = 5)

VizDimLoadings(seurat_subset,dims = 1:5,reduction = "pca")
DimPlot(seurat_subset,reduction = "pca")

DimHeatmap(seurat_subset,dims = 1,cells = 500,balanced = TRUE)
DimHeatmap(seurat_subset,dims = 1:10,cells = 500,balanced = TRUE)

seurat_subset<-JackStraw(seurat_subset,num.replicate = 100)
seurat_subset<-ScoreJackStraw(seurat_subset,dims = 1:20)
JackStrawPlot(seurat_subset, dims = 1:20)


seurat_subset<-FindNeighbors(seurat_subset,dims = 1:20)%>% FindClusters(resolution = 1)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:20)
umap_integrated1 <-DimPlot(seurat_subset, reduction = "umap")
umap_integrated1
ggsave("UMAP_cluster.png", umap_integrated1, width = 10, height = 8, dpi = 300)
head(seurat_subset)


a<- VlnPlot(seurat_subset,features = c("Rbfox3"))
b<- VlnPlot(seurat_subset,features = c("Sparc"))
c<- a+b
ggsave("differ_Nurons_not.png", c, width = 15, height = 8, dpi = 300)


seurat_subset@meta.data[["seurat_clusters"]] <- as.character(seurat_subset@meta.data[["seurat_clusters"]])
clusters_to_extract <- c("1", "2", "3", "6", "7", "8", "9", "10", "11", "13", "15", "16", "17")
Neurons_sce <- subset(seurat_subset, idents = clusters_to_extract)
table(Neurons_sce@meta.data[["seurat_clusters"]])
Neurons_sce <- NormalizeData(Neurons_sce) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

ElbowPlot(Neurons_sce, ndims=50, reduction="pca")
Neurons_sce <- FindNeighbors(Neurons_sce, reduction = "pca", dims = 1:12)
Neurons_sce <- FindClusters(Neurons_sce, resolution = 1)
Neurons_sce@meta.data$seurat_clusters <- Neurons_sce@meta.data$seurat_clusters
table(Neurons_sce@meta.data$seurat_clusters)
Neurons_sce <- RunUMAP(Neurons_sce, reduction = "pca", dims = 1:12)
p1 <- DimPlot(Neurons_sce, reduction = "umap", label = TRUE, raster = FALSE)
p1
ggsave("Neurons_plot1.png", p1, width = 8, height = 6, dpi = 300)
head(Neurons_sce@meta.data)


Neurons_sce <- readRDS("D:/results/Neurons/Neurons_sce_annotated.RDS")


features <- c("Tac1","Gpx3","Hpca","Mrgprd","Sst","Nefh","Pvalb","Cadps2","Fam19a4","Th",
              "Htr3a","Cplx2","Nptx1","Hapln4","Pvalb")
dot1<- FeaturePlot(
  Neurons_sce,
  features = features
)
dot1
ggsave("Neurons_dot.png", dot1, width = 25, height = 20, dpi = 300)

dot2<- FeaturePlot(
  Neurons_sce,
  features = c("Rbfox3","Sparc")
)
dot2

receptor <- c("Egfr","Tgfbr2")

dot3<- FeaturePlot(
  Neurons_sce,
  features = receptor
)
dot3
ggsave("Neurons_receptor_dot.png", dot3, width = 15, height = 6, dpi = 300)


x <- c("Tac1", "Gpx3",
       "Tac1", "Hpca",
       "Mrgprd",
       "Sst",
       "Nefh",
       "Pvalb",
       "Cadps2",
       "Fam19a4", "Th",
       "Fam19a4")
p2 <- VlnPlot(
  Neurons_sce,
  features = features,  
  pt.size = 0,  
  group.by = "seurat_clusters",  
  stack = TRUE,  
  same.y.lims = TRUE, 
  flip = FALSE  
) + NoLegend()  
p2
ggsave("Neurons_vlnPlot.png", p2, width = 12, height = 8, dpi = 300)


Neurons_sce@meta.data$cell_type <- "Othercells"

Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters %in% c(1, 6, 11)] <- "NP"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters %in% c(3, 7, 15)] <- "PEP1"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters == 2] <- "PEP2"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters %in% c(12,18)] <- "SST"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters == 0] <- "NF1"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters %in% c(5,9,17)] <- "NF2"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters == 14] <- "NF3"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters %in% c(8,16)] <- "cLTMR1"
Neurons_sce@meta.data$cell_type[Neurons_sce@meta.data$seurat_clusters == 10] <- "p_cLTMR2"


print(table(Neurons_sce@meta.data$cell_type))
#   NP      Othercells    PEP1        PEP2
#    1268       3848       1066        512
#   cLTMR1        NF1        NF2        NF3        NP    Othercells   p_cLTMR2       PEP1       PEP2        SST
#      537        733        955        190       1268        696        337       1066        512        400

p_umap_celltype <- DimPlot(
  Neurons_sce,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
)
p_umap_celltype

ggsave("UMAP_celltype_updated.png", p_umap_celltype, width = 10, height = 8, dpi = 300)
saveRDS(Neurons_sce, file = "Neurons_sce_annotated.rds")



selected_cell_types <- c("NP", "PEP1", "PEP2")
Neurons_selected <- subset(Neurons_sce, subset = cell_type %in% selected_cell_types)
print(table(Neurons_selected@meta.data$cell_type))

saveRDS(Neurons_selected, file = "Neurons_selected.rds")

count_matrix <- GetAssayData(Neurons_selected, assay = "RNA", layer = "counts")
meta_data <- Neurons_selected@meta.data
saveRDS(count_matrix, "Neurons_selected_count_matrix.rds")
saveRDS(meta_data, "Neurons_selected_meta_data.rds")


Neurons_selected<- NormalizeData(Neurons_selected) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)


ElbowPlot(Neurons_selected, ndims=50, reduction="pca")
Neurons_selected <- FindNeighbors(Neurons_selected, reduction = "pca", dims = 1:11)
Neurons_selected <- FindClusters(Neurons_selected, resolution = 1)
Neurons_selected@meta.data$seurat_clusters <- Neurons_selected@meta.data$seurat_clusters
table(Neurons_selected@meta.data$seurat_clusters)
Neurons_selected <- RunUMAP(Neurons_selected, reduction = "pca", dims = 1:11)

p3 <- DimPlot(Neurons_selected, reduction = "umap", label = TRUE, raster = FALSE, group.by = "cell_type")
p3
ggsave("NP_PEP1_PEP2.png", p3, width = 10, height = 8, dpi = 300)
head(Neurons_selected@meta.data)

getwd()
setwd("D:/results/Neurons")


dir.create("NP_PEP1_PEP2_plots", showWarnings = FALSE)
generate_and_save_plot <- function(gene, type) {
  if (!gene %in% rownames(GetAssayData(Neurons_selected, layer = "data"))) {
    cat("Not find", gene, "genes\n")
    return()
  }
  gene_expression <- GetAssayData(Neurons_selected, layer = "data")[gene, ]
  if (all(gene_expression == 0)) {
    cat("gene", gene, "exoression level is 0\n")
    return()
  }
  feature_plot <- FeaturePlot(Neurons_selected,
                              features = gene,
                              reduction = "umap",
                              cols = c("lightgrey", "red"))
  print(feature_plot)

  filename <- file.path("NP_PEP1_PEP2_plots", paste0(type, "_", gene, ".png"))
  ggsave(filename,
         plot = feature_plot,
         width = 8,
         height = 6,
         dpi = 300,
         bg = "white")
}


markers <- c("Twsg1","Hamp2","Hamp","Ret","Gfra4","Tgfbr2","Gfral")
for (gene in markers) {
  generate_and_save_plot(gene, "Gdf15")
}


#GnRH1 no Kiss1
markers <- c("Kiss1","Oxt","Gnrhr","Kiss1r","Mapk3","Mapk1","Tac2")
for (gene in markers) {
  generate_and_save_plot(gene, "GnRH1")
}

#HBEGF
markers <- c("Mmp9","Mmp17","Mmp14")
for (gene in markers) {
  generate_and_save_plot(gene, "HBEGF")
}

#Hemopexin
markers <- c("Plg","Timp2","Lrp1")
for (gene in markers) {
  generate_and_save_plot(gene, "Hemopexin")
}

#HSPA1A
markers <- c("Tlr2","Tlr4","Ybx1","Hspa1b")
for (gene in markers) {
  generate_and_save_plot(gene, "HSPA1A")
}

#IL1R2
markers <- c("Tnf","Il6ra","Il1rap","Il1rn")
for (gene in markers) {
  generate_and_save_plot(gene, "IL1R2")
}

#IL1RN
markers <- c("Tnf","Il6","Il6ra")
for (gene in markers) {
  generate_and_save_plot(gene, "IL1RN")
}

