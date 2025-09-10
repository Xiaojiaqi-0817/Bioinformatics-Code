library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(gtable)


dir_name <- list.files("D:/data/GSE255686")
dir_name
dir_name=list.files("D:/data/GSE255686")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("D:/data/GSE255686", dir_name[i], sep = "/"))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}


for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)]
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)
  scRNAlist[[i]] <- sc
  rm(sc)
}

scRNAlist_filtered <- list()

for (i in seq_along(scRNAlist)) {
  current_scRNA <- scRNAlist[[i]]
  identities <- unique(current_scRNA$orig.ident)
  for (identity in identities) {
    scRNA_subset <- subset(current_scRNA, subset = orig.ident == identity)
    if (identity == "28") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 50 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "29") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 6400 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "30") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 20 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "32") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7300 &
                                 mt_percent < 50 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "33") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "34") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "36") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "37") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7800 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "38") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7800 &
                                 mt_percent < 50 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "40") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7500 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "41") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "42") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7500 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "44") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 25 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "45") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 50 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    } else if (identity == "46") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 30 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "48") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7500 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "49") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7000 &
                                 mt_percent < 40 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }else if (identity == "50") {
      scRNA_filtered <- subset(scRNA_subset,
                               subset = nFeature_RNA > 0 & nFeature_RNA < 7500 &
                                 mt_percent < 30 &
                                 HB_percent < 3 &
                                 nCount_RNA < quantile(nCount_RNA, 0.97) &
                                 nCount_RNA > 1000)
    }

    scRNAlist_filtered[[paste0(identity, "_filtered")]] <- scRNA_filtered
  }
}


scRNAlist_merged_filtered <- Reduce(function(x, y) merge(x, y), scRNAlist_filtered)


print(table(scRNAlist_merged_filtered@meta.data$orig.ident))
scRNAlist <- scRNAlist_merged_filtered


rm(scRNAlist_merged_filtered)
rm(scRNAlist_filtered)
ls()
print(table(scRNAlist@meta.data$orig.ident))
print(table(scRNA_harmony@meta.data$orig.ident))


scRNAlist <- NormalizeData(scRNAlist) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = T)
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
a
ggsave("Batch_before.pdf", plot = a, width = 20, height =15)



top15 <- head(VariableFeatures(scRNAlist), 15)
plot1 <- VariableFeaturePlot(scRNAlist)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3)


feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )


scRNAlist <- JoinLayers(scRNAlist)


g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAlist))


s_genes = cc.genes$s.genes


scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAlist=CellCycleScoring(object = scRNAlist,
                           s.features = s_genes,
                           g2m.features = g2m_genes,
                           set.ident = TRUE)
scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
p4=VlnPlot(scRNAlist, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
            ncol = 2, pt.size = 0.1)
p4
scRNAlist@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()



scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
b
ggsave("Batch_after.pdf", plot = b, width = 20, height =15)

pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated
ggsave("pca_harmony_integrated.pdf", plot = pca_harmony_integrated, width = 20, height =15) 


ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:17) %>% FindClusters(resolution = 1)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:17)

umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
umap_integrated1
umap_integrated2

ggsave("UMAP_orig.ident111.png", umap_integrated1, width = 10, height = 8, dpi = 300)
ggsave("UMAP_cluster111.png", umap_integrated2, width = 10, height = 8, dpi = 300)

getwd()
setwd("D:/results/")
#save(scRNA_harmony,scRNAlist,file = "scdata2.Rdata")
#load("scdata2.Rdata")
table(scRNA_harmony@meta.data$seurat_clusters)



BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
library(celldex)
BiocManager::install("celldex", dependencies = TRUE)
immgen_ref <- celldex::ImmGenData()

sce1 <- scRNA_harmony
sce_for_SingleR <- GetAssayData(sce1, layer = "data")
clusters <- sce1@meta.data$seurat_clusters

pred.immgen <- SingleR(test = sce_for_SingleR,
                       ref = immgen_ref,
                       labels = immgen_ref$label.main,
                       method = "cluster",
                       clusters = clusters,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")


table(pred.immgen$labels)

celltype = data.frame(ClusterID = rownames(pred.immgen),
                      celltype = pred.immgen$labels,
                      stringsAsFactors = FALSE)

sce1@meta.data$singleR = celltype[match(clusters, celltype$ClusterID), 'celltype']

P9 <- DimPlot(sce1, reduction = "umap", group.by = "singleR")
P9
ggsave("UMAP_singleR.png", P9, width = 10, height = 8, dpi = 300)

sce1@meta.data$singleR <- celltype[match(clusters, celltype$ClusterID), 'celltype']

cell_counts <- table(sce1@meta.data$singleR)

print(cell_counts)

cluster_annotation <- data.frame(
  ClusterID = rownames(pred.immgen),  
  SingleR_Label = pred.immgen$labels,
  stringsAsFactors = FALSE
)

print(cluster_annotation)
write.csv(cluster_annotation, "Cluster_SingleR_Annotation.csv", row.names = FALSE)

saveRDS(sce1, file = "sce1_with_singleR_annotation_without_48h.rds")


