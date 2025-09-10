library(clusterProfiler)
library(stats)
library(fgsea)
library(enrichplot) 
library(ggplot2) 
library(grid)

data0 <- read.csv("D:\\Code\\files\\bulkRNAseq_GSEA.csv", header = T,row.names = 1)
colnames(data0)

control_group <- data0[, 1:4]  
experiment_group <- data0[, 5:8] 
experiment_mean <- rowMeans(experiment_group)
control_mean <- rowMeans(control_group)
experiment_sd <- apply(experiment_group, 1, sd)
control_sd <- apply(control_group, 1, sd)
experiment_sd <- pmax(experiment_sd, 0.2 * experiment_mean)
control_sd <- pmax(control_sd, 0.2 * control_mean)
S2N <- (experiment_mean - control_mean) / (experiment_sd + control_sd)
data0$S2N <- S2N

geneSet_go <- read.gmt("D:\\Code\\files\\m5.go.bp.v2024.1.Mm.symbols.gmt")
colnames(geneSet_go)
subset(geneSet_go, term == 	"GOBP_CIRCADIAN_RHYTHM")
sum(geneSet_go$term == "GOBP_CIRCADIAN_RHYTHM")


geneList <-data0$S2N                 
names(geneList) <- rownames(data0)      
geneList <- sort(geneList, decreasing = T)  
head(geneList)
set.seed(123) 
GSEA_enrichment <- GSEA(geneList,                 
                        TERM2GENE = geneSet_go, 
                        pvalueCutoff = 0.05,      
                        minGSSize = 15,           
                        maxGSSize = 500,         
                        eps = 0,                  
                        pAdjustMethod = "BH")     

result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)
colnames(result)
print(subset(result, ID == 	"GOBP_CIRCADIAN_RHYTHM"))


gseaplot2(GSEA_enrichment, "GOBP_CIRCADIAN_RHYTHM", color = "green3",
          base_size = 15,
          rel_heights = c(1.5, 0.5, 0.5),
          subplots = 1:3,
          ES_geom = "line",
          pvalue_table = F) 


NES <- textGrob("NES = 1.607", x = 0.65, y = 0.9, just = c("left", "top"),
                       gp = gpar(fontsize = 26, fontface = "plain", col = "black"))
grid.draw(NES)

p.adjust <- textGrob("p.adjust = 0.013", x = 0.65, y = 0.82, just = c("left", "top"),
                   gp = gpar(fontsize = 26, fontface = "plain", col = "black"))
grid.draw(p.adjust)

#grid.newpage()





