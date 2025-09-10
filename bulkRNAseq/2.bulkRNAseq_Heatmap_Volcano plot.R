#--------------------------------Heatmap-----------------------------
##install.packages("circlize")
library(openxlsx)
library(ComplexHeatmap)  
library(circlize)
library(grid)

data0 <- read.csv("D:\\Code\\files\\bulkRNAseq_Heatmap.csv", header = T,row.names = 1)
data0 <- as.matrix(data0) 
colnames(data0)
data <- t(scale(t(data0)))


mycolors<-colorRamp2(
  c(-2, 0, 2), 
  c("#3882AB", "white", "#BC1B0F")
)
mycolors(seq(-3, 3))


group_annotation <- HeatmapAnnotation(
  Group = factor(c(rep("Control", 4), rep("407nm", 4)), 
                 levels = c("Control", "407nm")),
  col = list(Group = c( "Control" = "#0002FB","407nm" = "#F90005")),
  show_legend = TRUE,
  which = "column",
  annotation_label = "Group", 
  show_annotation_name = FALSE
)


Heatmap(data, 
        name="DEGsHeatmap",
        col= mycolors,
        na_col = "black", 
        rect_gp = gpar(col = NA),
        cluster_rows = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_row_dend = FALSE, 
        cluster_columns = TRUE,
        clustering_distance_columns = "pearson",
        clustering_method_columns = "complete",
        show_column_dend =TRUE,
        column_dend_height = unit(15, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation = group_annotation,
        heatmap_width = unit(7, "cm"),  
        heatmap_height = unit(20, "cm"),
        show_heatmap_legend =TRUE,
        heatmap_legend_param = list(
          title = "z-score",
          direction = "vertical", # c("vertical", "horizontal")
          title_position = "leftcenter-rot",
          title_gp = gpar(fontsize = 12, fontface = "plain", col = "black"),
           legend_width = unit(4, "cm"),
          legend_height = unit(6, "cm")
        )       
        
)
label_grob <- textGrob("n = 1153", x = 0.7, y = 0.7, just = c("right", "top"),
                       gp = gpar(fontsize = 12, fontface = "bold", col = "black"))
grid.draw(label_grob)
#grid.newpage()




#---------------------------Volcano_plot----------------------------

library(dplyr)
data <- read.csv("D:\\Code\\files\\bulkRNAseq_Volcano.csv", header = T,row.names = 1)
colnames(data)
data$log10FDR <- -log10(data$padj)
colnames(data)


data <- data %>% 
  mutate(DEG = case_when(
    log2FoldChange > 1 & padj <= 0.05 ~ "Up",
    log2FoldChange < -1 & padj <= 0.05 ~ "Down",
    TRUE ~ "NoSignificant"  
  ))
table(data$DEG)
data$DEG <- factor(data$DEG, levels = c("Up", "Down","NoSignificant"))
table(data$DEG)


max(data$log2FoldChange)
min(data$log2FoldChange)
max(data$log10FDR, na.rm = TRUE)


library(ggplot2)
library(ggprism)
library(ggrepel)
ggplot(data, aes(x =log2FoldChange, y=log10FDR, colour=DEG)) +
  geom_point(shape = 16,alpha=0.9, size=2) +  
  scale_color_manual(values=c('#EF0032','#060390','gray')) + 
  #ggtitle("Volcano") +
  xlim(c(-6.1, 6.1)) +  
  ylim(c(-1, 90))+ 
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) + 
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  
  labs(x="log2(FoldChange)", y="-log10(FDR)") +  
  theme_prism(border = T)+ 
  theme_bw()+
  theme(
    legend.position="right",
    legend.title = element_blank(),  
    legend.text = element_text(size = 12,color = "black", face = "plain"), 
    axis.title.x = element_text(size = 12, color = "black", face = "bold",margin = margin(t = 7)), 
    axis.title.y = element_text(size = 12, color = "black", face = "bold",margin = margin(r = 10)), 
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),
    plot.margin = margin(20, 20, 20, 20) , 
    panel.border = element_rect(color = "black", linewidth=1.2) 
  )


library(grid)
Up <- textGrob("Up = 531", x = 0.63, y = 0.87, just = c("center", "center"),
               gp = gpar(fontsize = 18, fontface = "bold", col = "#EF0032"))
grid.draw(Up)
Down <- textGrob("Down = 622", x = 0.25, y = 0.87, just = c("center", "center"),
                 gp = gpar(fontsize = 18, fontface = "bold", col = "#060390"))
grid.draw(Down)
#grid.newpage()

