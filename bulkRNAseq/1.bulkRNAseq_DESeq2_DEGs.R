#install.packages('BiocManager')  
#BiocManager::install('DESeq2')
library(DESeq2)
library(dplyr)
library(openxlsx)
library(tibble)

input_readcount <- read.xlsx("D://Code//files//bulkRNAseq_readcount.xlsx", sheet = 1, colNames = TRUE, rowNames = FALSE)
colnames(input_readcount)
readcount <- input_readcount %>%
  dplyr::select(c(
    "gene_symbol",
    "read_count_Control_1.", "read_count_Control_2.",
    "read_count_Control_3.", "read_count_Control_4.",
    "read_count_E_407nm_1.", "read_count_E_407nm_2.",
    "read_count_E_407nm_3.", "read_count_E_407nm_4."
  )) %>%
  column_to_rownames("gene_symbol")
colnames(readcount)

#----------------------------------DESeq2----------------------------------

readcount <- readcount[rowSums(readcount == 0) == 0, ]
readcount <-round(readcount)
colnames(readcount)
condition <- factor(c(rep("control",4),rep("EXP",4)))
group<- data.frame(row.names=colnames(readcount), condition)

dds <- DESeqDataSetFromMatrix(countData = readcount, colData = group, design = ~ condition)
head(dds)
dds1 <- DESeq(dds, fitType = 'parametric', minReplicatesForReplace = Inf, parallel = FALSE) 
plotDispEsts(dds1)

res <- results(dds1,pAdjustMethod = "BH",independentFiltering = FALSE)
summary(res) 
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
colnames(res1)
res1 <- data.frame(gene_symbol = rownames(res1), res1, row.names = NULL)
colnames(res1)
#----------------------------------DEGs------------------------------------
res1_up<- res1[which(res1$log2FoldChange > 1 & res1$padj < 0.05),]     
res1_down<- res1[which(res1$log2FoldChange < -1 & res1$padj < 0.05),]    
res1_total <- rbind(res1_up,res1_down)

colnames(res1)
colnames(input_readcount)
merged_df <- merge(res1, input_readcount, by = "gene_symbol", all = FALSE, sort = FALSE) 
colnames(merged_df)
write.xlsx(merged_df, file = "D://Code//files//bulkRNAseq_readcount_DESeq2.xlsx", rowNames = FALSE)

colnames(res1_total)
colnames(input_readcount)
merged_df_total <- merge(res1_total, input_readcount, by = "gene_symbol", all = FALSE, sort = FALSE) 
colnames(merged_df_total)
write.xlsx(merged_df_total, file = "D://Code//files//bulkRNAseq_readcount_DEGs.xlsx", rowNames = FALSE)

