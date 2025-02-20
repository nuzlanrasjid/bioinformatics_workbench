# script to perform DGE
# setwd("D:/Rstudio gene expression analysis/DESeq2")

#load libraries
library("tidyverse")
library("DESeq2")
library("airway")
library("dplyr")
library("pheatmap")
library("ggplot2")
library("ggrepel")
library("RColorBrewer")
library("apeglm")

#read counts data
count_data <- read.csv("count_data.csv", header = TRUE, row.names = 1)
colnames(count_data)
head(count_data)

#read in column data
coldata <- read.csv("sample_data.csv", header = TRUE, row.names = 1)
colnames(coldata)
head(coldata)

#make sure row names in coldata matches to column names in count_data
all(colnames(count_data) %in% rownames(coldata))

#make sure they in the same order
all(colnames(count_data) == rownames(coldata))

#set factor levels
coldata$Treatment <- factor(coldata$Treatment)
coldata$Sequencing <- factor(coldata$Sequencing)

#set DESeq2 matrix and inport count data and coldata information
dds <- DESeqDataSetFromMatrix(countData = count_data,
                       colData = coldata,
                       design = ~ Sequencing + Treatment)

#setting reference for treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("untreated", "treated"))

#keeping rows
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

#NOTE: collapse technical replicates

#perform statistical analysis of DGE
dds <- DESeq(dds)
res <- results(dds)

#exploring results
summary(res)
res

# adjusting p-value
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#change the result into data frame
res0.01 <- as.data.frame(res0.01)
class(res0.01)
res <- as.data.frame(res)
class(res)

#reorder res0.01 dataset by increasing P value
reorderres0.01 <- res0.01[order(res0.01$pvalue),]
head(reorderres0.01)

# checking specific gene differentially expressed, ex: FBgn0003360
res0.01["FBgn0003360",]
res0.01["FBgn0029167",]

# filter based on p-adjusted value
filtered <- res %>% filter(res$padj < 0.05)

#filter based on fold changes, with a threshold of 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

#save DESeq2 result, both of the original file (res) and filtered one (hits)
write.csv(res, "res.all.csv")
write.csv(filtered, "filtereddata.filter.csv")

#save normalized data
normalized <- counts(dds,normalized=TRUE)
head(normalized)
write.csv(normalized,"normalized_data.csv")

#VISUALIZATION
#dispersion plot
plotDispEsts(dds)

#PCA plot
#variance stabilizing transformation
vds <- vst(dds, blind = FALSE)

#use transformed values to generate PCA plot
plotPCA(vds,intgroup=c("Sequencing", "Treatment"))

#Heatmaps
#R package: pheatmap

#Heatmap of sample-to-sample distance matrix (with clustering) based on the normalized counts.
#1. Generate distant matrix
sampledis <- dist(t(assay(vds)))
sampledismat <- as.matrix(sampledis)

#2. set color scheme
color <- colorRampPalette(rev(brewer.pal(9, "Reds")))(300)

#3.Generate heatmap
pheatmap(sampledismat, clustering_distance_rows = sampledis,
         clustering_distance_cols = sampledis, col=color)


#Heatmap of log transformed normalized counts. We will use the top 10 genes.
#1. ordering the top 10 genes
top10 <- res[order(res$padj),][1:10,]
top10 <- row.names(top10)
top10  

#2. performing log transformation
rld <- rlog(dds,blind = FALSE)

#3. Generate the heatmap
pheatmap(assay(rld)[top10,], cluster_rows = FALSE, show_rownames = TRUE,
         cluster_cols = FALSE)
pheatmap(assay(rld)[top10,], )

#4. adding annotation (optional)
annot <- as.data.frame(colData(dds)[,c('Sequencing', 'Treatment')])
pheatmap(assay(rld)[top10,], cluster_rows = FALSE, show_rownames = TRUE,
         cluster_cols = FALSE, annotation_col = annot)


#Heatmap of Z scores. We will use the top 10 genes
#1. calculate z score
cal_z_score <- function(x) {(x-mean(x)) / sd (x)}
allzscore <- t(apply(normalized, 1, cal_z_score))
# 2. subset top 10 genes
subset_zscore <- allzscore[top10,]
# 3. Generating heatmap
pheatmap(subset_zscore)

#MA plot
plotMA(dds, ylim=c(-2,2))
# removing nose MA plot
resLFC <- lfcShrink(dds, coef = "Treatment_treated_vs_untreated", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))
#change resLFC to dataframe
resLFC <- as.data.frame(resLFC)

#volcano plot
#labelling the genes
resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05] <- "DOWN"
resLFC$delabel <- NA

ggplot(data=resLFC,aes(x=log2FoldChange,y =- log10(pvalue),col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values=c('violet','blue','red'))+
  theme(text=element_text(size=20))
