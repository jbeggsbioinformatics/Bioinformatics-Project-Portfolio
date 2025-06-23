load("miRNA_DESeq2_analysis.Rdata")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force=TRUE)

library(DESeq2)

BiocManager::install("RUVSeq", force=TRUE)

library(RUVSeq)

BiocManager::install("tidyverse", force=TRUE)

library(tidyverse)

BiocManager::install("pheatmap", force=TRUE)

library(pheatmap)

BiocManager::install("ggplot2", force=TRUE)

library(ggplot2)

BiocManager::install("RColorBrewer", force=TRUE)

library(RColorBrewer)

BiocManager::install("repr", force=TRUE)

library(repr)

BiocManager::install("ashr", force=TRUE)

library(ashr)


BiocManager::install("gridExtra", force=TRUE)

library(gridExtra)

BiocManager::install("rlang", force=TRUE)

library(rlang)

BiocManager::install("cowplot", force=TRUE)

library(cowplot)

BiocManager::install("reshape2", force=TRUE)

library(reshape2)

BiocManager::install("ggrepel", force=TRUE)

library(ggrepel)

BiocManager::install("ggplot2", force=TRUE)

library(ggplot2)

BiocManager::install("reshape", force=TRUE)

library(reshape)

setwd("~/Dropbox/miRNA")

CBF_PCA<-function(data, groups, useLabels=F, labels = "", pcs=c(1,2), type='scores', scale=T, legendName="Treatment"){
  
  #this vector calculates the colours needed based on the number of groups
  #colores<-rainbow(length(unique(groups)))
  colores = brewer.pal(6,"Set1")
  
  # INPUTS:
  #
  #  data - data.frame or matrix
  #   - data to analyse with variables in columns and samples in rows
  #  groups - factor
  #   - a grouping factor that determines the colours of the points in scores plots
  #  useLabels (optional) - boolean - default=FALSE
  #   - if TRUE will draw labels next to each point
  #  labels (optional) - default=""
  #   - labels to be drawn next to each point. If useLabels=TRUE and labels is empty will use rownames(data) as labels.
  #  pcs (optional) - a vector of 2 integers - default=c(1,2)
  #   - principal components to be plotted
  #  type (optional) - string - default='scores'
  #   - which plot to produce at the end. Possible: 'scores', 'loadings', 'varAcc'.
  #  scale (optional) - boolean - default=TRUE
  #   - determines if variables are to be scaled (recommended)
  #  legendName (optional) - string - default='Groups'
  #   - sets the name for the legend of colours (set by groups)
  #
  # OUTPUTS:
  # a ggplot object. If not assigned to a variable will show the plot.
  
  
  if(scale){
    pc<-prcomp(data, scale = T)
  } else {
    pc <- prcomp(data)
  }
  
  
  if(type=='scores'){
    if(useLabels & length(labels) != nrow(data)){
      print("Warning: The labels not given or given incorrectly. Using rownames.")
      labels <- rownames(data)
    }
    
    pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
    
    if(useLabels) pcdf$labels<-labels
    
    perc_accounted<-summary(pc)$importance[2,pcs]*100
    
    .e <- environment()
    p <- ggplot(data=pcdf, aes(x=pc1, y=pc2), environment=.e) + 
      geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
      scale_fill_manual(values=colores,name=legendName)
    
    if(useLabels)  p <- p + geom_text_repel(aes(label = labels))
    
    p <- p+ 
      xlab(paste("PC",pcs[1], " (", round(perc_accounted[1],2), "%)", sep=""))+
      ylab(paste("PC",pcs[2], " (", round(perc_accounted[2],2), "%)", sep=""))+
      ggtitle("PCA plot")+
      theme_bw(base_size=20)+
      theme(legend.position="bottom")
    
    p
    
  } else if(type=='loadings'){
    
    if(useLabels & length(labels) != nrow(pc$rotation)){
      print("Warning: loadings labels not given or given incorrectly. Using the column names.")
      labels <- colnames(data)
    }
    
    pcdf<-data.frame(load1=pc$rotation[,pcs[1]], load2=pc$rotation[,pcs[2]], var=labels)
    
    label_offset_x <- 0.035 * (range(pcdf$load1)[2] - range(pcdf$load1)[1])
    label_offset_y <- 0.035 * (range(pcdf$load2)[2] - range(pcdf$load2)[1])
    
    .e <- environment()
    
    p <- ggplot(data=pcdf, aes(x=load1, y=load2), environment=.e) + geom_point()
    
    if(useLabels) p <- p + geom_text_repel(aes(x=load1,y=load2),label=labels)
    
    
    p <- p+
      xlab(paste("Loadings for PC",pcs[1],sep=""))+
      ylab(paste("Loadings for PC",pcs[2],sep=""))+
      ggtitle("PCA loadings plot")+
      theme_bw(base_size=20)
    p
    
  } else if(type=='varAcc'){
    perc_accounted <- (pc$sdev/sum(pc$sdev)*100)
    perc_with_cumsum <- data.frame(pc = as.factor(1:length(perc_accounted)),
                                   perc_acc = perc_accounted,
                                   perc_cumsum = cumsum(perc_accounted))
    p<- ggplot(data = perc_with_cumsum, aes(x=pc, y=perc_cumsum))+
      geom_bar(stat='identity', col='black', fill='white')+
      geom_hline(yintercept = 95, col='red')+
      geom_hline(yintercept = 0, col='black')+
      xlab('PC')+
      ylab('% Variance')+
      ggtitle('% Variance accounted for by principle components')+
      theme_bw()
    print(p)
    
  }else {
    cat(sprintf("\nError: no type %s", type))
  }
  
  return(p)
}

counts = read.csv("miR.Counts.csv", header=T, stringsAsFactors = F)
#convert to matrix with gene ID as rownames
count_matrix = as.matrix(counts[,2:12])
rownames(count_matrix) = counts$miRNA

# Read in sample metadata
meta = read.csv("metadata.csv", header=T, stringsAsFactors = F)
#make sampleID the rownames
rownames(meta) = meta$Sample

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta,
                              design = ~Group)

dds$Group <- factor(dds$Group, levels = c("HC","RA"))

# filter out any genes that are not expressed i.e. count of 0 in all samples
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds)

res <- results(dds, independentFiltering=FALSE)


#view the first few results (res) in adj p-val order (most significant at the top)
head(res[order(res$padj),])

write.csv(res[order(res$padj),], "RAvHC_results.csv")

CBF_PCA(t(normalized_counts), meta[,2], useLabels=F, pcs=c(1,2), legendName = "Group")

write.csv(normalized_counts, "RAvHC_normalisedcounts.csv")

sig_genes = rownames(res)[which(res$padj < 0.05 )]
sig_counts = normalized_counts[rownames(normalized_counts) %in% sig_genes, ]

#create log counts to make visualisation more sensible
log_sig_counts = log2(sig_counts + 1)

#make the headings of log_sig_counts just the HC or RA number
colnames(log_sig_counts) = c("HC20","RA018","RA029","RA034","HC27",
                             "HC28","HC31","HC32","RA004","RA008","RA017")

#plot the heatmap
pheatmap(
  mat = log_sig_counts,
  scale = "row",  
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 0,
  border_color = NA,
  height = height,
  fontsize_row = 10,
  fontsize_col = 10
)

results_df = as.data.frame(res[order(res$padj),])

#remove rows with NA (this is as a result of extreme outliers in a handful of genes)
results_df = na.omit(results_df)

#add a new column 
results_df$significant <- "NO"

# if Adj pvalue < 0.05, set as "YES" 
results_df$significant[results_df$padj < 0.05] <- "YES"

#add -log10 padj for ease of charting
results_df$logpadj = -log10(results_df$padj)

ggplot(results_df, aes(x = log2FoldChange, y = logpadj)) +
  geom_point(aes(color = significant), size=3) +
  scale_color_manual(values = c("NO" = "darkgrey", "YES" = "red")) +
  labs(title = paste0("Volcano Plot: RA v HC"), x = "logFC", y = "-log10(p-value)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12))

library(reshape)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)

#convert log counts to dataframe
log_sig_counts_df = as.data.frame(log_sig_counts)
log_sig_counts_df$GeneID = rownames(log_sig_counts)


# melt data into a long format:
miR_melted <- melt(log_sig_counts_df)

# assign some sensible column names:
colnames(miR_melted) <- c("miR","Group","expression")

# make the replicates all have the same sample name:
miR_melted$Group <- gsub("[0-9]","", miR_melted$Group)


#boxplot function
baw_plot <- function(d){
  ggplot(d, aes(Group, expression, fill = Group)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    facet_wrap(~miR) +
    theme_bw() +
    ylab("expression (log2)") +
    theme(axis.text.x = element_text(size = 10, face="bold"),
          axis.text.y = element_text(size = 10, face="bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=10, face = "bold"),
          legend.position = "none",
          legend.title = element_blank())+
    geom_jitter(width=0.2, shape=16, size=1.5)
}


# lay all the plots out in grid form, with each page having 3 x 4 plots (change accordingly):
b <- marrangeGrob(a, nrow = 4, ncol = 3)

# write them out to a .pdf file:
ggsave("miR-sig-boxplots.pdf", b, width = 7, height = 10)

save.image("miR-sig-boxplots.pdf")

load("miR-sig-boxplots.pdf")

#extract rows for just the genes of interest with names these can be any names these are just an example
my.genes = as.data.frame (c("hsa-miR-151a-3p","hsa-miR-543","hsa-miR-126-3p","hsa-miR-409-3p"))
my.gene.exp = as.data.frame(normalized_counts[rownames(normalized_counts) %in% my.genes[,1],])
my.gene.exp$GeneID = rownames(my.gene.exp)
colnames(my.gene.exp) = c("HC20","RA018","RA029","RA034","HC27","HC28",
                          "HC31","HC32","RA004","RA008","RA017","GeneID")

# melt data into a long format:
miR_melted <- melt(my.gene.exp)

# assign some sensible column names:
colnames(miR_melted) <- c("miR","Group","expression")

# make the replicates all have the same sample name:
miR_melted$Group <- gsub("[0-9]","", miR_melted$Group) 

#boxplot function
baw_plot <- function(d){
  ggplot(d, aes(Group, expression, fill = Group)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    facet_wrap(~miR) +
    theme_bw() +
    ylab("expression (log2)") +
    theme(axis.text.x = element_text(size = 10, face="bold"),
          axis.text.y = element_text(size = 10, face="bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=10, face = "bold"),
          legend.position = "none",
          legend.title = element_blank())+
    geom_jitter(width=0.2, shape=16, size=1.5)
}

#apply the function to the melted data
a <- dlply(miR_melted, "miR", baw_plot)

# lay all the plots out in grid form, with each page having 3 x 4 plots:
b <- marrangeGrob(a, nrow = 4, ncol = 3)

# write them out to a .pdf file:
ggsave("filename.pdf", b, width = 7, height = 10)

facet_wrap(~miR)

save.image("filename.pdf")


sessionInfo()

#extract rows for just the genes of interest with names these can be any names these are just an example
my.genes = as.data.frame (c("hsa-miR-126-3p","hsa-miR-543","hsa-miR-151a-3p","hsa-miR-411-5p", "hsa-miR-493-3p", "hsa-miR-144-3p", "hsa-miR-126", "hsa-miR-31-5p", "hsa-miR-493-5p", "hsa-miR-146a-5p", "hsa-miR-224-5p", "hsa-miR-654-3p", "hsa-miR-126-5p"))
my.gene.exp = as.data.frame(normalized_counts[rownames(normalized_counts) %in% my.genes[,1],])
my.gene.exp$GeneID = rownames(my.gene.exp)
colnames(my.gene.exp) = c("HC20","RA018","RA029","RA034","HC27","HC28",
                          "HC31","HC32","RA004","RA008","RA017","GeneID")

# melt data into a long format:
miR_melted <- melt(my.gene.exp)

# assign some sensible column names:
colnames(miR_melted) <- c("miR","Group","expression")

# make the replicates all have the same sample name:
miR_melted$Group <- gsub("[0-9]","", miR_melted$Group) 

#boxplot function
baw_plot <- function(d){
  ggplot(d, aes(Group, expression, fill = Group)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    facet_wrap(~miR) +
    theme_bw() +
    ylab("expression (log2)") +
    theme(axis.text.x = element_text(size = 10, face="bold"),
          axis.text.y = element_text(size = 10, face="bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=10, face = "bold"),
          legend.position = "none",
          legend.title = element_blank())+
    geom_jitter(width=0.2, shape=16, size=1.5)
}

#apply the function to the melted data
a <- dlply(miR_melted, "miR", baw_plot)

# lay all the plots out in grid form, with each page having 3 x 4 plots:
b <- marrangeGrob(a, nrow = 4, ncol = 3)

# write them out to a .pdf file:
ggsave("Identified-miRNA.pdf", b, width = 7, height = 10)

facet_wrap(~miR)

setwd("~/Dropbox/miRNA")

read.csv("ComBatRAvsHCallsamplesAlltoptable.csv", header=T)

read.csv("combat_norm_cpm.csv", header=T)

data <- read.csv("combat_norm_cpm.csv", header = TRUE)

row.names(data) <- data[, 1]

str(data)

summary(data)


head(data[order(data$padj),])

df <- read.csv("combat_norm_cpm.csv")

names(df)[1] <- "GeneID"

# melt data into a long format:
mRNA_melted <- melt(df)

# assign some sensible column names:
colnames(mRNA_melted) <- c("mRNA","Group","expression")

# make the replicates all have the same sample name:
mRNA_melted$Group <- gsub("[0-9]","", mRNA_melted$Group) 

baw_plot <- function(d){
  ggplot(d, aes(Group, expression, fill = Group)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    facet_wrap(~df) +
    theme_bw() +
    ylab("expression (log2)") +
    theme(axis.text.x = element_text(size = 10, face="bold"),
          axis.text.y = element_text(size = 10, face="bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=10, face = "bold"),
          legend.position = "none",
          legend.title = element_blank())+
    geom_jitter(width=0.2, shape=16, size=1.5)
}

#apply the function to the melted data
a <- dlply(mRNA_melted, "mRNA", baw_plot)

# lay all the plots out in grid form, with each page having 3 x 4 plots (change accordingly):
b <- marrangeGrob(a, nrow = 4, ncol = 3)

# write them out to a .pdf file:
ggsave("mRNABox.pdf", b, width = 7, height = 10)
