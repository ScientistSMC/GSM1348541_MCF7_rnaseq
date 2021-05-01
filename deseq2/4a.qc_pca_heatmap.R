# ======================================================================
# QC: PCA & heatmap 
#
# CUSTOMIZE:
#    intgroups  = what sampleTable columns define experimental groups
#    ntop       = how many genes to use for qc_heatmap
# input: rld
# ======================================================================
##rm(list=ls())

# Load DESeq2 Results
load("3b.dds_rlog.gene.NA.RData")

# Create Directory for Results
figureDir="analysis_v1"

# Load Libraries
library(DESeq2)
library(genefilter)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Making sure versioned ENSBL IDs are in table for relating gene names to results
geneKeyMap$ENSG_ID <- rownames(geneKeyMap)

# CUSTOMIZE
ntop <- 5000


# Enumerating groups for comparison
intgroups <- c("Condition")

# ======================================================================
# PCA (pdf/html) N most variable
# ======================================================================


fig_title <- unlist(strsplit(x=c("QC PCA",ntop,"most variable",paste0(comment(dds)["by"],"s"),paste0("overlap=",comment(dds)["soMode"])),split=" "))
pcaData <- plotPCA(rld, intgroup = intgroups[1], ntop=ntop, returnData=TRUE)
pcaPercentVar <- round(100 * attr(pcaData , "percentVar"), digits=1)

#  pull in other metadata
pcaData = data.frame(pcaData,sampleTable[rownames(pcaData),])

# format plot
gg=ggplot(pcaData , aes(PC1, PC2, color=Condition, label=Sample)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
  ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
  guides(color=guide_legend(title="Condition")) +
  ggtitle("PCA-5000 Most Variable Genes") +
  theme_bw() +
  theme(plot.title=element_text(size=14, hjust=0.5))
  
  

# OPEN PDF ----------------------------------------
pdfFile = paste(figureDir, "/figures/",paste(fig_title,collapse="_"),".pdf",sep="")
if(!file.exists(dirname(pdfFile))) { dir.create(dirname(pdfFile), recursive=TRUE); }
pdf(pdfFile, height=7.5, width=8.5)
print(gg)
print(gg+geom_label_repel(show.legend=F, point.padding = 0.3))
dev.off()
# CLOSE PDF ---------------------------------------
# OPEN PNG ----------------------------------------
pngFile = paste(figureDir, "/figures/",paste(fig_title,collapse="_"),".png",sep="")
png(pngFile, height=750, width=850, units="px"); print(gg); dev.off()
pngFile = paste(figureDir, "/figures/",paste(fig_title,collapse="_"),".labeled.png",sep="")
png(pngFile, height=750, width=850, units="px"); print(gg+geom_label_repel(show.legend=F)); dev.off()
# CLOSE PDF ---------------------------------------

# ===================================================================================
# heatmap find most variable genes
# ===================================================================================
library("genefilter")
ntop = 100
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),ntop)

#
# build heatmap data
#
mat <- assay(rld)[ topVarGenes, ]
#colnames(mat) = rld$DNAvarMatch
# Not sure why but adding column names results in an error when plotting pheatmap;
# so I had to manually specifiy colnames in pheatmap
s = geneKeyMap[rownames(rld)[topVarGenes],"symbol"]
rownames(mat) = ifelse(is.na(s),rownames(mat),s)
# normalize rows individually
mat <- mat - rowMeans(mat)

# format metadata
df <- data.frame("Condition"=colData(rld)[,intgroups])
rownames(df) <- rownames(colData(rld))
head(df)


# PLOT heatmap 
library("pheatmap")

fig_title <- paste("QC_heatmap_of",ntop,"most_variable_genes", sep="_")

# OPEN PDF ----------------------------------------
figureFile = gsub(paste(figureDir, "/figures/",paste(fig_title,collapse="_"),".pdf",sep=""),pattern="[|]",replacement = "-")
if(!file.exists(dirname(figureFile))) { dir.create(dirname(figureFile), recursive=TRUE); }
pdf(figureFile, height=10, width=5)

pheatmap(mat, annotation_col=df,fontsize_row=5, main=paste("100 Most Variable\n Genes Entire Cohort"), legend = F, annotation_legend=T, labels_col=c(rld$Sample))

dev.off()
# CLOSE PDF  --------------------------------------



