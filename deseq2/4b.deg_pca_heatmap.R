# ===================================================================================
# find most DE genes
# ===================================================================================
## rm(list=ls())

# Load DESeq2 Results
load("3b.dds_rlog.gene.NA.RData")

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

# find most sig DEG genes

ntop=75
alpha_var=0.05 # default

# removing NA values
resValid <- res[!is.na(res$padj),]
resValid$deg <- log10(resValid$pvalue)

# Calculating # of DEG
dim(resValid[resValid$padj<=0.05,])
# 2343

# Upregeulated (in estrogen compared to control)
dim(resValid[resValid$padj<=0.05 & resValid$log2FoldChange>0,])
# 1546

# Downregeulated (in estrogen compared to control)
dim(resValid[resValid$padj<=0.05 & resValid$log2FoldChange<0,])
# 797

# ======================================================================
# VOLCANO PLOT
# ======================================================================

EnhancedVolcano(res,
                lab = res$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                ylab=bquote(~-Log[10]~italic(q-value)),
                pointSize=2,
                labSize=4,
                #DrawConnectors=T,
                colAlpha=.75,
                pCutoff = 0.05,
                FCcutoff=1,
                title="Estrogen Stimulated MCF7 cells vs Control",
                ylim=c(0,270),
                xlim=c(-3,5.3),
                legendPosition = "top",
                caption="",
                subtitle="")

ggplot(res2, aes(x=log2FoldChange, y=(-log10(padj)), label=gene_symbol)) +
  geom_point(data=filter(res2,padj<=0.01), color="red", alpha=0.4)+
  geom_point(data=filter(res2,padj>0.01), color="blue", alpha=0.4)+
  geom_point(data=filter(res2,padj>=0.05), color="black", alpha=0.4)+
  geom_text_repel(data=filter(res2, padj<1e-100 | (log2FoldChange<=-1.3 & padj<0.05) | log2FoldChange>=2 & padj<0.05)) +
  scale_y_continuous(expand = c(0.01, 0.0), trans = "sqrt", limits=c(0,300)) +
  scale_x_continuous(limits=c(-3,5.3), breaks=c(seq(-3,5,1)), expand=c(0.0,0)) +
  #scale_color_manual(name="", labels=c("p>0.05","0.05>p>0.01","p<0.01"), values=c("black","blue","red"))+
  labs(x=expression(log[2]~fold~change), y=bquote(~-Log[10]~italic(q-value)), title="Estrogen Stimulated MCF7 cells vs Control")+
  theme_bw()+theme(legend.position = "right",
                   legend.text=element_text(size=14),
                   legend.title=element_text(size=14),
                   axis.text.x=element_text(color = "black",size=18), 
                   axis.title.y = element_text(face="bold", size=18),
                   axis.title.x = element_text(face="bold", size=18),
                   axis.text.y=element_text(size=20, color = "black"),
                   plot.title = element_text(face="bold", color = "black", size=18, hjust=0.5))


topDEGenes <-  row.names(head(resValid[order(resValid$deg),],ntop))

# add reference genes, if you have a list in ref_genes
#refGenes <- which(rownames(assay(rld)) %in% ref_genes)
#topGenes <- c(refGenes,topVarGenes)
topGenes <- topDEGenes

figureDir="analysis_v1"

# ======================================================================
# build heatmap data
# ======================================================================

mat <- assay(rld)[topGenes, ]
# Getting gene names

s = geneKeyMap[rownames(mat),"symbol"]
rownames(mat) = ifelse(is.na(s),rownames(mat),s)
# normalize rows individually
mat <- mat - rowMeans(mat)
colnames(mat) <- gsub(".Aligned.sortedByCoord.out.bam","",colnames(mat))

# format metadata
df <- data.frame("Condition"=colData(rld)[,"Condition"])
rownames(df) <- rownames(colData(rld))
rownames(df) <- gsub(".Aligned.sortedByCoord.out.bam","",rownames(df))

# PLOT heatmap 
library("pheatmap")

fig_title <- unlist(strsplit(x=c("Heatmap of",ntop,"most significantly differentially expressed genes between control and estrogen stimulated samples"),split=" "))

# OPEN PDF ----------------------------------------
figureFile = gsub(paste(figureDir,"/figures/",paste(fig_title,collapse="_"),".pdf",sep=""),pattern="[|]",replacement = "-")
if(!file.exists(dirname(figureFile))) { dir.create(dirname(figureFile), recursive=TRUE); }
pdf(figureFile, height=9.0, width=7.5)

pheatmap(mat, annotation_col=df,fontsize_row=5, main="Heatmap of 75 most significantly \ndifferentially expressed genes between \ncontrol and estrogen stimulated samples", legend = T, annotation_legend=T, labels_col=c(rld$DNAvarMatch))
#pheatmap(mat, annotation_col=df[1:30,], main=paste(fig_title, collapse=" "), legend = F, annotation_legend=T)

dev.off()
# CLOSE PDF  --------------------------------------

