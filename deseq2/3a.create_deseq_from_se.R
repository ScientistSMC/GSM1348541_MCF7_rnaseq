# 3.create_deseq_from_se.R


load("2z.se.gene.IntersectionStrict.RData")
library(DESeq2)


#----------------------------------------------------------------------
#
# create DESeq2 object
#
#----------------------------------------------------------------------
# http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

print(paste("start:", date()))

# -----------------------------------------------------
# DDS
#
# counts from a SE object
# -----------------------------------------------------
figureDir = "./figures_v1"  # match comment(se)["sampleFile"]

##### Estrogen Treatment vs Control #####
dds = DESeqDataSet(se, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >= 4, ]  # DESeq2 1.3.6 Pre-Filtering; remove rows with low expression
dds <- DESeq(dds)  # 2-3 minutes to build model
comment(dds) = c(rowSumMin = "4", comment(se))
print(paste("done dds:", date()))

res <- results(dds, contrast=c("Condition","Estrogen","Control"))

# Adding gene symbols/names to results
s = geneKeyMap[rownames(res),"symbol"]
res$gene_symbol <- ifelse(is.na(s),rownames(res),s)

write.csv(res, "./results/2021-04-29_res.estrogen_vs_control.csv")

# -----------------------------------------------------
# RLOG
#
# compute RLOG transformation for QC (thus blind=FALSE)
# -----------------------------------------------------
rld<- rlog(dds, blind=FALSE)
comment(rld) = comment(dds)
print(paste("done rlog:", date()))

counts <- assay(rld)
colnames(counts) <- gsub(".Aligned.sortedByCoord.out.bam","",colnames(counts))

write.csv(counts, "./results/2021-04-27_rld_count_data.csv")

# -----------------------------------------------------
# save.image
#
# -----------------------------------------------------
save.image(file=paste("3b","dds_rlog",comment(dds)["by"],comment(dds)["soMode2"],"RData",sep="."))
print(paste("done save:", date()))

# -----------------------------------------------------
#quit(save="no")
# -----------------------------------------------------
