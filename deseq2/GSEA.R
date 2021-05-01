library("tidyverse")
library("fgsea")
library("KEGGREST")
library("biomaRt")
library("pathview")
library("ggrepel")


#### Loading Data ####
res <- read.csv("./results/2021-04-29_res.estrogen_vs_control.csv", row.names=1, stringsAsFactors = FALSE)
res2 <- res[complete.cases(res$padj),]

#### GSEA with KEGG ####
data(examplePathways)
data(exampleRanks)

########### Get KEGG data ###########
# List of Pathways:
keggPathways <- keggList("pathway","hsa")
keggPathways <- data.frame(pathId=names(keggPathways),pathName=keggPathways)
keggPathways$pathId <- gsub("path:hsa","",keggPathways$pathId)
keggPathways$pathway <- rownames(keggPathways)

# Creating a dataframe of NON disease pathways only
nondisease_keggpathways <- keggPathways[grep("^05", keggPathways$pathId, invert=T),]

# Creating a dataframe of disease pathways only
disease_keggpathways <- keggPathways[grep("^05", keggPathways$pathId),]

# Pathway Genes:
keggPathwayGenes <- keggLink("pathway","hsa")
keggPathwayGenes <- data.frame(path=keggPathwayGenes,gene=names(keggPathwayGenes))
keggPathwayGenes$gene <- gsub("hsa:","",keggPathwayGenes$gene)
keggPathwayGenes$path <- gsub("pathway:","",keggPathwayGenes$path)
keggGeneSet <- list()
keggUniquePaths <- unique(keggPathwayGenes$path)

# remove disease pathways
keggUniquePaths <- keggUniquePaths[grep("path:hsa05", keggUniquePaths, invert=T)]
for(path in keggUniquePaths){
  keggGeneSet[[path]] <- unique(keggPathwayGenes$gene[keggPathwayGenes$path==path])
}

# Get table to convert ensbl ids to entrez ids
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensbl_entrez <- getBM(attributes = c("ensembl_gene_id","ensembl_gene_id_version", "entrezgene_id"), bmHeader = T, mart = ensembl)
ensbl_entrez <- ensbl_entrez[!(is.na(ensbl_entrez$'NCBI gene (formerly Entrezgene) ID')),]
  # sum(duplicated(ensbl_entrez$`Gene stable ID version`))
  # [1] 277
  # sum(duplicated(ensbl_entrez$`NCBI gene ID`))
  # [1] 2862
  # sum(duplicated(ensbl_entrez$`HGNC symbol`))
  # [1] 3268

# dim(female_vp_data)
# [1] 39203     7
# sum(res2$ENSG_ID %in% ensbl_entrez$`Gene stable ID version`)
# [1] 36320

res2$ENSG_ID <- rownames(res2)
res2$ENSG_ID_nv <- gsub("\\..*$","",res2$ENSG_ID)

# Rank gene orders (first add entez ids by matching on ensbl ids, 
  # then remove 58 entrez genes with dupicate entries, next order by 
  # p-val and fold change to get ranks)

ranked_data <- left_join(res2, ensbl_entrez, by = c("ENSG_ID_nv" = "Gene stable ID"))
ranked_data <- ranked_data[order(-abs(ranked_data$stat)),]
ranked_data <- ranked_data[complete.cases(ranked_data$`NCBI gene (formerly Entrezgene) ID`),]

ranked_data2 <- res2[order(-abs(res2$stat)),]
ranked_data3 <- data.frame("Gene"=ranked_data2$gene_symbol, "Rank"=ranked_data2$stat)
write.table(ranked_data3,file="2021_04_30_ranked_data_gsea.rnk",quote=F,sep="\t",row.names=F)
# write.csv(ranked_data3, "GSM_ranked_genes.rnk")

geneStats <- abs(ranked_data$stat)
names(geneStats) <- ranked_data$'NCBI gene (formerly Entrezgene) ID'

set.seed(18)

# Gene Set Enrichment Analysis:
KeggGeneGSEA <- fgsea(keggGeneSet,geneStats,nperm=10000,minSize=2,maxSize=Inf, scoreType = "pos")

test <- left_join(KeggGeneGSEA, keggPathways, by="pathway")
test$pathName <- gsub(" - Homo sapiens \\(human\\)", "", test$pathName)
colnames(test)[colnames(test)=="size"] <- "Size"

write.csv(test[,c(1:7,10)], "KeggGeneGSEA_results.csv", row.names = F)

### KEGG Pathway Analysis Plots ###
ggplot(filter(test, padj<=0.01 & pathName != "Metabolic pathways") , aes(x=reorder(pathName, -padj), y=-log10(padj), color=NES)) +
  labs(x="KEGG Pathway", y="-log10(q-value)", title="Significantly Enriched Pathways (adjusted p-value <0.01)") +
  geom_point(aes(size=Size)) + 
  scale_colour_gradient(low="green",high="red") +
  scale_size_area(limits = c(0,215)) +
  scale_y_continuous(limits = c(0,2.7), expand=c(0,0), breaks=c(seq(0,3,0.3))) +
  geom_segment( aes(x=pathName, xend=pathName, y=0, yend=-log10(padj)), lwd=1.2) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", vjust = 0.5), plot.title = element_text(size=10)) 


#### KEGG Enrichment Plots ####

plotEnrichment(keggGeneSet[["path:hsa04915"]],
               geneStats) + labs(title="Estrogen signaling pathway")



###### Pathview PLOTS #####
# Gene data:
path_Genes <- ranked_data[ranked_data$`NCBI gene (formerly Entrezgene) ID` %in% keggGeneSet[[paste0("path:hsa04390")]] & ranked_data$padj <= 0.05,]
geneData <- as.numeric(path_Genes$log2FoldChange)
names(geneData) <- as.character(path_Genes$`NCBI gene (formerly Entrezgene) ID`)

pathview(gene.data=geneData,
         pathway.id="04390",species="hsa",
         kegg.native=T,key.pos="topright",
         limit=list(gene=c(-1,1)),
         high=list(gene="green"),mid=list(gene="grey50"),
         low=list(gene="red"))







