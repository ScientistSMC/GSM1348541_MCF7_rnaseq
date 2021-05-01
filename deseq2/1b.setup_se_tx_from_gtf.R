# 1b.setup_se_tx.R

#--------------------------------------------------
# read ensembl gene annotations
#--------------------------------------------------
library("GenomicFeatures")
library("mygene")

refOrg = "GRCh38.primary_assembly.genome"

# mm10
if( refOrg == "mm10" ) {
  refOrgAbbrev = "mm10"
  refOrgName = paste0("Mus musculus (",refOrg,")")
  refOrgTaxId = 10090
  # htseq-count was run with iGenomes UCSC/mm10/2015
  ref_gtf <- "/scratch/share/public_datasets/ngs/transcriptomes/igenomes/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"  
  
  # other gtf options
  #ref_ensebl84_gtf <- "ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.chr.gtf.gz"  # chr in chroms
} 
if( refOrg == "GRCh38.primary_assembly.genome" ) {
  # hg19
  refOrgAbbrev="hg38"
  refOrgName = paste0("Homo sapiens (",refOrg,")")
  refOrgTaxId = 9606
  ref_gtf <- "../ref/gencode_hg38/gencode.v31.annotation.gtf"
}
# cuffnorm was run with merged.gtf
#exp_gtf <- "../../cuffmerge/merged.gtf" # cuffnorm


#
# load the gtf we want
#
#(txdb_ref <- makeTxDbFromGFF(ref_gtf, format="gtf", circ_seqs=character())
txdbFile = paste0(basename(ref_gtf),".sqlite")
if( !file.exists(txdbFile) ) {
  (txdb <- makeTxDbFromGFF(ref_gtf, format="gtf", circ_seqs=character(), organism=refOrgName, taxonomyId=refOrgTaxId))
  comment(txdb) = c(ref_gtf=ref_gtf, txdbFile=txdbFile, refOrgName=refOrgName, taxonomyId=refOrgTaxId)
  saveDb(txdb, txdbFile)

} else {
  txdb = loadDb(txdbFile)     
  comment(txdb) = c(ref_gtf=ref_gtf, txdbFile=txdbFile, refOrgName=refOrgName, taxonomyId=refOrgTaxId)
}
# build the key2 symbol map
# go get the ENSG -> GeneSymbol map; the info is in the GTF, but we can't access it through the txdb!
# takes a few minutes
library("mygene")
baseGeneIDs = sub(keys(txdb), pattern = "[.][0-9]+$", replace="")
geneIdMap = getGenes(baseGeneIDs) 
#Querying chunk 1..61
geneIdMap[ifelse(is.na(geneIdMap$symbol=="XIST"), F, geneIdMap$symbol=="XIST"),]

# get unique symbols per gene ID
# http://stats.stackexchange.com/a/6682
l <-tapply(geneIdMap$symbol,geneIdMap$query,function(l)l)
#n <- max(sapply(l,length))
geneId2Symbol = t(sapply(l,function(x){res <- paste(x,collapse="|");c(x[1],res)}))
rm(l,baseGeneIDs)
colnames(geneId2Symbol)= c("symbol", "symbols")

geneKeyMap = data.frame(row.names=keys(txdb), baseId=sub(keys(txdb),pattern = "[.][0-9]+$", replacement = ""))
geneKeyMap = data.frame(geneKeyMap, geneId2Symbol[geneKeyMap$baseId,],stringsAsFactors = F)
#save(file="1.geneKeyMap.RData", geneKeyMap)

# build index to first instance of each baseID
cat("start", date())
xgi = function(x) { match(x=x,table=geneIdMap$query) }
tGeneI = sapply(geneKeyMap$baseId,FUN=xgi)
cat("done", date())
# use that to add "name" and "entrezgene" ids
geneKeyMap$name = geneIdMap[tGeneI,]$name
geneKeyMap$entrezgene = geneIdMap[tGeneI,]$entrezgene
rm(tGeneI,xgi)
#save(file="1.geneKeyMap.withName.RData", geneKeyMap)

# build gene->exon map
(ebg <- exonsBy(txdb, by="gene"))
comment(ebg)=c(by="gene", comment(txdb))
# build transcript->exon map
#(ebt <- exonsBy(txdb, by="tx"))
#comment(ebt)=c(by="tx", comment(txdb))
# build transcript->exon map with TX_NAME
(ebtn <- exonsBy(txdb, by="tx", use.names=T))
comment(ebtn)=c(by="tx", comment(txdb), use.names=T)
# 
# other sources of annotation
#library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#TxDb.Mmusculus.UCSC.mm10.knownGene

save.image("1c.setupMetaAndGenes.RData")

