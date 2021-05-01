# 1.setup_meta.R


#--------------------------------------------------
# load sample info
#--------------------------------------------------
#

sampleFile <- file.path("../sample_data","GSM1348541_MCF7_rnaseq_sample_meta_data.csv")
(sampleTable <- read.delim(sampleFile,sep=",",comment="#",header=T, fill=T, stringsAsFactors=F))

rownames(sampleTable)= sampleTable$Sample
comment(sampleTable) =c(sampleFile=sampleFile, loadedOn = date())

