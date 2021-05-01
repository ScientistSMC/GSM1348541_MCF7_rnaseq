# 2.compute_summerizeOverlaps.R
#
#--------------------------------------------------
# libraries
#--------------------------------------------------
load("1z.RData") # load setup info
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")

#----------------------------------------------------------------------
# compute read counts
#
# http://www.bioconductor.org/help/workflows/rnaseqGene/#count
#---------------------------------------------------------------------- 
library("GenomicAlignments")
library("BiocParallel")


# use one core
register(SerialParam())

# several cores on one machine - by hand, SGE or SLURM
NCORES = strtoi(Sys.getenv("SLURM_CPUS_PER_TASK",unset=Sys.getenv("NSLOTS",unset="1")))
register(MulticoreParam(workers = 7), default = TRUE)
registered()[1]

#
# BAM(sort=pos)
# cheaha 15 CPU @ 2 G.RAM => ~35 minutes
#
print(paste("start:", date()))
soMode = "IntersectionStrict" # "Union", "IntersectionStrict", or "IntersectionNotEmpty"
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode=soMode,
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
# asociate sample info with "se" (SummarizedExperiment) object
(colData(se) <- DataFrame(sampleTable[match(colnames(se), rownames(sampleTable)),]))
comment(se) = c(comment(ebg), comment(sampleTable), soMode=soMode)

save.image(file=paste("2z","se",comment(se)["by"],comment(se)["soMode"],"RData",sep="."))
print(paste("done:", date()))

#save.image("2.RData")

#######################################################################
quit(save="yes")
#######################################################################
