# 1c.setup_se_bam_paths.R
load("1c.setupMetaAndGenes.RData")
# --------------------------------------------------
# bam files: tophat sort=pos
#
# In Bioconductor > 2.12 paired-end files do not need to
#          be sorted by 'qname'. Instead set 'asMates=TRUE' in the
#          'BamFile' when using the 'readGAlignmentsList' function from
#          the 'GenomicAlignments' package.
#
# BamFile(.., asMates=T) => PE
#--------------------------------------------------
# build bam file paths - STAR
filenames <- list.files(path = "../STAR_hg38_alignment/trimmed", pattern = ".Aligned.sortedByCoord.out.bam$", all.files = FALSE,
                         full.names=T, recursive=F,
                         ignore.case=F)

# validate if they exist
print(paste("bams, pos   sort: ", sum(FALSE==(file.exists(filenames))),    " missing"))

library("stringr")
library("dplyr")
# make "bam" file reader
library("Rsamtools")
bamfiles  <- BamFileList(filenames, yieldSize=2000000, asMates=FALSE)

# check chrom names
seqinfo(bamfiles[1])

# get FQ names from BAM file's Command Line
bamFastqs = data.frame(t(sapply(filenames,
               FUN = function (bamfile) {
                 bamHeader = scanBamHeader(bamfile);
                 bamCO=bamHeader[[1]]$text$"@CO";
                 alignArgs = unlist(strsplit(x=bamCO, split=" "));
                 fqs = grep(x=alignArgs,".fq|.fastq", value=T); 
                 return(c(bamFilename=bamfile, bamCO=bamCO,bamFastq1=fqs[1], bamFastq2=fqs[2]))
               }
)), stringsAsFactors=F)

# Creating a new column in the data frame where the sample id is extracted from the bam file name to be 
# able to match the bamfastq data to the sampleTable data and creating a column with a "simpliefied" bam file name
bamFastqs$sampleid <- gsub("../STAR_hg38_alignment/trimmed/|.Aligned.sortedByCoord.out.bam", "",bamFastqs$bamFilename)
bamFastqs$simplebamFilename <- gsub("../STAR_hg38_alignment/trimmed/", "",bamFastqs$bamFilename)

# Now joining the two data frames and making sure the data matches line by line correctly
sampleTable2 = full_join(sampleTable,bamFastqs, by=c("Sample" ="sampleid"))
comment(sampleTable2) = comment(sampleTable)
sampleTable=sampleTable2
rownames(sampleTable) <-sampleTable$simplebamFilename
sampleTable$simplebamFilename <- gsub(".Aligned.sortedByCoord.out.bam","",sampleTable$simplebamFilename)
rm("sampleTable2")


# save fastq filenames in sampleTable
save.image("1z.RData")
