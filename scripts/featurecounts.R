# library
library(Rsubread)
library(Rgb)

# set arguements
args = commandArgs(trailingOnly=TRUE)

# get bam list file 
list_of_files <- strsplit(args[1], split = " ")[[1]]
file <- featureCounts(files=list_of_files, annot.ext=args[2], isGTFAnnotationFile=TRUE, 
                      GTF.featureType="exon", GTF.attrType="gene_id", isPairedEnd=TRUE)
saveRDS(file, file = "count.rds")
