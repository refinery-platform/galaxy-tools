# Name: region_motif_intersect.r
# Description: Takes a bed file of target regions and counts intersections
# of each motif (built in rdata database) and target regions.
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 14/07/02
# Note: This script is meant to be invoked with the following command
# R --slave --vanilla -f ./region_motif_intersect.r --args <workingdir> <inbed> <outtab>
# Dependencies: regions.r

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
  input_list <- list(...)
  return(paste(input_list, sep="", collapse=""))
}

# Set common and data directories
args <- commandArgs()
workingDir = args[7]
commonDir = "/Users/jeremyliu1/galaxy-dist/tools/my_tools"
motifDB = concat(workingDir, "/tools/my_tools/region_motif_db/pouya")
#motifDB = concat(workingDir, "/tools/my_tools/region_motif_db/JolmaJaspar")

# Set input and reference files, comment to toggle commmand line arguments
inBed = args[8]
outTab = args[9]

# Load dependencies
source(concat(commonDir, "/region_motif_lib/regions.r"))

# Auxiliary function to read in BED file
read_bed <- function(file) {
  return(read.table(file, sep="\t", stringsAsFactors=FALSE))
}

# Auxiliary function to write a data frame to BED file.
write_bed <- function(df, file) {
  write.table(df, file, quote=FALSE, sep="\t",
                  row.names=FALSE, col.names=FALSE)
}

startTime = Sys.time()
cat("Running ... Started at:", format(startTime, "%a %b %d %X %Y"), "...\n")

# Parsing motif database to generate a list of motif names
motifFilesFull = dir(motifDB, "RData", full.names=TRUE)
motifFilesPart = dir(motifDB, "RData", full.names=FALSE)
motifNames = gsub(".RData", "", motifFilesPart, ignore.case=TRUE)
names(motifFilesFull) = motifNames

# Loading input bed file, convert integer columns to numeric, name columns
regionsDF = read_bed(inBed)
dfTemp = sapply(regionsDF, is.integer)
regionsDF[dfTemp] = lapply(regionsDF[dfTemp], as.numeric)
names(regionsDF)[names(regionsDF) == "V1"] = "chr"
names(regionsDF)[names(regionsDF) == "V2"] = "start"
names(regionsDF)[names(regionsDF) == "V3"] = "end"

# Prepare list of chromosomes and call regions
chrs = unique(regionsDF$chr)
names(chrs) = chrs
callRegions = lapply(chrs, function(chr) {
  chrIndexes = which(regionsDF$chr == chr)
  cbind(regionsDF$start[chrIndexes], regionsDF$end[chrIndexes])  
})

# Counting intersections of regions in bed file and motifs in database
counts = sapply(motifFilesFull, function(motif) {
  load(motif) # pos data structure
  length(which(distance.to.closest.region.of.poslist(pos,callRegions)==0))
})

# Outputting intersection counts to tab delineated file
write.table(counts, outTab, quote=FALSE, sep="\t",
            row.names=TRUE, col.names=FALSE)
cat("Done. Job started at:", format(startTime, "%a %b %d %X %Y."),
    "Job ended at:", format(Sys.time(), "%a %b %d %X %Y."), "\n")
