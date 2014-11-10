# Name: region_motif_intersect.r
# Description: Takes a bed file of target regions and counts intersections
# of each motif (built in rdata database) and target regions.
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 14/07/02
# Note: This script is meant to be invoked with the following command
# R --slave --vanilla -f ./region_motif_intersect.r --args <workingdir> <db> <inbed> <outtab>
# <workingdir> is working directory of galaxy installation
# <db> types: "t" test, "p" pouya, "j" jaspar jolma, "m" mouse
# Dependencies: none

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
  input_list <- list(...)
  return(paste(input_list, sep="", collapse=""))
}

# Set common and data directories
args <- commandArgs()
workingDir = args[7]
dbDir = concat(workingDir, "/region_motif_db")
dbCode = args[8]
if (dbCode == "t") {
  motifDB = concat(dbDir, "/pouya_test_motifs.bed.bgz")
} else if (dbCode == "p") {
  motifDB = concat(dbDir, "/pouya_motifs.bed.bgz")
} else if (dbCode == "j") {
  motifDB = concat(dbDir, "/jaspar_jolma_motifs.bed.bgz")
} else if (dbCode == "m") {
  motifDB = concat(dbDir, "/mm9_motifs.bed.bgz")
} else {
  motifDB = concat(dbDir, "/pouya_motifs.bed.bgz")
}

# Set input and reference files, comment to toggle commmand line arguments
inBed = args[9]
outTab = args[10]

# Auxiliary function to read in BED file
read_bed <- function(file) {
  return(read.table(file, sep="\t", stringsAsFactors=FALSE))
}

startTime = Sys.time()
cat("Running ... Started at:", format(startTime, "%a %b %d %X %Y"), "...\n")

# Load dependencies
cat("Loading dependencies...\n")
suppressPackageStartupMessages(library(Rsamtools, quietly=TRUE))

# Initializing hash table (as env) with motif names and loading tabix file
cat("Loading motif database and initializing hash table...\n")
motifTable = new.env()
motifTbx <- TabixFile(motifDB)

# Loading input bed file, convert integer columns to numeric, name columns
cat("Loading region file...\n")
regionsDF = read_bed(inBed)
dfTemp = sapply(regionsDF, is.integer)
regionsDF[dfTemp] = lapply(regionsDF[dfTemp], as.numeric)
names(regionsDF)[names(regionsDF) == "V1"] = "chr"
names(regionsDF)[names(regionsDF) == "V2"] = "start"
names(regionsDF)[names(regionsDF) == "V3"] = "end"

# Filtering regions to exclude chromosomes not in motif database
cat("Determining intersection counts...\n")
motifTbxChrs = seqnamesTabix(motifTbx)
regionsDFFilter = subset(regionsDF, chr %in% motifTbxChrs)

# Loading regions into GRanges object and scanning motif tabix database
# Region end is incremented by 1 since scanTabix querying is inclusive for
# position start but exclusive for position end.
param = GRanges(regionsDFFilter$chr, IRanges(regionsDFFilter$start, 
                end=regionsDFFilter$end + 1))
regionsIntersects = scanTabix(motifTbx, param=param)

# Parsing result list and updating motif count hash table
cat("Parsing result list...\n")
for(regionIntersects in regionsIntersects) {
  for(regionIntersect in strsplit(regionIntersects, " ")) {
    intersectMotif = strsplit(regionIntersect, "\t")[[1]][4]
    if(is.null(motifTable[[intersectMotif]])) {
      motifTable[[intersectMotif]] = 1
    } else {
      motifTable[[intersectMotif]] = motifTable[[intersectMotif]] + 1
    }
  }
}

# Converting motif count hash table to an integer vector for output
counts = integer(length = length(ls(motifTable)))
names(counts) = ls(motifTable)
for(motifName in ls(motifTable)) {
  counts[motifName] = as.integer(motifTable[[motifName]])
}

# Outputting intersection counts to tab delineated file
cat("Outputting to file...\n")
write.table(counts, outTab, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
cat("Done. Job started at:", format(startTime, "%a %b %d %X %Y."),
    "Job ended at:", format(Sys.time(), "%a %b %d %X %Y."), "\n")
