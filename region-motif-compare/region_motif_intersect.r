# Name: region_motif_intersect.r
# Description: Takes a bed file of target regions and counts intersections
# of each motif (in separately installed tabix database) and target regions.
# Author: Jeremy Liu
# Email: jeremy.liu@yale.edu
# Date: 15/02/11
# Note: This script is meant to be invoked with the following command
# R --slave --vanilla -f ./region_motif_intersect.r --args <db_bgz> <db_tbi> <inbed> <outtab>
# Dependencies: region_motif_data_manager

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
  input_list <- list(...)
  return(paste(input_list, sep="", collapse=""))
}

# Retrive motif database path
args <- commandArgs()
motifDB_bgz = unlist(strsplit(args[7], ','))[1] # Handles duplicate entries in data table
motifDB_tbi = unlist(strsplit(args[8], ','))[1] # Just takes the first one

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
motifTbx <- TabixFile(motifDB_bgz)

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
