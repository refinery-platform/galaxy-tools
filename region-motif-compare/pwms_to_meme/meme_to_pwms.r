# Name: 
# Description: 
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 15/03/24
# Note: This script can be invoked with the following command
# R --slave --vanilla -f ./.r --args <meme_file> <rdata_file_name>

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
    input_list <- list(...)
    return(paste(input_list, sep="", collapse=""))
}

args <- commandArgs()
infile = args[7]
outfile = args[8]

#load("jaspar.jolma.pwms.from.seq.RData")

#lines = scan("jaspar.jolma.pwms.from.seq.meme.txt", what="character", sep="\n")

lines = scan(infile, what="character", sep="\n")
indices = which(grepl("MOTIF", lines))
names(indices) = lapply(indices, function(i) {
    nameline = lines[i]
    name = substr(nameline, 7, nchar(nameline))
    })

pwms = sapply(indices, function(i) {
    infoline = unlist(strsplit(lines[i+1], " "))
    alength = as.numeric(infoline[4])
    width = as.numeric(infoline[6])
    subset = lines[(i+2):(i+2+width-1)]
    motiflines = strsplit(subset, " ")
    motif = t(do.call(rbind, motiflines))
    motif = apply(motif, 2, as.numeric)
    }, simplify=FALSE, USE.NAMES=TRUE)
# dimnames(motif) = list(c("a", "c", "g", "t"), NULL)
# NOTE: does not store (a, c, g, t) as row names of each pwm matrix
