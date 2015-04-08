# Name: pwms_to_meme.r
# Description: 
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 15/03/24
# Note: This script can be invoked with the following command
# R --slave --vanilla -f ./.r --args <pwms_r_data_file> <outfile>

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
    input_list <- list(...)
    return(paste(input_list, sep="", collapse=""))
}

args <- commandArgs()
infile = args[7]
outfile = args[8]

load(infile)

write("MEME version 4\n", file=outfile, append=FALSE)
write("ALPHABET= ACGT\n", file=outfile, append=TRUE)
write("strands: + -\n", file=outfile, append=TRUE)
write("Background letter frequencies\nA 0.250 C 0.250 G 0.250 T 0.250\n", file=outfile, append=TRUE)

for (motif in names(pwms)) {
    write(concat("MOTIF ", motif), file=outfile, append=TRUE)
    dims = dim(pwms[motif][[1]])
    write(concat("letter-probability matrix: alength= ", dims[1], " w= ", dims[2], 
        " nsites= N/A E= N/A"), file=outfile, append=TRUE)
    write.table(t(pwms[motif][[1]]), file=outfile, append=TRUE, quote=FALSE,
        row.names=FALSE, col.names=FALSE)
    write("", file=outfile, append=TRUE)
    #print(t(pwms[motif][[1]]))
}
