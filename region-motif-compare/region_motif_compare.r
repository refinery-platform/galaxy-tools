# Name: region_motif_compare.r
# Description: Reads in two count files and determines enriched and depleted
# motifs (or any location based feature) based on poisson tests and gc
# corrections. All enrichment ratios relative to overall count / gc ratios.
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 14/07/03
# Note: This script is meant to be invoked with the following command
# R --slave --vanilla -f ./enhancer_motif_compare.r --args <intab1> <intab2> <enriched_tab> <depleted_tab>
# Dependencies: plotting.r

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
  input_list <- list(...)
  return(paste(input_list, sep="", collapse=""))
}

# Set common and data directories
args <- commandArgs()
workingDir = args[7]
commonDir = concat(workingDir, "/tools/my_tools")
pwmFile = concat(workingDir, "/tools/my_tools/region_motif_db/pouya.pwms.from.seq.RData")
pwmFile2 = concat(workingDir, "/tools/my_tools/region_motif_db/jasparjolma..pwms.from.seq")

# Set input and reference files
inTab1 = args[8]
inTab2 = args[9]
enrichTab = args[10]
depleteTab = args[11]

# Load dependencies
source(concat(commonDir, "/region_motif_lib/plotting.r"))

# Auxiliary function to read in tab file and prepare the data
read_tsv <- function(file) {
  data = read.table(file, sep="\t", stringsAsFactors=FALSE)
  names(data)[names(data) == "V1"] = "motif"
  names(data)[names(data) == "V2"] = "counts"
  return(data)
}

# Auxiliary function for comparing two motif pwms
# Only used if including pwms in the outfile files
pwm.cors <- function(pwm1,pwm2) {
  nc = 21
  if(ncol(pwm1)<nc | ncol(pwm2)<nc) stop("nc<21")
  m1=max(sapply(-5:5,function(i) {
    ind1 = (1:nc) + i
    ind = which(ind1>0 & ind1<=nc)
    ind1 = ind1[ind]
    ind2 = ind1 -i
    cor(as.numeric(pwm1[,ind1]),as.numeric(pwm2[,ind2])) 
  }))
  pwm2 = pwm2[4:1,21:1]
  m2=max(sapply(-10:2,function(i) {
    ind1 = (1:nc) + i
    ind = which(ind1>0 & ind1<=nc)
    ind1 = ind1[ind]
    ind2 = ind1 -i
    cor(as.numeric(pwm1[,ind1]),as.numeric(pwm2[,ind2])) 
  }))
  max(m1,m2)
}

startTime = Sys.time()
cat("Running ... Started at:", format(startTime, "%a %b %d %X %Y"), "...\n")

# Loading motif position weight matrix (pwm) file and input tab file
cat("Loading and reading input region motif count files...\n")
load(pwmFile) # pwms data structure
##########
#temp = pwms
#load(pwmFile2)
#pwms = append(temp, pwms)
##########
region1DF = read_tsv(inTab1)
region2DF = read_tsv(inTab2)
region1Counts = region1DF$counts
region2Counts = region2DF$counts
names(region1Counts) = region1DF$motif
names(region2Counts) = region2DF$motif
# DATASET NAMES NOT INCLUDED IN INPUT, WHERE SHOULD THIS BE ENCODED?

# Generate gc content matrix
gc = sapply(pwms, function(i) mean(i[2:3,3:18]))

# Apply poisson test, calculate p and q values, and filter significant results
cat("Applying poisson test...\n")
rValue = sum(region2Counts) / sum(region1Counts)
pValue = sapply(seq(along=region1Counts), function(i) {
  poisson.test(c(region1Counts[i], region2Counts[i]), r=1/rValue)$p.value
})
qValue = p.adjust(pValue, "fdr")
indices = which(qValue<0.1 & abs(log2(region1Counts/region2Counts/rValue))>log2(1.5))

# Apply gc correction, recalculate p and q values, filter significant results
cat("Applying gc correction, rerunning poisson test...\n")
ind = which(region1Counts>5)
xlab = "region1_count"
ylab = "region2_count"
# This function is in plotting.r
# THIS FUNCTION CREATES Rplots.pdf IN THE CURRENT DIRECTORY, REMOVE?
gc = gc[names(region2Counts)]
lo = plot.scatter(gc,log2(region2Counts/region1Counts),draw.loess=T, 
                xlab="gc content of motif",ylab=paste("log2(",ylab,"/",xlab,")"),
                cex.lab=2.2,cex.axis=1.8,ind=ind)
gcCorrection = 2^approx(lo$loess,xout=gc,rule=2)$y
pValueGC = sapply(seq(along=region1Counts),function(i) {
  poisson.test(c(region1Counts[i],region2Counts[i]),r=1/gcCorrection[i])$p.value
})
qValueGC=p.adjust(pValueGC,"fdr")
indicesGC = which(qValueGC<0.1 & abs(log2(region1Counts/region2Counts*gcCorrection))>log2(1.5))

# Trim results, compile statistics and output to file
# Only does so if significant results are computed
# REMOVE THIS CONDITIONAL?
if(length(indicesGC) > 0) {
  # Calculate expected counts and enrichment ratios
  cat("Calculating statistics...\n")
  nullExpect = region1Counts * gcCorrection
  enrichment = region2Counts / nullExpect

  # Reorder selected indices in ascending pvalue
  cat("Reordering by ascending pvalue...\n")
  indicesReorder = indicesGC[order(pValueGC[indicesGC])]

  # Loop through and find repeated or similar motifs, determined by pwm gc content
  # This is used for inserting motif diagrams into output
  # cat("Finding repeated or similar motif pwms...\n")
  # pwm = pwms[indicesReorder]
  # similarExists = rep(F,length(indicesReorder))
  # if(length(indicesReorder)>1) {
  #   for(i in 2:length(indicesReorder)) {
  #     for(j in 1:(i-1)) {
  #       if(pwm.cors(pwm[[i]],pwm[[j]])>0.8) {
  #         similarExists[i]=T
  #         break
  #       }
  #     }
  #   }
  # }

  # Combine data into one data frame and output to two files
  cat("Splitting and outputting data...\n")
  outDF = data.frame(motif=names(pValueGC), p=as.numeric(pValueGC), q=qValueGC, 
                     stringsAsFactors=F, region_1_count=region1Counts, 
                     null_expectation=round(nullExpect,2), region_2_count=region2Counts,
                     enrichment=enrichment)[indicesReorder,]
  names(outDF)[which(names(outDF)=="region_1_count")]=xlab
  names(outDF)[which(names(outDF)=="region_2_count")]=ylab
  indicesEnrich = which(outDF$enrichment>1)
  indicesDeplete = which(outDF$enrichment<1)
  outDF$enrichment = ifelse(outDF$enrichment>1,
                            round(outDF$enrichment,3),
                            paste("1/",round(1/outDF$enrichment,3)))
  write.table(outDF[indicesEnrich,], file=enrichTab, quote=FALSE, 
              sep="\t", append=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(outDF[indicesDeplete,], file=depleteTab, quote=FALSE, 
              sep="\t", append=FALSE, row.names=FALSE, col.names=TRUE)
}
cat("Done. Job started at:", format(startTime, "%a %b %d %X %Y."),
    "Job ended at:", format(Sys.time(), "%a %b %d %X %Y."), "\n")
