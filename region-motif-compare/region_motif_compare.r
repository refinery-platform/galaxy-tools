# Name: region_motif_compare.r
# Description: Reads in two count files and determines enriched and depleted
# motifs (or any location based feature) based on poisson tests and gc
# corrections. All enrichment ratios relative to overall count / gc ratios.
# Author: Jeremy liu
# Email: jeremy.liu@yale.edu
# Date: 15/02/11
# Note: This script can be invoked with the following command
# R --slave --vanilla -f ./region_motif_compare.r --args <workingdir> <pwm_file> 
#        <intab1> <intab2> <enriched_tab> <depleted_tab> <plots_png>
# <workingdir> is the directory where plotting.r is saved
# Dependencies: region_motif_data_manager, plotting.r, 

# Auxiliary function to concatenate multiple strings
concat <- function(...) {
	input_list <- list(...)
	return(paste(input_list, sep="", collapse=""))
}

# Supress all warning messages to prevent Galaxy treating warnings as errors
options(warn=-1)

# Set common and data directories
args <- commandArgs()
workingDir = args[7]
pwmFile = unlist(strsplit(args[8], ','))[1]  # If duplicate entires, take first one

# Set input and reference files
inTab1 = args[9]
inTab2 = args[10]
enrichTab1 = args[11]
enrichTab2 = args[12]
plotsPng = args[13]

# Load dependencies
source(concat(workingDir, "/plotting.r"))

# Auxiliary function to read in tab file and prepare the data
read_tsv <- function(file) {
	data = read.table(file, sep="\t", stringsAsFactors=FALSE)
	names(data)[names(data) == "V1"] = "motif"
	names(data)[names(data) == "V2"] = "counts"
	return(data)
}

startTime = Sys.time()
cat("Running ... Started at:", format(startTime, "%a %b %d %X %Y"), "...\n")

# Loading motif position weight matrix (pwm) file 
cat("Loading motif postion weight matrices...\n")
lines = scan(pwmFile, what="character", sep="\n", quiet=TRUE)
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

# Loading input tab files
cat("Loading and reading input region motif count files...\n")
region1DF = read_tsv(inTab1)
region2DF = read_tsv(inTab2)
region1Counts = region1DF$counts
region2Counts = region2DF$counts
names(region1Counts) = region1DF$motif
names(region2Counts) = region2DF$motif

# Processing count vectors to account for missing 0 count motifs, then sorting
cat("Performing 0 count correction and sorting...\n")
allNames = union(names(region1Counts), names(region2Counts))
region1Diff = setdiff(allNames, names(region1Counts))
region2Diff = setdiff(allNames, names(region2Counts))
addCounts1 = rep(0, length(region1Diff))
addCounts2 = rep(0, length(region2Diff))
names(addCounts1) = region1Diff
names(addCounts2) = region2Diff
newCounts1 = append(region1Counts, addCounts1)
newCounts2 = append(region2Counts, addCounts2)
region1Counts = newCounts1[sort.int(names(newCounts1), index.return=TRUE)$ix]
region2Counts = newCounts2[sort.int(names(newCounts2), index.return=TRUE)$ix]

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

# Setting up output diagnostic plots, 4 in 1 png image
png(plotsPng, width=800, height=800)
xlab = "region1_count"
ylab = "region2_count"
lim = c(0.5, 5000)
layout(matrix(1:4, ncol=2))
par(mar=c(5, 5, 5, 1))

# Plot all motif counts along the linear correlation coefficient
plot.scatter(region1Counts+0.5, region2Counts+0.5, log="xy", xlab=xlab, ylab=ylab,
						 cex.lab=2.2, cex.axis=1.8, xlim=lim, ylim=lim*rValue)
abline(0, rValue, untf=T)
abline(0, rValue*2, untf=T, lty=2)
abline(0, rValue/2, untf=T, lty=2)
	
# Plot enriched and depleted motifs in red, housed in second plot    
plot.scatter(region1Counts+0.5, region2Counts+0.5, log="xy", xlab=xlab, ylab=ylab,
						 cex.lab=2.2, cex.axis=1.8, xlim=lim, ylim=lim*rValue)
points(region1Counts[indices]+0.5, region2Counts[indices]+0.5, col="red")
abline(0, rValue, untf=T)
abline(0, rValue*2, untf=T, lty=2)
abline(0, rValue/2, untf=T, lty=2)

# Apply and plot gc correction and loess curve
cat("Applying gc correction, rerunning poisson test...\n")
ind = which(region1Counts>5)
gc = gc[names(region2Counts)] # Reorder the indices of pwms to match input data
lo = plot.scatter(gc,log2(region2Counts/region1Counts),draw.loess=T, 
								xlab="gc content of motif",ylab=paste("log2(",ylab,"/",xlab,")"),
								cex.lab=2.2,cex.axis=1.8,ind=ind) # This function is in plotting.r
gcCorrection = 2^approx(lo$loess,xout=gc,rule=2)$y

# Recalculate p and q values, and filter for significant entries
pValueGC = sapply(seq(along=region1Counts),function(i) {
	poisson.test(c(region1Counts[i],region2Counts[i]),r=1/gcCorrection[i])$p.value
})
qValueGC=p.adjust(pValueGC,"fdr")
indicesGC = which(qValueGC<0.1 & abs(log2(region1Counts/region2Counts*gcCorrection))>log2(1.5))

# Plot gc corrected motif counts 
plot.scatter(region1Counts+0.5, (region2Counts+0.5)/gcCorrection, log="xy", 
						 xlab=xlab, ylab=paste(ylab,"(normalized)"), cex.lab=2.2, cex.axis=1.8,
						 xlim=lim, ylim=lim)
points(region1Counts[indicesGC]+0.5, 
			 (region2Counts[indicesGC]+0.5)/gcCorrection[indicesGC], col="red")
abline(0,1)
abline(0,1*2,untf=T,lty=2)
abline(0,1/2,untf=T,lty=2)

# Trim results, compile statistics and output to file
# Only does so if significant results are computed
if(length(indicesGC) > 0) {
	# Calculate expected counts and enrichment ratios
	cat("Calculating statistics...\n")
	nullExpect = region1Counts * gcCorrection
	enrichment = region2Counts / nullExpect

	# Reorder selected indices in ascending pvalue
	cat("Reordering by ascending pvalue...\n")
	indicesReorder = indicesGC[order(pValueGC[indicesGC])]

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
	write.table(outDF[indicesEnrich,], file=enrichTab1, quote=FALSE, 
							sep="\t", append=FALSE, row.names=FALSE, col.names=TRUE)
	write.table(outDF[indicesDeplete,], file=enrichTab2, quote=FALSE, 
							sep="\t", append=FALSE, row.names=FALSE, col.names=TRUE)
}

# Catch display messages and output timing information 
catchMessage = dev.off()
cat("Done. Job started at:", format(startTime, "%a %b %d %X %Y."),
		"Job ended at:", format(Sys.time(), "%a %b %d %X %Y."), "\n")
