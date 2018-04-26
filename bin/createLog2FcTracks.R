args <- commandArgs()

help <- function(){
    cat("createFcTracks.R :
- Create fold cahnge tracks for two bigWig files.\n")
    cat("Usage: \n")
    cat("--conFile  : Control bigWig file path (FC denominator)                [required]\n")
    cat("--expFile  : Experiment bigWig file path (FC numerator)               [required]\n")    
    cat("--outName  : Prefix for your output file (no strand or .bw extension) [default = expFile]\n")    
    cat("--assembly : Genome build hg19, mm9, mm10, dm3 ect.                   [default = hg19]\n")
    cat("--pseudo   : pseudo count to prevent x/0                              [default = 0.001]\n")
    cat("\n")
    q()
}

## Save values of each argument

if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    conFile     <- sub('--conFile=', '', args[grep('--conFile=', args)])
    expFile     <- sub('--expFile=', '', args[grep('--expFile=', args)])
    outName     <- sub('--outName=', '', args[grep('--outName=', args)])
    assembly    <- sub('--assembly=', '', args[grep('--assembly=', args)])
    pseudo     <- sub('--pseudo=', '', args[grep('--pseudo=', args)])
}

print(paste("conFile:", conFile))
print(paste("expFile:", expFile))

if (identical(outName,character(0))){
   outName <- sub(".bw$", "log2FC.bw", expFile)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(pseudo,character(0))){
   pseudo <- 0.001
}else{
   pseudo <- as.numeric(pseudo)
}

library(rtracklayer)
library(GenomicRanges)

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster

con.bw  <- import.bw(conFile)
exp.bw  <- import.bw(expFile)

con.cov <- coverage(con.bw, weight='score')
exp.cov <- coverage(exp.bw, weight='score')    

log2FC <- log2( (exp.cov + pseudo) / (con.cov + pseudo) )
seqinfo(log2FC) <- seqinfo(organism)[seqlevels(log2FC)]
export.bw(log2FC, outName)

print("done")
