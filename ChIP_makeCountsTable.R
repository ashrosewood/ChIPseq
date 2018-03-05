args <- commandArgs()

help <- function(){
    cat("ChIPseq_makeCountsTable.R :
- From a bed file or GRanges object make a counts table. This assumes you mapped with bowtie 1 and have no multimappers!
- Otherwise, you need to fitler your bam file before if the bam file has multimapped reads you do not want.
- Assumes that this is an unstranded library.
- Two files are returned the counts for the regions and library sizes for the bamfiles.\n")
    cat("Usage: \n")
    cat("--bamDir     : bamFiles directory to make counts table                                [required]\n")
    cat("--Pattern    : grep pattern of samples that you want included in table                [default = .*.bam]\n")
    cat("--regions    : Path to rda GRanges object, or bed file to make counts of              [required]
                       bed file must have chr, start, end, name\n")    
    cat("--assembly   : Genome (hg19, mm9, dm3, mm10 or sacCer3.                               [hg19]\n")
    cat("--numCores   : Number of cores used should not be higher than the number of bam files [required]\n")
    cat("--outName    : Path and prefix to output file (.counts.txt will be appended to end)   [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    assembly   <-sub('--assembly=', '', args[grep('--assembly=', args)])
    bamDir     <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
    Pattern    <-sub('--Pattern=', '', args[grep('--Pattern=', args)])
    regions    <-sub('--regions=', '',args[grep('--regions=',args)])
    numCores   <- sub('--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub('--outName=', '', args[grep('--outName=',args)])
}

if (identical(assembly,character(0))){
    assembly <- "hg19"
}

if (identical(Pattern,character(0))){
    Pattern <- ""
}

cat("bamDir:",bamDir , sep="\n")

library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)
library(GenomicRanges)

print(assembly)
if (assembly == "hg19") organismStr <- "Hsapiens"
if (assembly == "mm9") organismStr <- "Mmusculus"
if (assembly == "mm10") organismStr <- "Mmusculus"
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if (assembly == "hg19") {
    organism <- Hsapiens
}
if (assembly == "mm9") {
    organism <- Mmusculus
}
if (assembly == "mm10") {
    organism <- Mmusculus
}
if (assembly == "sacCer3") {
    organism <- Scerevisiae
}
if (assembly == "dm3") {
    organism <- Dmelanogaster
}

## load regions
if( length(grep(".bed$", regions)) > 0 ){
    print("regions are in bed format")
    Model <- import.bed(regions)
    names(Model) <- Model$name
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}else if( length(grep(".rda$", regions)) > 0 ){
    print("regions are in rda GRanges format")
    Model          <- get(load(regions))
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}

if(exists("Model")==FALSE){
    print("regions not provided in bed or rda format")
}

BFL<- list.files(bamDir,pattern=".bam", full.names=TRUE)
## if no pattern it will keep them all
BFL <- BFL[grep(Pattern,BFL,invert=FALSE)]
BFL

#setwd(bamDir)
#BFL <-dir(pattern='.*.bam$',recursive=FALSE)
#BFL <- 
   
## this param only reads in unique reads
                                        #param <- ScanBamParam(what='mapq',tag='NH',tagFilter=list(NH=1))
param <- ScanBamParam(what='mapq',tag='NH')

counter <-function(BF,mapq=10){
    aln                        <- readGAlignments(BF,param=param) # Read in gapped single end bam file
    seqlevels(aln)             <- sub("MT","M", seqlevels(aln))
    seqlevels(aln)             <- sub("Mito","M", seqlevels(aln))
    if( length(grep("^chr",seqlevels(aln)))==0 ){
        seqlevels(aln)         <- sub("^","chr", seqlevels(aln))
    }
    strand(aln)                <- "*"
    seqlevels(aln,force=TRUE)  <- seqlevels(aln)[grep("chrM",seqlevels(aln), invert=TRUE)] # throw out chrM
    ovlp                       <- countOverlaps(aln, Model) #Count how many genes each alignment has; remove reads that map to 2 genes
    aln                        <- aln[ovlp==1 ]             # Remove reads overlapping more than one gene
    counts                     <- countOverlaps(Model, aln) # Count the number of reads for each gene
    counts                                                  # Return counts for each bam
}

counts        <- do.call(cbind, mclapply(BFL, counter, mc.cores=numCores, mc.preschedule=FALSE))
counts        <- as.data.frame(counts)
names(counts) <- sub(".bam","", basename(BFL))

counts <- cbind(as.data.frame(Model),counts)

write.table(counts, file=sub("$", ".counts.txt", outName),sep="\t", quote=F, col.names=TRUE, row.names=FALSE)

## add counts for libsize
libsize <-function(BF){
    aln <- readGAlignments(BF,param=param)         # Read in gapped single end bam file
    seqlevels(aln) <- sub("MT","M", seqlevels(aln))
    seqlevels(aln) <- sub("Mito","M", seqlevels(aln))
    if(length(grep("^chr",seqlevels(aln)))==0){
        seqlevels(aln) <- sub("^","chr", seqlevels(aln))
    }
    ## throw out chrM
    seqlevels(aln,force=TRUE) <- seqlevels(aln)[grep("chrM",seqlevels(aln), invert=TRUE)]
    length(aln)
}

lib.counts <- do.call(cbind,mclapply(BFL,libsize,mc.cores=numCores,mc.preschedule=FALSE))
lib.counts <- as.data.frame(lib.counts)
names(lib.counts) <- sub(".bam","", basename(BFL))

write.table(lib.counts, file=sub("$", ".libraryCounts.txt", outName),sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
