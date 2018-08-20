args <- commandArgs()

help <- function(){
    cat("ChIPseq_makeCountsTable.R :
- From a bed file or GRanges object make a counts table. This assumes you mapped with bowtie 1 and have no multi-mappers!
- Otherwise, you need to filter your bam file before if the bam file has multi-mapped reads you do not want.
- Assumes that this is an unstranded library.
- Two files are returned the counts for the regions after resizing and library sizes for the bam files.\n")
    cat("Usage: \n")
    cat("--bamDir     : bamFiles directory to make counts table                                 [required]\n")
    cat("--Pattern    : grep pattern of samples that you want included in table                 [default = .*.bam]\n")
    cat("--regions    : Path to rda GRanges object, or bed file to make counts of               [required]
                       bed file must have chr, start, end, name\n")
    cat("--type       : Tss, Tes, or Peaks (peaks will use entire provided region if window = 0)[required]\n")
    cat("--window     : total size of region around types above (use for type = Peaks)          [default = 0]\n")
    cat("--upStream   : distance upstream of the region to take (use for Tss and Tes)           [default = 50]\n")
    cat("--downStream : distance downstream of the region to take (use for Tss and Tes)         [default = 50]\n")
    cat("--assembly   : Genome (hg19, mm9, dm3, mm10 or sacCer3.                                [hg19]\n")
    cat("--numCores   : Number of cores used should not be higher than the number of bam files  [required]\n")
    cat("--outName    : Path and prefix to output file (.counts.txt will be appended to end)    [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    bamDir     <-sub('--bamDir=', '', args[grep('--bamDir=', args)])
    Pattern    <-sub('--Pattern=', '', args[grep('--Pattern=', args)])
    regions    <-sub('--regions=', '',args[grep('--regions=',args)])
    Type       <- sub( '--type=', '', args[grep('--type=', args)] )
    Window     <- sub( '--window=', '', args[grep('--window=', args)] )
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    assembly   <-sub('--assembly=', '', args[grep('--assembly=', args)])
    numCores   <- sub('--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub('--outName=', '', args[grep('--outName=',args)])
}

if (identical(assembly,character(0))){
    assembly <- "hg19"
}

if (identical(Pattern,character(0))){
    Pattern <- ""
}

if (Type == "Peaks" ){
    print("Type is Peaks")
    if (identical(Window,character(0))){
        print("Window not provided use peak regions")
        Window <- 0
    }else{
        Window <- as.numeric(Window)
        print(paste("Window provided as", Window))
    }
}

if (Type == "Tss" | Type == "Tes" ){
    print(paste("Type is", Type))
    if (identical(upStream,character(0))){
        print("upStream not provided use 50bp")
        upStream <- 50
    }else{
        upStream <- as.numeric(upStream)
        print(paste("upStream provided as", upStream))
    }
    if (identical(downStream,character(0))){
        print("downStream not provided use 50bp")
        downStream <- 50
    }else{
        downStream <- as.numeric(downStream)
        print(paste("downStream provided as", downStream))
    }
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

if( Type=="Tss" ){
    print("model is Tss")
    if( upStream > 0 | downStream > 0){
        print(paste("resizing Tss regions to upstream", upStream, "and downstream", downStream))
        Model.win <- promoters(Model, upstream=upStream, downstream=downStream)
    }else{
        print(paste("resizing Tss regions to upstream 0 and downstream 1"))
        Model.win <- promoters(Model, upstream=0, downstream=1)
    }    
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream), outName)
    idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    if( upStream > 0 | downStream > 0){
        print(paste("resizing Tes regions to upstream", upStream, "and downstream", downStream))
        Model.win      <- promoters(tes, upstream=upStream, downstream=downStream)
    }else{
        print(paste("resizing Tes regions to upstream 0 and downstream 1"))
        Model.win <- tes
    }
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream), outName)
    idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Peaks" ){
    print("model is Peaks")
    if( Window > 0 ){
        print(paste("resizing peak regions to a window of", Window))
        Model.win      <- resize(Model, width=Window, fix="center")
    }else{
        print("using provided peak regions")
        Model.win <- Model
    }
    fname <- sub("$", paste0("_", Type, Window), outName)
    idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}


BFL<- list.files(bamDir,pattern=".bam", full.names=TRUE)
## if no pattern it will keep them all
BFL <- BFL[grep(Pattern,BFL,invert=FALSE)]
BFL

#param <- ScanBamParam(what='mapq',tag='NH',tagFilter=list(NH=1))
param <- ScanBamParam(what='mapq',tag='NH')

counter <-function(BF,mapq=10){
    aln                        <- readGAlignments(BF,param=param)                          # Read in gapped single end bam file
    seqlevels(aln)             <- sub("MT","M", seqlevels(aln))
    seqlevels(aln)             <- sub("Mito","M", seqlevels(aln))
    if( length(grep("^chr",seqlevels(aln)))==0 ){
        seqlevels(aln)         <- sub("^","chr", seqlevels(aln))
    }
    strand(aln)                <- "*"
    seqlevels(aln,force=TRUE)  <- seqlevels(aln)[grep("chrM",seqlevels(aln), invert=TRUE)] # throw out chrM
    ovlp                       <- countOverlaps(aln, Model.win)                            #Count how many genes each alignment has; remove reads that map to 2 genes
    aln                        <- aln[ovlp==1 ]                                            # Remove reads overlapping more than one gene
    counts                     <- countOverlaps(Model.win, aln)                            # Count the number of reads for each gene
    counts                                                                                 # Return counts for each bam
}

counts        <- do.call(cbind, mclapply(BFL, counter, mc.cores=numCores, mc.preschedule=FALSE))
counts        <- as.data.frame(counts)
names(counts) <- sub(".bam","", basename(BFL))

counts <- cbind(as.data.frame(Model.win),counts)

write.table(counts, file=sub("$", ".counts.txt", fname),sep="\t", quote=F, col.names=TRUE, row.names=FALSE)

## add counts for libsize
libsize <-function(BF){
    aln <- readGAlignments(BF,param=param)              # Read in gapped single end bam file
    seqlevels(aln) <- sub("MT","M", seqlevels(aln))
    seqlevels(aln) <- sub("Mito","M", seqlevels(aln))
    if(length(grep("^chr",seqlevels(aln)))==0){
        seqlevels(aln) <- sub("^","chr", seqlevels(aln))
    }
    seqlevels(aln,force=TRUE) <- seqlevels(aln)[grep("chrM",seqlevels(aln), invert=TRUE)] # throw out chrM
    length(aln)
}

lib.counts <- do.call(cbind,mclapply(BFL,libsize,mc.cores=numCores,mc.preschedule=FALSE))
lib.counts <- as.data.frame(lib.counts)
names(lib.counts) <- sub(".bam","", basename(BFL))

write.table(lib.counts, file=sub("$", ".libraryCounts.txt", fname),sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
