args <- commandArgs()

help <- function(){
    cat("metaPeakMatrix.R :
- From a provided GRnages obnject of peaks of interest, approximate to the same length.
- Then make a matrix with upstream and downstream regions.
- This script is for one bigWig containing the both the plus and minus coverage.
- This script also assumes that there is no strand info. \n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest           [required]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                           [default = hg19 ]\n")    
    cat("--bwFile      : path to bigWig File                                             [required]\n")
    cat("--upStream    : how many bases upstream of the Tss                              [default = 2000 ]\n")
    cat("--downStream  : how many bases downstream of the Tes                            [default = 2000 ]\n")
    cat("--approxLen   : what size to approximate the regions to                         [default = 10000 ]\n")
    cat("--numCores    : number of cores to use                                          [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                   [default = basename(bigWigFile) ]\n")
    cat("--Bins        : Number of bins to take the average coverage after approximating [default = 50 ]
                            genes to the same length.                                                                 \n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFile     <- sub( '--bwFile=', '', args[grep('--bwFile=', args)])
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    approxLen  <- sub( '--approxLen=', '', args[grep('--approxLen=', args)] )
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    Bins       <- sub( '--Bins=', '',args[grep('--Bins=',args)])

}

#setwd("/projects/b1025/arw/analysis/kevin/leukemia_SEC")
#regions <- "tables/heatmaps/lmr/MLL_AFF1ct_DMSO_SEM_802_moi_Peaks0.countsAFF1ct_DMSO_SEM_802.lmr.bed"
#assembly <- "hg19"
#bwFile <- "data_chipseq/AFF1ct_DMSO_SEM_802.bw"
#outName <-paste("tables/heatmaps/lmr", sub(".bw", "", basename(bwFile)), sep="/")

if (identical(Bins,character(0))){
    Bins      <- 50
}else{
    Bins      <- as.numeric(Bins)
}
print(c("Bins:", Bins))

if (identical(approxLen,character(0))){
   approxLen <- 10000
}else{
   approxLen <- as.numeric(approxLen)
}
print(c("approxLen:", approxLen))

if (identical(upStream,character(0))){
   upStream <- 2000
}else{
   upStream <- as.numeric(upStream)
}

if (identical(downStream,character(0))){
   downStream <- 2000
}else{
   downStream <- as.numeric(downStream)
}

if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bw", "", basename(bwFile))
}

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens"
    species <- "Homo sapiens"}
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus"
    species <- "Mus musculus"}
if (assembly == "sacCer3") { organismStr <- "Scerevisiae"
    species <- "Saccharomyces cerevisiae"}
if (assembly == "dm3") { organismStr <- "Dmelanogaster"
    species <- "Drosophila melanogaster"}
if (assembly == "rn6") { organismStr <- "Rnorvegicus"
    species <- "Rattus norvegicus"}

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)
library(rtracklayer)
library(GenomicRanges)
library(parallel)

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
    Model          <- import.bed(regions)
    names(Model)   <- Model$name
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
    print("change strand from * to +")
    strand(Model)  <- '+'
}else if( length(grep(".rda$", regions)) > 0 ){
    print("regions are in rda GRanges format")
    Model          <- get(load(regions))
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
    print("change strand from * to +")
    strand(Model)  <- '+'
}

TSS <- promoters(Model, upstream = upStream, downstream = 0)
TES <- resize(Model, fix = 'end', width = 1)
TES <- promoters(TES, upstream = 0, downstream = downStream)
Body <- Model

head(ranges(TSS))
head(ranges(Body))
head(ranges(TES))

fname          <- sub("$", paste0("_metaPeak", "up", upStream, "down", downStream, ".rda"), outName)

Bin <- function(bw,model,Size){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    seqlengths(bw.peak.cov) <- seqlengths(organism)[seqlevels(bw.peak.cov)]
    cat("approx coverage in region to set size\n")
    Cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            r <- bw.peak.cov[[seqname]][start:end]
            if(strand == '-'){r <- rev(r)}
            r <- approx(r,n=Size)$y
            return(r)
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })
    cov           <- as.data.frame(t(Cov))
    if( Bins > 1){
        cat("Bin matrix\n")
        window.cov <- function(row){
            window <- as.integer(ncol(cov)/Bins)
            window.coverage <- lapply(0:(window-1), function(jump)
                rowMeans(row[as.numeric(jump*Bins+1):as.numeric((jump*Bins)+Bins)])
                )
            t(as.matrix(unlist(window.coverage)))
        }
        win <- mclapply(1:nrow(cov), function(i)
            window.cov(cov[i,]),mc.cores=Cores)    
        bin.mat      <- do.call(rbind, mclapply(win, as.numeric, mc.cores=Cores))        
        df           <- data.frame(bin.mat)
        rownames(df) <- model$name
    }else{
        df           <- cov
        rownames(df) <- model$name
    }
    df
}

tssMat  <- Bin(bw=bwFile, model=TSS,  Size=upStream)    
bodyMat <- Bin(bw=bwFile, model=Body, Size=approxLen)    
tesMat  <- Bin(bw=bwFile, model=TES,  Size=downStream)    

stopifnot(rownames(tssMat)==rownames(tesMat))
stopifnot(rownames(tssMat)==rownames(bodyMat))

df <- data.frame(cbind(tssMat, bodyMat, tesMat))
rownames(df) <- rownames(tssMat)
save(df,file=fname)
