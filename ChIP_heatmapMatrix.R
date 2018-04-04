args <- commandArgs()

help <- function(){
    cat("ChIP_heatmapMatrix.R :
- From a provided GRnages obnject make binned heatmap matrix of defined windows around peaks, tss, or tes.
- For smaller windows 5-10kb, I recommend binning in 25bp. For larger windows, I recommend 50bp bins. 
- Window is an option if type is Peaks, and upStream and downStream are exclusive to type either Tss or Tes.\n")
 
    cat("Usage: \n")
    cat("--regions     : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--type        : Peaks, Tss, or Tes                                             [required]\n")
    cat("--window      : total size of region around types above (use for type = Peaks) [default = 5kb]\n")
    cat("--upStream    : distance upstream of the region to take (use for Tss and Tes)  [default = 2.5kb]\n")
    cat("--downStream  : distance downstream of the region to take (use for Tss and Tes)[default = 2.5kb]\n")
    cat("--bins        : number to bin coverage the window should be divisable by this  [default = 25bp]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")    
    cat("--bwFile      : path to bigWig File                                            [required]\n")
    cat("--numCores    : number of cores to use                                         [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                  [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions   <- sub( '--regions=', '', args[grep('--regions=', args)] )
    Type      <- sub( '--type=', '', args[grep('--type=', args)] )
    Window    <- sub( '--window=', '', args[grep('--window=', args)] )
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    Bins      <- sub( '--bins=', '', args[grep('--bins=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFile    <- sub( '--bwFile=', '', args[grep('--bwFile=', args)])
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])

}

#setwd("/projects/b1025/arw/analysis/kevin/SEC/")
#regions  <- "tables/Pol2_293T_DMSO_817_rep1.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#Window   <- 4000
#Bins     <- 25
#Type     <- "Tes"
#bwFile <- "data_chipseq/PolII_468_293T_817_rep1.bw"
#Cores <- 10
#outName <-paste("tables/heatmaps", sub(".bw", "", basename(bwFile)), sep="/")
    
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

if (Type == "Peaks" ){
    print("Type is Peaks")
    if (identical(Window,character(0))){
        print("Window not provided use 5kb")
        Window <- 5000
    }else{
        Window <- as.numeric(Window)
        print(paste("Window provided as", Window))
    }
}

if (Type == "Tss" | Type == "Tes" ){
    print(paste("Type is", Type))
    if (identical(upStream,character(0))){
        print("upStream not provided use 2.5kb")
        upStream <- 2500
    }else{
        upStream <- as.numeric(upStream)
        print(paste("upStream provided as", upStream))
    }
    if (identical(downStream,character(0))){
        print("downStream not provided use 2.5kb")
        downStream <- 2500
    }else{
        downStream <- as.numeric(downStream)
        print(paste("downStream provided as", downStream))
    }
}

if (identical(Bins,character(0))){
    Bins <- 25
}else{
    Bins <- as.numeric(Bins)
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

## read in transcripts
if( length(grep(".bed$", regions)) > 0 ){
    Model <- import.bed(regions)
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}else{
    Model          <- get(load(regions))
    seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]
}

if( Type=="Tss" ){
    print("model is Tss")
    Model.win      <- promoters(Model, upstream=upStream, downstream=downStream)
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".rda"), outName)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    Model.win      <- promoters(tes, upstream=upStream, downstream=downStream)
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".rda"), outName)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Peaks" ){
    print("model is Peaks")
    Model.win      <- resize(Model, width=Window, fix="center")
    fname <- sub("$", paste0("_", Type, Window, ".rda"), outName)
    idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}

matBin <- function(bw, model){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    seqlengths(bw.peak.cov) <- seqlengths(organism)[seqlevels(bw.peak.cov)]
    cat("get coverage for peak region\n")
    cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            r <- bw.peak.cov[[seqname]][start:end]
            if(strand == '-'){r <- rev(r)}
            return(r)
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })
    cat("convert list to matrix\n")
    mat <- do.call(rbind, mclapply(cov, as.numeric,mc.cores=Cores))
    x <- mat
    cov <- data.frame(x)
    cat("bin the matrix\n")
    window.cov <- function(row){
        window <- as.integer(ncol(cov)/Bins)
        window.coverage <- lapply(0:(window-1), function(jump)
            rowMeans(row[as.numeric(jump*Bins+1):as.numeric((jump*Bins)+Bins)])
            )
        t(as.matrix(unlist(window.coverage)))
    }
    win <- mclapply(1:nrow(cov), function(i)
        window.cov(cov[i,]),mc.cores=Cores)    
    bin.mat <- do.call(rbind, mclapply(win, as.numeric, mc.cores=Cores))
    df <- data.frame(bin.mat)
    rownames(df) <- model$name
    save(df,file=fname)
    print("done")
}

matBin(bw=bwFile, model=Model.win)
