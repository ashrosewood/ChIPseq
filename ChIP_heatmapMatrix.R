args <- commandArgs()

help <- function(){
    cat("heatMatrix.R :
- From a provided GRnages obnject make binned heatmap matrix of defined windows around peaks, tss, or tes.
- For smaller windows 5-10kb, I recommend binning in 25bp. For larger windows, I recommend 50bp bins.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest         [required]\n")
    cat("--type        : Peaks, Tss, or Tes                                            [required]\n")
    cat("--window      : total size of region around types above                       [default = 5kb]\n")
    cat("--bins        : number to bin coverage the window should be divisable by this [default = 25bp]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                         [default = hg19]\n")    
    cat("--bwFile      : path to bigWig File                                           [required]\n")
    cat("--numCores    : number of cores to use                                        [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                 [default = basename(bigWigFile) ]\n")
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
    Bins      <- sub( '--bins=', '', args[grep('--bins=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFile    <- sub( '--bwFile=', '', args[grep('--bwFile=', args)])
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])

}

#regions <- "tables/Pol2_293T_DMSO_817_rep1.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#Window <- 4000
#Bins <- 25
#Type <- "Tss"
#bwFile <- "data_chipseq/PolII_468_293T_817_rep1.bw"
#Cores <- 10
#outName <-paste("tables/heatmaps", sub(".bw", "", basename(bwFile)), sep="/")
    
if (identical(Cores,character(0))){
   Cores <- 6
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bw", "", basename(bwFile))
}

if (identical(Window,character(0))){
    Window <- 5000
}else{
    Window <- as.numeric(Window)
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

Model          <- get(load(regions))
seqinfo(Model) <- seqinfo(organism)[seqlevels(Model)]

if( Type=="Tss" ){
    print("model is Tss")
    Model.win      <- promoters(Model, upstream=Window/2, downstream=Window/2)
    Model.win$name <- Model.win$gene_id
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    Model.win      <- promoters(tes, upstream=Window/2, downstream=Window/2)
    Model.win$name <- Model.win$gene_id
}else if( Type=="Peaks" ){
    print("model is Peaks")
    Model.win      <- resize(Model, width=Window, fix="center")
    Model.win$name <- Model.win$name
}

matBin <- function(bw, model){
    fname <- sub("$", paste0("_", Type, Window, ".rda"), outName)
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
