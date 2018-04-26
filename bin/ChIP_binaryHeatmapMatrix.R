args <- commandArgs()

help <- function(){
    cat("ChIP_binaryHeatmapMatrix.R :
- From a provided GRnages object make a binary heatmap matrix of defined windows around peaks, tss, or tes.
- Unlike ChIP_heatmapMatrix.R that uses coverage, this script gives a 1 or a zero for the presence of a peak.
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
    cat("--bedFile     : path to peaks bed File                                         [required]\n")
    cat("--numCores    : number of cores to use                                         [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                  [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions    <- sub( '--regions=', '', args[grep('--regions=', args)] )
    Type       <- sub( '--type=', '', args[grep('--type=', args)] )
    Window     <- sub( '--window=', '', args[grep('--window=', args)] )
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    Bins       <- sub( '--bins=', '', args[grep('--bins=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bedFile    <- sub( '--bedFile=', '', args[grep('--bedFile=', args)])
    Cores      <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])

}

#setwd("/projects/b1025/arw/analysis/kevin/SEC/")
#regions  <- "/projects/b1025/arw/analysis/kevin/SEC/tables/MYC_H2171_1021_nonTss.rda"
#bedFile  <- "/projects/b1025/arw/analysis/kevin/SEC/data_chipseq/AFF1_SW1271_1021.macsPeaks.bed"
#assembly <- "hg19"
#Window <- 5000
#Bins     <- 25
#Type     <- "Peaks"
#Cores <- 10
#outName <-paste("tables/heatmaps/MYC_H2171", sub(".macsPeaks.bed", "", basename(bedFile)), sep="/")
    
if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bed", "", basename(bedFile))
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
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".binary.rda"), outName)
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    Model.win      <- promoters(tes, upstream=upStream, downstream=downStream)
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".binary.rda"), outName)
}else if( Type=="Peaks" ){
    print("model is Peaks")
    Model.win      <- resize(Model, width=Window, fix="center")
    fname <- sub("$", paste0("_", Type, Window, ".binary.rda"), outName)
}

matBin <- function(BED,model,name){
    cat("importing:", BED, sep="\n")
    bed                       <- import.bed(BED)
    seqinfo(bed)              <- seqinfo(Hsapiens)[seqlevels(bed)]
    # must have the same chromosomes as the model!
    #seqlevels(bed,force=TRUE)   <- seqlevels(bed)[grep("chrY|chrM",seqlevels(bed), invert=TRUE)]
    #seqlevels(model,force=TRUE) <- seqlevels(model)[grep("chrY|chrM",seqlevels(model), invert=TRUE)]
    cat("calc coverage\n")
    bd.cov          <- coverage(bed)
    seqinfo(bd.cov) <- seqinfo(Hsapiens)[seqlevels(bd.cov)]
    cat("get coverage for peak region\n")
    cov             <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end){
            r <- bd.cov[[seqname]][start:end]
            r[seq(1,length(r),by=Bins)]
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end
        )})
    cat("convert list to matrix\n")
    mat             <- do.call(rbind, mclapply(cov, as.numeric,mc.cores=Cores))
    df              <- data.frame(mat)
    rownames(df)    <- model$name
    save(df,file=fname)
    print("done")
}

matBin(BED=bedFile, model=Model.win)
