args <- commandArgs()

help <- function(){
    cat("ChIP_CalcPausingIndex.R :
- From a provided GRnages object of transcripts calculate the average coverage in a defined promoter region and the remaining gene body region.
- It first filters for the best transcript with the highest coverage in the defined transcript tss region (from the tss to + txTssDown ).
- If a peakFile is provided, then take the transcripts with tss regions overlapping peaks.
- This script is set up for ChIP-seq and does not take the stranded coverage into account.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest    [required]\n")    
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                    [default = hg19]\n")    
    cat("--promUp      : number of nucleotides upstream of the tss                [required]\n")
    cat("--promDown    : number of nucleotides downstream the tss                 [required]\n")
    cat("--bodyLen     : if you want a smaller set body region (ex 2000)          [default = entire gene length w/o promoter]\n")
    cat("--bwFiles     : path to bigWig Files                                     [required]\n")
    cat("--bwPattern   : grep pattern for bigWigs to use quotes (ex. PolII.*293T) [optional]
                         if not provided uses all bw in path\n")
    cat("--numCores    : number of cores to use                                   [default = 10 ]\n")
    cat("--outName     : prefix to your out file names                            [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-h",args)) || !is.na(charmatch("-help",args)) ){
    help()
} else {
    regions   <- sub( '--regions=', '', args[grep('--regions=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    promUp    <- as.numeric( sub( '--promUp=', '', args[grep('--promUp=', args)] ))
    promDown  <- as.numeric( sub( '--promDown=', '', args[grep('--promDown=', args)] ))
    bodyLen  <- sub( '--bodyLen=', '', args[grep('--bodyLen=', args)] )
    bwFiles   <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])

}

## select all bigWig files 
bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
bws <- bws[grep(bwPattern,bws,invert=FALSE)]

if (identical(Cores,character(0))){
   Cores <- 6
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".filteredProteinCodingTx.rda|.bed", "", basename(regions))
}

## make output directory if does not exist
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

print(paste("regions:", regions))
print(paste("promUp:", promUp))
print(paste("promDown:", promDown))
print(paste("outName:", outName))
print(paste("bwFiles:", bwFiles))
print(paste("bwPattern:", bwPattern))

print(paste("bwFiles:", bws))

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
    ## this is for enseble version 75
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    anno <- get(load(annoFile))
    }
if (assembly == "mm9") {
    organism <- Mmusculus
}
if (assembly == "mm10") {
    organism <- Mmusculus
}
if (assembly == "sacCer3") {
    organism <- Scerevisiae
    annoFile <- "/projects/b1025/anno/biomaRt/sacSer3.Ens_78.biomaRt.geneAnno.Rdata"
    anno <- get(load(annoFile))
}
if (assembly == "dm3") {
    organism <- Dmelanogaster
}

###############################
## set up tss and body regions
###############################

## read in transcripts
if( length(grep(".bed$", regions)) > 0 ){
    boo <-  read.table(regions)
    names(boo) <- c("seqnames", "start", "end", "gene_id", "tx_name", "strand")
    gnModel <- as(boo, "GRanges")
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
}else{
    gnModel          <- get(load(regions))
    seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]
}

## define TSS
Tss <- promoters(gnModel ,upstream=promUp ,downstream=promDown)
ranges(Tss)

if (identical(bodyLen,character(0))){
    Body <- resize(gnModel, fix='end', width=width(gnModel)-promDown )
}else{
    bodyLen <- as.numeric(bodyLen)
    gnModel <- resize(gnModel,fix='start',width= bodyLen + promDown)
    Body <- resize(gnModel,fix='end',width=width(gnModel) - promDown )
}
ranges(Body)

###############################
## calc average coverage in tss and body regions
###############################

## set up function
Bin <- function(bw,model){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    cat("get coverage for peak region\n")
    mean.cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end){
            mean(bw.peak.cov[[seqname]][start:end])
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end)
    })
    mean.cov <- data.frame(mean.cov)
    rownames(mean.cov) <- model$gene_id
    colnames(mean.cov) <- sub(".bw", "", basename(bw))
    mean.cov
}

tssCov <- do.call(cbind,mclapply(bws,model=Tss,Bin,mc.cores=1))
bodyCov <- do.call(cbind,mclapply(bws,model=Body,Bin,mc.cores=1))

colnames(tssCov) <- sub("$", "_tss", colnames(tssCov))
colnames(bodyCov) <- sub("$", "_body", colnames(bodyCov))

stopifnot(rownames(tssCov)==rownames(bodyCov))

## combine tss and body dataframes
df <- cbind(as.data.frame(gnModel), tssCov, bodyCov)

## save the gnModel as a granges object and a tab delimited file
write.table(df
           ,file=paste0(outName, ".pausingIndexAverageCoverages.txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )
