args <- commandArgs()

help <- function(){
    cat("metaGeneMatrix.R :
- From a provided GRnages obnject of genes of interest make metaGene matrix of defined windows of upstream of the TSS
  and down stream the Tes.
- This script is for one bigWig containing the both the plus and minus coverage.\n")
    cat("Usage: \n")
    cat("--regions     : GenomicRanges object with the transcripts of interest         [required]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                         [default = hg19 ]\n")    
    cat("--bwFile      : path to bigWig File                                           [required]\n")
    cat("--upStream    : how many bases upstream of the Tss                            [default = 2000 ]\n")
    cat("--downStream  : how many bases downstream of the Tes                          [default = 2000 ]\n")
    cat("--numCores    : number of cores to use                                        [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                 [default = basename(bigWigFile) ]\n")
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
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])

}

#regions <- "tables/Pol2_293T_DMSO_817_rep1.filteredProteinCodingTx.rda"
#assembly <- "hg19"
#bwFile <- "data_chipseq/PolII_468_293T_817_rep1.bw"
#upStream <- 2000
#downStream <- 2000
#Cores <- 10
#outName <-paste("tables/metaGene", sub(".bw", "", basename(bwFile)), sep="/")

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

TSS <- promoters(Model, upstream = upStream, downstream = 0)
TES <- resize(Model, fix = 'end', width = 1)
TES <- promoters(TES,upstream=0,downstream=downStream)
Body <- Model

head(ranges(TSS))
head(ranges(Body))
head(ranges(TES))

fname <- sub("$", ".metaGene.rda", outName)
 
Bin <- function(bw,model,Size){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    seqlengths(bw.peak.cov) <- seqlengths(organism)[seqlevels(bw.peak.cov)]
    cat("approx coverage in region to set size\n")
    cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end,strand){
            r <- bw.peak.cov[[seqname]][start:end]
            if(strand == '-'){r <- rev(r)}
            r <- approx(r,n=Size)$y
            return(r)
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end,as.character(strand))
    })
    cov <- as.data.frame(t(cov))
    rownames(cov) <- model$gene_id
    cov
}

tssMat  <- Bin(bw=bwFile, model=TSS,  Size=100)    
bodyMat <- Bin(bw=bwFile, model=Body, Size=400)    
tesMat  <- Bin(bw=bwFile, model=TES,  Size=100)    

stopifnot(rownames(tssMat)==rownames(tesMat))
stopifnot(rownames(tssMat)==rownames(bodyMat))

df <- data.frame(cbind(tssMat, bodyMat, tesMat))
rownames(df) <- rownames(tssMat)
save(df,file=fname)
