args <- commandArgs()

help <- function(){
    cat("ChIP_geneProperties.R :
-Annotate all transcripts from a txdb with ChIPseq occupancy at the promoter,
-transcript total coverage and total coverage with each of the specified number of sections.
-You can optionally also provide a peak file that and the output table will contain a column (0/1) for overlap presence.
-The output table will also contian refseq_mrna id for the gene, the ensembl_gene_id biotype and the external gene name.
-This script is set up for ChIP-seq and does not take the stranded coverage into account.\n")
    cat("Usage: \n")
    cat("--txdbFile    : transcript data base (.txdb)                   [required]\n")
    cat("--assembly    : genome (hg19, mm9, mm10, dm3, sacCer3)         [hg19]\n")
    cat("--promUp      : number of nucleotides upstream of the tss      [required]\n")
    cat("--promDown    : number of nucleotides downstream the tss       [required]\n")
    cat("--numSections : number of sections to break the gene body into [default = 4 ]\n")
    cat("--bwFile      : bigWigFile                                     [required]\n")
    cat("--peakFile    : bed file                                       [optional]\n")
    cat("--numCores    : number of cores to use                         [default = 6 ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    txdbFile    <- sub( '--txdbFile=', '', args[grep('--txdbFile=', args)] )
    assembly    <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    promUp      <- as.numeric( sub( '--promUp=', '', args[grep('--promUp=', args)] ))
    promDown    <- as.numeric( sub( '--promDown=', '', args[grep('--promDown=', args)] ))
    numSections <- sub( '--numSections=', '', args[grep('--numSections=', args)])
    bwFile      <- sub( '--bwFile=', '', args[grep('--bwFile=', args)])
    peakFile    <- sub( '--peakFile=', '', args[grep('--peakFile=', args)])
    Cores       <- sub( '--numCores=', '',args[grep('--numCores=',args)])

}

if (identical(Cores,character(0))){
   Cores <- 6
}else{
    Cores <- as.numeric(Cores)
}

## set defaults
if (identical(numSections,character(0))){
   numSections <- 4
}else{
    numSections <- as.numeric(numSections)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
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
library(GenomicFeatures)
library(parallel)
library(biomaRt)

if (assembly == "hg19") {
    organism <- Hsapiens
    ## this is for enseble version 75
    bm <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    if(!(file.exists(annoFile))) {  
        anno <- getBM(mart=bm, attributes=c('ensembl_gene_id'
                                           ,'external_gene_name'
                                           ,'gene_biotype'
                                           ,'refseq_mrna'
                                           ,'refseq_ncrna'
                                           ,'entrezgene'
                                           ,'description'))
        anno[anno$refseq_mrna=="",  "refseq_mrna"] <- NA
        anno[anno$refseq_ncrna=="", "refseq_ncrna"] <- NA
        save(anno, file=annoFile)
    }else{
        anno <- get(load(annoFile))
    }}
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

##############
## tss average coverage
##############

## load model
txdb <- loadDb(txdbFile)

## filter out chroms
seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("_|\\d+.1$|^M$",seqlevels(txdb), invert=TRUE)]
seqlevels(txdb) <- sub("MT","M", seqlevels(txdb))
seqlevels(txdb) <- sub("Mito","M", seqlevels(txdb))

if(length(grep("^chr",seqlevels(txdb)))==0){
    seqlevels(txdb) <- sub("^","chr", seqlevels(txdb))
}

## make gnModel to combine with counts for only looking at the browser
## rpkms if calculated with need to use the exons below
gnModel          <- transcriptsBy(txdb, 'gene')
gnModel          <- unlist(gnModel)## Gets the genomic region convered by transcripts
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]

## add gene info to gnModel
gnModel$gene_id            <- names(gnModel)
iv                         <- match(names(gnModel), anno$ensembl_gene_id)
gnModel$external_gene_name <- anno[iv, "external_gene_name"]
gnModel$gene_biotype       <- anno[iv, "gene_biotype"]
gnModel$refseq_mrna        <- anno[iv, "refseq_mrna"]
gnModel$refseq_ncrna       <- anno[iv, "refseq_ncrna"]

## filter for lincRNA and protein coding RNA
gnModel <- gnModel[gnModel$gene_biotype %in% c("protein_coding", "lincRNA")]

## define tss region and remove names because they are redundant
TSS            <- promoters(gnModel, upstream=promUp, downstream=promDown)
TSS$gene_id    <- names(TSS)
names(TSS)     <- NULL
names(gnModel) <- NULL

sumCov <- function(bw,model){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    cat("get coverage for peak region\n")    
    mean.cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end){
            sum(bw.peak.cov[[seqname]][start:end])
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end)
    })
    mean.cov <- data.frame(mean.cov)
    rownames(mean.cov) <- model$tx_name
    mean.cov
}

## sum the coverage in the defined promoter region
txTssCov <- do.call(cbind,mclapply(bwFile,model=TSS,sumCov,mc.cores=1))
elementMetadata(gnModel)[[paste0("tssUp", promUp, "Down", promDown, "cov")]] <- txTssCov$mean.cov

## entire tx coverage
entireTx <- do.call(cbind,mclapply(bwFile,model=gnModel,sumCov,mc.cores=1))
gnModel$entireTxCov <- entireTx[,1]

########################
## section the gene
########################

for (i in 1:numSections){
    # Set up GB sections, positive strand
    print(paste("Section",i,"positive"))
    gbSection                                        <- gnModel
    gbSection$lengths                                <- width(gbSection)
    ## plus strand
    start(gbSection[which(strand(gbSection) =="+")]) <- start(gnModel[which(strand(gnModel) =="+")]) + (i-1)*(round(gbSection[which(strand(gnModel) =="+")]$lengths/numSections))
    end(gbSection[which(strand(gbSection)=="+")])    <- start(gnModel[which(strand(gnModel)=="+")]) + i*round(gbSection[which(strand(gnModel) =="+")]$lengths/numSections)
    ## negative strand; in granges for minus strand the start is the end
    print(paste("Section",i,"negative"))
    start(gbSection[which(strand(gnModel) =="-")])   <- end(gnModel[which(strand(gnModel) =="-")]) - i*(gbSection[which(strand(gnModel) =="-")]$lengths/numSections)
    end(  gbSection[which(strand(gnModel)=="-") ])   <- end(gnModel[which(strand(gnModel) =="-")]) - (i-1)*(gbSection[which(strand(gnModel) =="-")]$lengths/numSections)
    ## calc the total cov in that region
    sectionSum                                       <- do.call(cbind,mclapply(bwFile,model=gbSection,sumCov,mc.cores=1))
    elementMetadata(gnModel)[[paste0("section", i, "totalCov")]] <- sectionSum$mean.cov
}

if ( ! identical( peakFile, character(0)) ){
    peaks                                               <- import.bed(peakFile)
    ol                                                  <- as.data.frame(findOverlaps(gnModel, peaks))
    gnModel$peakOverLap                                 <- 0
    mcols(gnModel[unique(ol$queryHits)])["peakOverLap"] <- 1
}


df <- as.data.frame(gnModel)

write.table(df, file=sub("bw$", paste0(numSections, "geneSectionsSumCov.txt"), basename(bwFile)), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
