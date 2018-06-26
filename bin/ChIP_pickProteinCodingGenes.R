args <- commandArgs()

help <- function(){
    cat("pickBestTSSfromTxdbandBW.withStrandMultiSections.R :
- This script generates a filtered gene list from a Txdb based on ChIP-seq in the tss region and filters for ensembl protein coding
  genes that are also refseq validated.
- It first filters for the best transcript with the highest coverage in the defined transcript tss region (from the tss to + txTssDown )
     I recommend not taking anything upstream as this could be from an alternate transcript.
- If a peakFile is provided, then take the transcripts only overlapping peaks.
- Filter for a minimal gene length.
- Filter for a minimal distance to the nearest gene after selecting the best transcript.
- Require a max rpm in the tss region of at least 1.
- This script is set up for ChIP-seq and does not take the stranded coverage into account\n")
    cat("Usage: \n")
    cat("--txdbFile    : transcript data base (.txdb)                             [required]\n")
    cat("--assembly    : genome (hg19, mm9, mm10, dm3, sacCer3)                   [default = hg19]\n")
    cat("--txTssDown   : number of nucleotides downstream tss for tx picking only [default = 50]\n")    
    cat("--minLength   : min transcript length                                    [default = 2000]\n")
    cat("--minDist     : min distance between genes                               [default = 2000]\n")
    cat("--bwFile      : bigWigFile                                               [required]\n")
    cat("--peakFile    : bed file                                                 [optional]\n")
    cat("--numCores    : number of cores to use                                   [default = 10 ]\n")
    cat("--outName     : prefix to your out file names                            [default = basename(bigWigFile) ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    txdbFile    <- sub( '--txdbFile=', '', args[grep('--txdbFile=', args)] )
    assembly    <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    txTssDown   <- sub( '--txTssDown=', '', args[grep('--txTssDown=', args)])
    minLength   <- sub( '--minLength=', '', args[grep('--minLength=', args)])
    minDist     <- sub( '--minDist=', '', args[grep('--minDist=', args)])
    numSections <- sub( '--numSections=', '', args[grep('--numSections=', args)])
    bwFile      <- sub( '--bwFile=', '', args[grep('--bwFile=', args)])
    peakFile    <- sub( '--peakFile=', '', args[grep('--peakFile=', args)])
    Cores       <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName     <- sub( '--outName=', '',args[grep('--outName=',args)])
}

if (identical(txTssDown,character(0))){
   txTssDown <- 51
}else{
    txTssDown <- as.numeric(txTssDown)
}

if (identical(minLength,character(0))){
   minLength <- 2000
}else{
    minLength <- as.numeric(minLength)
}

if (identical(minDist,character(0))){
   minDist <- 2000
}else{
    minDist <- as.numeric(minDist)
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

## make output directory if does not exist
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

print(paste("bwFile:", bwFile))
print(paste("peakFile:", peakFile))
print(paste("txdbFile:", txdbFile))
print(paste("outName:", outName))

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
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    if(!(file.exists(annoFile))) {  
        bm <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
        anno <- getBM(mart=bm, attributes=c('ensembl_gene_id'
                                           ,'ensembl_transcript_id'
                                           ,'external_gene_name'
                                           ,'gene_biotype'
                                           ,'transcript_biotype'
                                           ,'refseq_mrna'
                                           ,'refseq_ncrna'
                                           ,'entrezgene'
                                           ,'description'))
        anno[anno$refseq_mrna=="",  "refseq_mrna"] <- NA
        anno[anno$refseq_ncrna=="", "refseq_ncrna"] <- NA
        save(anno, file=annoFile)
    }else{
        anno <- get(load(annoFile))
    }
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
    if(!(file.exists(annoFile))) {
        bm <- useEnsembl(biomart="ensembl", dataset="scerevisiae_gene_ensembl", host="dec2014.archive.ensembl.org")        
        anno <- getBM(mart=bm, attributes=c('ensembl_gene_id'
                                       ,'ensembl_transcript_id'
                                       ,'external_gene_name'
                                       ,'gene_biotype'
                                       ,'transcript_biotype'
                                       ,'description'))
        save(anno, file=annoFile)
    }else{
        anno <- get(load(annoFile))
    }
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

seqlevels(txdb,force=TRUE) <- seqlevels(txdb)[grep("chrY|chrM",seqlevels(txdb), invert=TRUE)]

## make gnModel to combine with counts for only looking at the browser
## rpkms if calculated with need to use the exons below
gnModel          <- transcriptsBy(txdb, 'gene')
gnModel          <- unlist(gnModel)## Gets the genomic region convered by transcripts
seqinfo(gnModel) <- seqinfo(organism)[seqlevels(gnModel)]

## add gene info to gnModel
if(assembly=="sacCer3"){
    gnModel$gene_id            <- names(gnModel)
    iv                         <- match(gnModel$tx_name, anno$ensembl_transcript_id)
    gnModel$external_gene_name <- anno[iv, "external_gene_name"]
    if( length(gnModel[gnModel$external_gene_name==""]) > 0){
        gnModel[gnModel$external_gene_name==""]$external_gene_name <- NA
    }
    gnModel$gene_biotype       <- anno[iv, "gene_biotype"]
    gnModel$tx_biotype         <- anno[iv, "transcript_biotype"]
    ## filter for protein coding genes
    gnModel                    <- gnModel[gnModel$tx_biotype=="protein_coding" & gnModel$gene_biotype=="protein_coding"]
}else{
    gnModel$gene_id            <- names(gnModel)
    iv                         <- match(gnModel$tx_name, anno$ensembl_transcript_id)
    gnModel$external_gene_name <- anno[iv, "external_gene_name"]
    if( length(gnModel[gnModel$external_gene_name==""]) > 0){
        gnModel[gnModel$external_gene_name==""]$external_gene_name <- NA
    }
    gnModel$gene_biotype       <- anno[iv, "gene_biotype"]
    gnModel$tx_biotype         <- anno[iv, "transcript_biotype"]
    gnModel$refseq_mrna        <- anno[iv, "refseq_mrna"]
    ## filter for only protein coding with refseq mrna id
    gnModel                    <- gnModel[!is.na(gnModel$refseq_mrna) & gnModel$tx_biotype=="protein_coding" & gnModel$gene_biotype=="protein_coding"]
}

###########################################
## pick the mosth highly occupied transcript tss
###########################################
TSS            <- promoters(gnModel, upstream=0, downstream=txTssDown)
TSS$gene_id    <- names(TSS)
names(TSS)     <- NULL

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
elementMetadata(gnModel)[["tssTotalCov"]] <- txTssCov$mean.cov

## sort and take most highly occupied tss
or         <- order( elementMetadata(gnModel)$tssTotalCov,decreasing=TRUE )
gnModel.or <- gnModel[or]
reps       <- duplicated( gnModel.or$gene_id )
gnModel    <- gnModel.or[!reps]

## calculate distance to the nearest gene that does not directly overlap
dists                     <- as.data.frame( distanceToNearest(gnModel, ignore.strand=TRUE) )
gnModel$distanceToNearest <- dists$distance
#gnModel$nearestTx         <- gnModel[dists$subjectHits]$tx_name

###########################################
## Filter for tss's with peaks if provided
###########################################

if ( ! identical( peakFile, character(0)) ){
    peaks                                               <- import.bed(peakFile)
    ol                                                  <- as.data.frame(findOverlaps(promoters(gnModel, upstream=0, downstream=txTssDown)
                                                                                     ,peaks)
                                                                         )
    gnModel$peakOverLap                                 <- 0
    mcols(gnModel[unique(ol$queryHits)])["peakOverLap"] <- 1
    gnModel                                             <- gnModel[gnModel$peakOverLap == 1 ]
    gnModel$peakOverLap                                 <- NULL
}

###########################################
## Filter for minimal length
###########################################
gnModel <- gnModel[width(gnModel) >= minLength]

###########################################
## Filter for minimal distance between genes
###########################################
gnModel <- gnModel[gnModel$distanceToNearest >= minDist]

#gnModel[gnModel$tx_name=="ENST00000529694"] example gene that directly overlaps 

maxCov <- function(bw,model){
    cat("importing:", bw, sep="\n")
    bw.peak <- import.bw(bw,RangedData=FALSE,selection = BigWigSelection(model))
    cat("calc coverage\n")
    bw.peak.cov <- coverage(bw.peak,weight='score')
    cat("get coverage for peak region\n")    
    mean.cov <- with(as.data.frame(model),{
        mcmapply(function(seqname,start,end){
            max(bw.peak.cov[[seqname]][start:end])
        }
       ,mc.cores=Cores
       ,as.character(seqnames),start,end)
    })
    mean.cov <- data.frame(mean.cov)
    rownames(mean.cov) <- model$tx_name
    mean.cov
}

## sum the coverage in the defined promoter region
txMaxCov <- do.call(cbind,mclapply(bwFile
                                  ,model=promoters(gnModel, upstream=0, downstream=txTssDown)
                                  ,maxCov,mc.cores=1)
                    )
elementMetadata(gnModel)[["tssMaxCov"]] <- txMaxCov$mean.cov

gnModel <- gnModel[gnModel$tssMaxCov > 1]

gnModel$tx_id <- NULL

## save the gnModel as a granges object and a tab delimited file
save(gnModel, file=paste0(outName, ".filteredProteinCodingTx.rda"))
write.table(as.data.frame(gnModel)
           ,file=paste0(outName, ".filteredProteinCodingTx.txt")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )

