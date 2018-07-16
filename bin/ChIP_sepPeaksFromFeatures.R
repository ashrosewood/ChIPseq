args <- commandArgs()

help <- function(){
    cat("ChIP_sepPeaksFromFeatures.R :
- Separate peaks into overlapping and non overlapping and save each in a bed file.
- Note only considers protein coding with refseq_mrna identifier.\n")
    cat("Usage: \n")
    cat("--txdbFile    : transcript data base (.txdb)              [required]\n")
    cat("--assembly    : genome (hg19, mm9, mm10, dm3, sacCer3)    [default = hg19]\n")
    cat("--Type        : Type of feature (gene or Tss)             [default = Tss]\n")
    cat("--upStream    : distance upstream of the regions          [default = 0]\n")
    cat("--downStream  : distance downstream of the regions        [default = 0]\n")
    cat("--peakFile    : bed file                                  [required]\n")
    cat("--outName     : prefix to your out file names             [default = peakFile%.bed ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if( !is.na(charmatch("-help",args)) || !is.na(charmatch("-h",args)) ){
    help()
} else {
    txdbFile   <- sub( '--txdbFile=', '', args[grep('--txdbFile=', args)] )
    assembly   <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    Type       <- sub( '--Type=', '', args[grep('--Type=', args)])
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)])
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)])
    peakFile   <- sub( '--peakFile=', '', args[grep('--peakFile=', args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
}

if (identical(upStream,character(0))){
    print("upStream not provided use 0bp")
    upStream <- 0
}else{
    upStream <- as.numeric(upStream)
    print(paste("upStream provided as", upStream))
}
if (identical(downStream,character(0))){
    print("downStream not provided use 0bp")
    downStream <- 0
}else{
    downStream <- as.numeric(downStream)
    print(paste("downStream provided as", downStream))
}

if (identical(Type,character(0))){
    print("Type not provided use Tss")
    Type <- "Tss"
}else{
    print(paste("Type provided as", Type))
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bed", "", peakFile)
}

## make output directory if does not exist
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

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
    annoFile <- "/projects/b1025/anno/biomaRt/hg19.Ens_75.biomaRt.geneAnno.Rdata"
    anno     <- get(load(annoFile))
}
if (assembly == "mm9") {
    organism <- Mmusculus
    annoFile <- "/projects/b1025/anno/biomaRt/mm9.Ens_67.biomaRt.geneAnno.Rdata"
    anno     <- get(load(annoFile))    
}
if (assembly == "mm10") {
    organism <- Mmusculus
    annoFile <- "/projects/b1025/anno/biomaRt/mm10.Ens_78.biomaRt.geneAnno.Rdata"
    anno     <- get(load(annoFile))    
}
if (assembly == "sacCer3") {
    organism <- Scerevisiae
    annoFile <- "/projects/b1025/anno/biomaRt/sacSer3.Ens_78.biomaRt.geneAnno.Rdata"
    anno     <- get(load(annoFile))
}
if (assembly == "dm3") {
    organism <- Dmelanogaster
    annoFile <- "/projects/b1025/anno/biomaRt/dm3.Ens_74.biomaRt.geneAnno.Rdata"
    anno     <- get(load(annoFile))
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
    gnModel                   <- gnModel[!is.na(gnModel$refseq_mrna) & !gnModel$refseq_mrna =="NA" & gnModel$tx_biotype=="protein_coding" & gnModel$gene_biotype=="protein_coding"]
}


if(Type=="Tss"){
    if(upStream > 0 | downStream > 0){
        model     <- promoters(gnModel, upstream=upStream, downstream=downStream)
    }else{
        model     <- promoters(gnModel, upstream=upStream, downstream=1)
    }
}else{
    if(upStream > 0 | downStream > 0){
        Tss   <- promoters(gnModel, upstream = upStream, downstream = 0 )
        model <- resize(Tss, fix="start", width = (width(gnModel) + upStream + downStream ))
    }else{
        model <- gnModel
    }
}

## find overlapping peaks
peaks     <- import.bed(peakFile)
ol        <- as.data.frame(findOverlaps(peaks,model))

## seperate
peaks.ol  <- as.data.frame(peaks[unique(ol$queryHits)])
peaks.non <- as.data.frame(peaks[-unique(ol$queryHits)])

## fix the start for 0 based
peaks.ol$start <- peaks.ol$start - 1
peaks.non$start <- peaks.non$start - 1

peaks.ol  <- peaks.ol[,c("seqnames", "start", "end", "name", "score")]
peaks.non <- peaks.non[,c("seqnames", "start", "end", "name", "score")]

## save files
write.table(peaks.ol, file=sub("$", paste0(".", Type, ".overlap.bed"), outName)
           ,row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

write.table(peaks.non, file=sub("$", paste0(".", Type, ".nonoverlap.bed"), outName)
           ,row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

print("done")
