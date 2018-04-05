args <- commandArgs()

help <- function(){
    cat("ChIP_genomicDistribution.R :
- From a provided GRnages object or peaks calculate the 
- Window is an option if you would like to resize the regions, and the regions are peaks.
- If the regions are Tss or Tes, use upstream and downstream options.\n")
 
    cat("Usage: \n")
    cat("--regions     : Summits bed file where the start and end are the summit position [required]\n")
    cat("--anno        : Bed file with annotation for the entire genome                   [required]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                            [default = hg19]\n")
    cat("--outName     : prefix to your out file names (No .extention)                    [required]\n")
    cat("--cols        : need the same number as samples separated by comma               [default = rainbow]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    regions   <- sub( '--regions=', '', args[grep('--regions=', args)] )
    anno      <- sub( '--anno=', '', args[grep('--anno=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols    <- sub( '--cols=', '',args[grep('--cols=',args)])
}

outName
regions
anno

#setwd("/projects/b1025/arw/analysis/yohhei/quiescent/")
#regions  <- "/projects/b1025/tango/TANGO-984.985.989/TANGO-989/peaks/BY4741A_Qui120hr_8WG16_121117_summits.bed"
#anno <- "/projects/b1025/arw/anno/sacCer3/genome_annotation.bed.gz"
#outName <-"debug"
#assembly <- "sacCer3"

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
library(plyr)
library(RColorBrewer)

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

Peaks  <- import.bed(regions)
Genome <- import.bed(anno)

ol <- as.data.frame(findOverlaps(Peaks, Genome))

ol$name <- Genome[ol$subjectHits]$name

counts <- count(ol, "name")
counts$percent <- counts$freq / nrow(ol) * 100
counts

if (identical(cols,character(0))){
    Cols <- rainbow(nrow(counts))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}

#rev(brewer.pal(nrow(counts), "RdYlBu") )
#colorRampPalette(c("gold", "darkorange", "maroon", "deeppink4", "darkorchid4"))(500)[seq(1,500,length.out=nrow(counts))]

## plot
pdf(file=paste0(outName, ".genomicDistributionPie.pdf"),width=6,height=6)
pie(counts$percent,col=Cols
   ,labels=sub("$","%",(signif(counts$percent,digits=3)))
   ,main=""
   ,border=TRUE
   ,cex=1.1,cex.main=1.1)
points(0,0,col="black",pch=21,cex=28, bg = "white")
legend("center", counts$name, fill=Cols, cex=1.1,box.col=NA)
text(0,-1,paste(length(Peaks),"Total Peaks",sep=" "), cex=1.12)
text(0,1,"Peaks Distribution", cex=1.12)
dev.off()
