args <- commandArgs()

help <- function(){
    cat("ChIP_makePeakVenns.R :
- For a directory of peak bed files, plot the overlap for each possible comparison.
- By default it looks for all .bed files, so if you have other types of bed files, use peakPattern.\n")
    cat("Usage: \n")
    cat("--peakDir     : path to peak bed files                                         [required]\n")
    cat("--peakPattern : grep pattern for (ex. PolII.*rep1)                             [optional]\n")    
    cat("--outDir      : path to save pdfs                                              [required]\n")
    cat("--colA        : color of the first circle                                      [default = red]\n")
    cat("--colB        : color of the second circle                                     [default = blue]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    peakDir     <- sub( '--peakDir=', '', args[grep('--peakDir=', args)])
    peakPattern <- sub( '--peakPattern=', '', args[grep('--peakPattern=', args)])
    outDir      <- sub( '--outDir=', '',args[grep('--outDir=',args)])
    colA        <- sub( '--colA=', '',args[grep('--colA=',args)])
    colB        <- sub( '--colB=', '',args[grep('--colB=',args)])
}

if (identical(colA,character(0))){
    colA <- "red"
}

if (identical(colA,character(0))){
    colB <- "blue"
}


if(!(file.exists( outDir ))) {
    print(paste("mkdir", outDir))
    dir.create(outDir,FALSE,TRUE)  
}

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(VennDiagram)

peaks <- list.files(peakDir, pattern=".bed", full.names=TRUE)
peaks <- peaks[grep(peakPattern, peaks, invert=FALSE)]
peaks

for (i in 1:length(peaks))
{
    oname = gsub(".macsPeaks.bed|.peaks.bed",".peak",basename(peaks[i]))
    oname <- gsub("-","_",oname)
    Peak  <- import.bed(peaks[i])
    assign(oname, Peak)
}

Peaks <- ls(pattern=".peak$")
print(Peaks)


for (i in 1:length(Peaks)) {
    for (j in 1:length(Peaks)){        
        if(i!=j){
            DrawVenn <- function(Con, Exp){
                print(paste("Peak file 1:", Con))
                print(paste("Peak file 2:", Exp))
                Con.name                    <- sub(".peak", "", Con)
                Exp.name                    <- sub(".peak", "", Exp)
                gr1                         <- get(Con)
                gr2                         <- get(Exp)
                seqlevels(gr1,force=TRUE)   <- seqlevels(gr1)[grep("chrY|chrM",seqlevels(gr1), invert=TRUE)]
                seqlevels(gr2,force=TRUE)   <- seqlevels(gr2)[grep("chrY|chrM",seqlevels(gr2), invert=TRUE)]
                a1                          <- length(gr1)
                a2                          <- length(gr2)
                a1
                a2
                ## find overlaps of at least 1
                ol12                        <- as.data.frame(findOverlaps(gr1,gr2))
                ## take the highest number of overlaps
                if( length(unique(ol12$queryHits)) > length(unique(ol12$subjectHits)) ){
                    a12                     <- length(unique(ol12$queryHits))
                }else{
                    a12                     <- length(unique(ol12$subjectHits))
                }
                if( a12 > a2 ){
                    a12                     <- a2
                }
                if( a12 > a1 ){
                    a12                     <- a1
                }
                a12
                pdf(paste(outDir, paste0(paste(Con, "vs", Exp, "peakOverlap", sep="_"), ".pdf"), sep="/"),width=5,height=5)
                draw.pairwise.venn(area1=a1,area2=a2,cross.area=a12,
                                   category = c(Con.name, Exp.name),
                                   fill = c(colA, colB),
                                   lty = "solid",#"blank",
                                   lwd=2,
                                   cex = 2,
                                   cat.cex = 2,
                                   #cat.pos=c(-7,7),
                                   cat.pos = c(180, 0),
                                   cat.dist=rep(0.09,2),
                                   rotation.degree=90,
                                   alpha=c(0.6,0.6))                
                dev.off()
            }
            DrawVenn(Con=Peaks[i], Exp=Peaks[j])
        }
    }
}

print("done")
