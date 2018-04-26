args <- commandArgs()

help <- function(){
    cat("extendReads.R :
- Extend reads of a bam file and save to bed.
You have the option to save separate files for each strands and to flip the strand.
If you only want to separate strands without extention provide an extLen of zero and sepStrands of 1
By default the unassembled contigs and chrM are not kept.\n")
    cat("Usage: \n")
    cat("--bamFile    : Sample bam file (must be a full path)                         [required]\n")    
    cat("--outName    : Prefix for your output file                                   [default = bamFile full prefix]\n")    
    cat("--assembly   : Genome build hg19, mm9, mm10, dm3 ect.                        [default = hg19]\n")
    cat("--extLen     : total number of bases to extend reads to                      [required]\n")
    cat("--flipStrand : If the strand should be flipped. Depends where the primer is. [default = 1]\n")
    cat("--sepStrands : If the strands should be saved to separate files (0/1)        [default = 1]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    bamFile     <- sub('--bamFile=', '', args[grep('--bamFile=', args)])
    outName     <- sub('--outName=', '', args[grep('--outName=', args)])
    assembly    <- sub('--assembly=', '', args[grep('--assembly=', args)])
    extLen      <- sub('--extLen=', '', args[grep('--extLen=', args)])
    sepStrands  <- sub('--sepStrands=', '', args[grep('--sepStrands=', args)])
    flipStrand   <- sub('--flipStrand=', '', args[grep('--flipStrand=', args)])
}

##bamFile <- list.files("data_ttseq", pattern="DMSO.*.bam", full.names=TRUE)[1]
##outName <- paste0("/projects/b1042/Shilatifard/arw/sicer.test/", sub(".bam", "", basename(bamFile)))
#extLen <- 150
#sepStrands <- 1

print(bamFile)

if (identical(outName,character(0))){
   outName <- sub(".bam", "", bamFile)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(sepStrands,character(0))){
   sepStrands <- 1
}else{
   sepStrands <- as.numeric(sepStrands)
}

if (identical(flipStrand,character(0))){
   flipStrand <- 1
}else{
    flipStrand <- as.numeric(flipStrand)
}

extLen <- as.numeric(extLen)

library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

print(assembly)
if ((assembly == "hg19") || (assembly == "hg38")) { organismStr <- "Hsapiens" }
if ((assembly == "mm9") || (assembly == "mm10")) { organismStr <- "Mmusculus" }
if (assembly == "sacCer3") organismStr <- "Scerevisiae"
if (assembly == "dm3") organismStr <- "Dmelanogaster"
print(organismStr)

assemblyLibrary <- paste("BSgenome.", organismStr, ".UCSC.", assembly, sep="")
print(assemblyLibrary)

library(assemblyLibrary,character.only=TRUE)

if ((assembly == "hg19") || (assembly == "hg38")) { organism <- Hsapiens }
if ((assembly == "mm9") || (assembly == "mm10")) { organism <- Mmusculus }
if (assembly == "sacCer3") organism <- Scerevisiae
if (assembly == "dm3") organism <-Dmelanogaster

bam2bw <- function(BF,organism){
    
    cat("opening:", BF, sep="\n")
    bd                       <- readGAlignments(BF, use.names=TRUE)
    seqlevels(bd,force=TRUE) <- seqlevels(bd)[grep("_|chrM",seqlevels(bd), invert=TRUE)]
    print(seqlevels(bd))

    
    
    if( flipStrand > 0 ){
        cat("change the strand:", BF, sep="\n")
        strand(bd) <- ifelse(strand(bd) == '+', '-', '+')
    }

    cat("convert to GRanges\n")
    mygr <- as(bd,"GRanges")
    if (extLen > 0){
        cat("extend reads:", BF, sep="\n")
        mygr <- resize(mygr, extLen, fix='start')
    }
    
    if (sepStrands > 0){
        
        cat("save separate strands\n")
        ## get plus coverage                                                             
        plus                  <- mygr[strand(mygr) == "+"]
        ## get minus coverage
        minus                 <- mygr[strand(mygr) == "-"]
        ## set the outfile name
        plus.outfile  <- sub("$", ".plus.bed", outName)
        minus.outfile <- sub("$", ".minus.bed", outName)

        ## export rpm to bigWig
        cat(paste("exporting to plus bed", plus.outfile, "\n", sep="\t"))
        export.bed(plus, plus.outfile)
        cat(paste("exporting to minus bed", minus.outfile, "\n", sep="\t"))
        export.bed(minus, minus.outfile)      
        cat("export complete:", BF, sep="\n")

    }else{

        cat("save both strands in one file\n")
        outfile <- sub("$", ".bed", outName)
        cat(paste("exporting to bed", outfile, "\n", sep="\t"))
        export.bed(mygr, outfile)
        cat("export complete:", BF, sep="\n")
        
    }
}

# for each element of our vector, call the bam2bw function
mclapply(bamFile,organism=organism,bam2bw,mc.cores=1,mc.preschedule=FALSE)
