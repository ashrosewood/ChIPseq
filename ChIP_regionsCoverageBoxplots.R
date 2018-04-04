args <- commandArgs()

help <- function(){
    cat("ChIP_regionsCoverageBoxplots.R :
- From a provided GRnages object or peaks calculate the average coverge and make a boxplot.
- Window is an option if you would like to resize the regions, and the regions are peaks.
- If the regions are Tss or Tes, use upstream and downstream options.\n")
 
    cat("Usage: \n")
    cat("--regions     : GenomicRanges or bed object with the transcripts of interest   [required]\n")
    cat("--type        : Peaks, Tss, or Tes                                             [required]\n")
    cat("--window      : total size of region around types above (use for type = Peaks) [default = 0]\n")
    cat("--upStream    : distance upstream of the region to take (use for Tss and Tes)  [default = 50]\n")
    cat("--downStream  : distance downstream of the region to take (use for Tss and Tes)[default = 50]\n")
    cat("--assembly    : genome assembly build (ex. hg19, dm3)                          [default = hg19]\n")
    cat("--bwFiles     : path to bigWig Files                                           [required]\n")
    cat("--bwPattern   : grep pattern for bigWigs to use quotes (ex. PolII.*293T)       [default = all .bw in path]\n")    
    cat("--numCores    : number of cores to use                                         [default = 10 ]\n")
    cat("--outName     : prefix to your out file names (No .extention)                  [default = basename(bigWigFile) ]\n")
    cat("--cols        : need the same number as samples separated by comma             [default = rainbow]\n")
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
    upStream   <- sub( '--upStream=', '', args[grep('--upStream=', args)] )
    downStream <- sub( '--downStream=', '', args[grep('--downStream=', args)] )
    assembly  <- sub( '--assembly=', '', args[grep('--assembly=', args)] )
    bwFiles   <- sub( '--bwFiles=', '', args[grep('--bwFiles=', args)])
    bwPattern <- sub( '--bwPattern=', '', args[grep('--bwPattern=', args)])
    Cores     <- sub( '--numCores=', '',args[grep('--numCores=',args)])
    outName   <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols    <- sub( '--cols=', '',args[grep('--cols=',args)])
}

outName
bwPattern
bwFiles

bws <- list.files(bwFiles,pattern=".bw", full.names=TRUE)
## if no pattern it will keep them all
bws <- bws[grep(bwPattern,bws,invert=FALSE)]
bws

#setwd("/projects/b1025/arw/analysis/yohhei/quiescent/")
#regions  <- "data_chipseq/Rpb3_3xFLAG_8WG16_120h_Qui_rep1_macsPeaks.bed"
#assembly <- "sacCer3"
#Window   <- 0
#Bins     <- 25
#Type     <- "Peaks"
#bwFiles <- "data_chipseq/"
#bwPattern <- "BY4741A_Qui120hr_8WG16_rep1|set1del_Cps60TAP_Qui120hr_8WG16_rep1|set2del_Qui120hr_8WG16_rep1|set1del_Qui120hr_8WG16_rep1|rad6del_Qui120hr_8WG16_rep1|bre1del_Qui120hr_8WG16_rep1|dot1del_Qui120hr_8WG16_rep1|upb8ubp10del_Qui120hr_8WG16_rep1|Rpb3_3xFLAG_8WG16_120h_Qui_rep1|Rpb3_3xFLAG_8WG16_120h_Qui_rep2"
#Cores <- 10
#outName <-"debug"
    
if (identical(Cores,character(0))){
   Cores <- 10
}else{
    Cores <- as.numeric(Cores)
}

if (identical(assembly,character(0))){
   assembly <- "hg19"
}

if (identical(outName,character(0))){
   outName <- sub(".bw", "", basename(bwFile))
}

if (Type == "Peaks" ){
    print("Type is Peaks")
    if (identical(Window,character(0))){
        print("Window not provided use peak regions")
        Window <- 0
    }else{
        Window <- as.numeric(Window)
        print(paste("Window provided as", Window))
    }
}

if (Type == "Tss" | Type == "Tes" ){
    print(paste("Type is", Type))
    if (identical(upStream,character(0))){
        print("upStream not provided use 50bp")
        upStream <- 50
    }else{
        upStream <- as.numeric(upStream)
        print(paste("upStream provided as", upStream))
    }
    if (identical(downStream,character(0))){
        print("downStream not provided use 50bp")
        downStream <- 50
    }else{
        downStream <- as.numeric(downStream)
        print(paste("downStream provided as", downStream))
    }
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
    if( upStream > 0 | downStream > 0){
        print(paste("resizing Tss regions to upstream", upStream, "and downstream", downStream))
        Model.win <- promoters(Model, upstream=upStream, downstream=downStream)
    }else{
        print(paste("resizing Tss regions to upstream 0 and downstream 1"))
        Model.win <- promoters(Model, upstream=0, downstream=1)
    }    
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".rda"), outName)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Tes" ){
    print("model is Tes")
    tes            <- resize(Model, width=1, fix="end")
    if( upStream > 0 | downStream > 0){
        print(paste("resizing Tes regions to upstream", upStream, "and downstream", downStream))
        Model.win      <- promoters(tes, upstream=upStream, downstream=downStream)
    }else{
        print(paste("resizing Tes regions to upstream 0 and downstream 1"))
        Model.win <- tes
    }
    Model.win$name <- Model.win$gene_id
    fname <- sub("$", paste0("_", Type, "up", upStream, "down", downStream, ".rda"), outName)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}else if( Type=="Peaks" ){
    print("model is Peaks")
    if( Window > 0 ){
        print(paste("resizing peak regions to a window of", Window))
        Model.win      <- resize(Model, width=Window, fix="center")
    }else{
        print("using provided peak regions")
        Model.win <- Model
    }
    fname <- sub("$", paste0("_", Type, Window, ".rda"), outName)
    idx <- GenomicRanges:::get_out_of_bound_index(Model.win)
    if(length(idx) > 0){
        print(paste("remove out of bounds", idx))
        Model.win <- Model.win[-idx]
    }
}

## make sure the output directory exists 
Dir <- dirname(outName)
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

###############################
## calculate average coverage
###############################
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

avgCov <- do.call(cbind,mclapply(bws,model=Model.win,Bin,mc.cores=1))

colnames(avgCov) <- sub("$", paste0("_", Type), colnames(avgCov))

## combine resized regions and average coverage
df <- cbind(as.data.frame(Model.win), signif(avgCov, 4))

## save the gnModel as a granges object and a tab delimited file
write.table(df
           ,file=paste(outName, Type, "AverageCoverages.txt", sep=".")
           ,sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
            )


###############################
## make boxplot
###############################
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(ggsignif)

df.long <- melt(avgCov)
df.long$variable <- sub(paste0("_", Type), "", df.long$variable)

if (identical(cols,character(0))){
    Cols <- rainbow(ncol(avgCov))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}

#minMax <- do.call(rbind, lapply(1:ncol(avgCov), function(x){
#     iQr  <- IQR(avgCov[,x])
#     Ymax <- round(as.numeric(quantile(avgCov[, x])[4] + 1.5 * iQr) *1.05)
#     Ymin <- floor(as.numeric(quantile(avgCov[, x])[2] - 1.5 * iQr) *1.05)
#     cbind(Ymin, Ymax)
#}
#))

minMax <- apply(avgCov, 2, function(x)boxplot.stats(x)$stats[c(1, 5)])

Ymin <- floor(log2(min(minMax[1,])))
Ymax <- ceiling(log2(max(minMax[2,])))

pdf(file=sub("$", paste0(Type, ".boxPlot.pdf"), outName),width=6,height=5)
print({
    p <-
        ggplot(df.long, aes(x=variable, y=log2(value), fill=variable)) + 
        stat_boxplot(geom ='errorbar')+
        geom_boxplot(notch=TRUE,outlier.colour=NA) +
        theme_classic() +
        scale_fill_manual("sample", values=Cols)+
        coord_cartesian(ylim = c(Ymin, Ymax)) +
        ylab("log2( Average Coverage ) ") + 
        xlab("") +
        ggtitle(Type) +
        theme(text = element_text(size=8),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color="black",size=8),
              plot.title=element_text(size=8))
})
dev.off()









