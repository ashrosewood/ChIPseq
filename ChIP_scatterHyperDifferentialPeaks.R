args <- commandArgs()

help <- function(){
    cat("ChIP_scatterHyperDifferentialPeaks.R :
- Make an ma plot and scatter plot from the *hyper.txt output from ChIP_callDifferentialPeaksHyper.R.
- It will calculated the number of reads that overlap the entire peak regions in the experiment and control.\n")
    cat("Usage: \n")
    cat("--dataTable  : Full path to hyper.txt                                             [required]\n")
    cat("--Control    : Name of control sample in table (must match table column names)    [required]\n")
    cat("--Experiment : Name of experiment sample in table (must match table column names) [required]\n")
    cat("--adjp       : adjusted p-value cut off for differential peaks                    [default = 0.01]\n")
    cat("--FC         : adjusted p-value cut off for differential peaks                    [default = 1.5]\n")
    cat("--upCol      : color of up-regulated peaks                                        [maroon]\n")
    cat("--downCol    : color of down-regulated peaks                                      [dodgerblue]\n")
    cat("--ncCol      : color of non-changing peaks                                        [gray87]\n")
    cat("--calcCorr   : add pearson correlation to the scatter plot (0/1)                  [default = 0]\n")
    cat("--outName    : Path to output file (no file extension ie .pdf)                    [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    dataTable   <- sub('--dataTable=', '', args[grep('--dataTable=', args)])
    Con         <- sub( '--Control=', '', args[grep('--Control=', args)] )
    Exp         <- sub( '--Experiment=', '', args[grep('--Experiment=', args)] )
    adjp        <- sub('--adjp=', '', args[grep('--adjp=', args)])
    FC          <- sub('--FC=', '', args[grep('--FC=', args)])
    upCol       <- sub('--upCol=', '', args[grep('--upCol=', args)])
    downCol     <- sub('--downCol=', '', args[grep('--downCol=', args)])
    ncCol       <- sub('--ncCol=', '', args[grep('--ncCol=', args)])
    outName     <- sub('--outName=', '',args[grep('--outName=',args)] )
    calcCorr    <- sub('--calcCorr=', '',args[grep('--calcCorr=',args)] )
}

if (identical(calcCorr,character(0))){
    calcCorr <- 0
}else{
    calcCorr <- as.numeric(calcCorr)
}

if (identical(adjp,character(0))){
    adjp <- 0.01
}else{
    adjp <- as.numeric(adjp)
}

if (identical(FC,character(0))){
   FC <- 1.5
}else{
    FC <- as.numeric(FC)
}

if (identical(upCol,character(0))){
    upCol <- "maroon"
}
if (identical(ncCol,character(0))){
    ncCol <- "grey87"
}
if (identical(downCol,character(0))){
    downCol <- "dodgerblue"
}

## for debugging
#setwd("/projects/b1025/arw/analysis/fei/paf1/")
#dataTable <- "tables/peak_calling/H3K4me3_1D_vs_0h_Aux_PAF1AID_DLD1_1015.phyper.txt"
#upCol="#B2182B60"
#downCol="#2166AC60"
#ncCol="#A0A0A060"
#Con="H3K4me3_0h_Aux_PAF1AID_DLD1_1015"
#Exp="H3K4me3_1D_Aux_PAF1AID_DLD1_1015"
#outName <- "plots/chipseq/H3K4me3_1D_vs_0h_Aux_PAF1AID_DLD1_1015"

## Column names with hypens are automatically changed to periods
## check for hyphens in Con and Exp
if(length(grep("", Con))==1){
    Con <- gsub("-",".",Con)
}

if(length(grep("", Exp))==1){
    Exp <- gsub("-",".",Exp)
}

## make plot directory
Dir <- dirname(outName)
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

scatterPlots <- function(dataTable){
    ## read in table using existing column names
    df <- read.table(dataTable, sep='\t', header=TRUE)
    ##-----------------Corr plot---------------##
    down <- df$p.under.adj < adjp & df$log2FC <= -log2(FC)
    up <- df$p.over.adj < adjp & df$log2FC >= log2(FC)
    ## set axis min and max
    MINs <- min(apply(df[,grep(".rpkm$", names(df))],2,  function(x) floor(log10(min(x[x>0])))))
    MAXs <- max(apply(df[,grep(".rpkm$", names(df))],2,  function(x) ceiling(log10(max(x)))))
    ## make scatter plot
    pdf(paste(outName, ".scatterPlot.pdf", sep=""),height=5,width=5)
    par(mar=c(5,5,5,5))
    plot(x=log10(df[,paste(Con,"rpkm", sep=".")])
          ,y=log10(df[,paste(Exp,"rpkm", sep=".")])
          ,col=ncCol
          ,pch=19
          ,cex=1
          ,xlim=c(-1, 2.5)
          ,ylim=c(-1, 2.5)
          ,main=gsub("_rep[0-9]|_Rep[0-9]", "",gsub(".all", "", paste(Exp, Con, sep="\nvs\n")))
          ,ylab=paste0("log2(", sub(".all", "", Exp), ")") 
          ,xlab=paste0("log2(", sub(".all", "", Con), ")") 
        ,cex.lab=1, cex.axis=1.1, cex.main=1, cex.sub=1)
    ## highlight dge
    points(log10(df[up,paste(Con,"rpkm", sep=".")]), log10(df[up,paste(Exp,"rpkm", sep=".")]), col=upCol, pch=19,cex=1)
    points(log10(df[down,paste(Con,"rpkm", sep=".")]), log10(df[down,paste(Exp,"rpkm", sep=".")]), col=downCol, pch=19,cex=1)
    ## total dge with cutoffs
    legend("bottomright", c(paste(sum(up), "FC >=", FC)
                           ,paste(sum(down), "FC <=", FC)
                           ,paste(nrow(df)-sum(up)-sum(down), "no change")
                            ),
           pch=c(19,19,19), col=c(upCol,downCol,ncCol),cex=0.8,title=paste("adjp < ",adjp))
    ## add corr values
    if(calcCorr > 0){
        legend("topleft", title="PCC", legend= round(cor(df[,paste(Con,"cpm", sep=".")],df[,paste(Exp,"cpm", sep=".")]), 4), cex=1)
    }
    abline(a=0,b=1, lty=2, col="grey70")
    dev.off()    
    ## make ma plot
    df$rpkm.avg <- log2(apply(df[,grep("rpkm", names(df))],1,mean))
    pdf(paste(outName, ".maPlot.pdf", sep=""),height=5,width=5)
    par(mar=c(5,5,5,5))
    plot(x=df$rpkm.avg
        ,y=df$log2FC
          ,col=ncCol
          ,pch=19
          ,cex=1
         # ,xlim=c(MINs, MAXs)
         # ,ylim=c(MINs, MAXs)
          ,main=gsub("_rep[0-9]|_Rep[0-9]", "",gsub(".all", "", paste(Exp, Con, sep="\nvs\n")))
          ,ylab=paste0("log2(", sub(".all", "", Exp), ")") 
          ,xlab=paste0("log2(", sub(".all", "", Con), ")") 
        ,cex.lab=1, cex.axis=1.1, cex.main=1, cex.sub=1)
    ## highlight dge
    points(df[up,"rpkm.avg"], df[up,"log2FC"], col=upCol, pch=19,cex=1)
    points(df[down,"rpkm.avg"], df[down,"log2FC"], col=downCol, pch=19,cex=1)
    ## total dge with cutoffs
    legend("bottomright", c(paste(sum(up), "FC >=", FC)
                            ,paste(sum(down), "FC <=", FC)
                           ,paste(nrow(df)-sum(up)-sum(down), "no change")
                            ),
           pch=c(19,19,19), col=c(upCol,downCol,ncCol),cex=0.8,title=paste("adjp < ",adjp))
    abline(h=0, lty=2, col="grey70")
    dev.off()    
}

scatterPlots(dataTable=dataTable)
