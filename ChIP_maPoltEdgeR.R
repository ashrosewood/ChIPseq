args <- commandArgs()

help <- function(){
    cat("ChIP_maPlotEdgeR.R :
- For edgeR output from ChIP_differentialPeaksEdgeR.R make MA plot of differentially expressed regions.
- Output file will be the the same name as the degFile but with a file extension MAplot.pdf.
- Color options can be hex or rcolors\n")
    cat("Usage: \n")
    cat("--degFile : edgeR table                 [ required ]\n")
    cat("--adjp    : FDR adjusted p-value cutoff [ default = 0.01 ]\n")
    cat("--FC      : fold change cutoff          [ default = 2 ]\n")
    cat("--upCol   : up regulated genes color    [ default = purple ]\n")
    cat("--downCol : down regulated genes color  [ default = green ]\n")    
    cat("--ncCol   : non-changing genes color    [ default = grey ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    degFile  <-sub('--degFile=', '', args[grep('--degFile=', args)])
    adjp     <-sub('--adjp=', '', args[grep('--adjp=', args)])
    FC       <-sub('--FC=', '', args[grep('--FC=', args)])
    upCol    <- sub('--upCol=', '', args[grep('--upCol=', args)])
    downCol  <- sub('--downCol=', '', args[grep('--downCol=', args)])
    ncCol    <- sub('--ncCol=', '', args[grep('--ncCol=', args)])
}

if (identical(adjp,character(0))){
   adjp<-0.01
}else{
    adjp <- as.numeric(adjp)
}
if (identical(adjp,character(0))){
   labelTop<-1
}
if (identical(FC,character(0))){
   FC <- 1
}else{
    FC <- as.numeric(FC)
}
if (identical(upCol,character(0))){
   upCol <- "purple"
}
if (identical(downCol,character(0))){
   downCol <- "green"
}
if (identical(ncCol,character(0))){
   ncCol <- "grey"
}

library(edgeR)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('rda',degFile)){
    deg <- get(load(file=degFile))
        deg <- deg[order(deg$adj.p),]
}
if(grepl('txt',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")
    rownames(deg)<-deg$X
    deg <- deg[order(deg$adj.p),]
}
head(deg)
dim(deg)
sampleNum <- (dim(deg)[2] - 15)/2

up <- deg$adj.p < adjp & deg$log2FC > log2(FC)
sum(up)

down <- deg$adj.p < adjp & deg$log2FC < -log2(FC)
sum(down)

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.edgeR.txt$|\\.edgeR.df.rda","",degFile)

pdfFile <- paste(comparison,adjplabel,"MAplot.pdf",sep=".")
print(pdfFile)
comparison <- gsub("^.*/","",comparison)

pdf(pdfFile,height=6,width=6)
plot(deg$log2CPM,deg$log2FC
    ,col=ncCol
    ,pch=19
    ,cex=1
    ,main=comparison,xlab="log2(CPM)",ylab="log2FC"
    ,cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(deg[up,"log2CPM"],deg[up,"log2FC"], col=upCol, pch=19,cex=1)
points(deg[down,"log2CPM"],deg[down,"log2FC"], col=downCol, pch=19,cex=1)
legend("topleft", title=paste("FC:", FC, "adjp:", adjp)
      ,legend=c(paste(sum(up))
               ,paste(sum(down))
               ,paste(nrow(deg) - sum(up) - sum(down))
                )
      ,pch=c(19,19,19), col=c(upCol,downCol,ncCol),cex=0.9,)
abline(h=0, lty="dashed", col="grey")
dev.off()
