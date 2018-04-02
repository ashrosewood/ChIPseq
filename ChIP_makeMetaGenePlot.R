args <- commandArgs()

help <- function(){
    cat("ChIP_makeMetaGenePlot.R :
- Make meta gene plot from metaGeneMatrix.R and save as pdf\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                          [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T) [optional]\n")    
    cat("--outName    : prefix for plot name. Type will be appended to this in script   [required]\n")
    cat("--cols       : need the same number as samples separated by comma              [default = viridis palette]\n")
    cat("--upStream   : number of bp upstream of Tss in matrix                          [required]\n")
    cat("--downStream : number of bp downstream of Tss in matrix                        [required]\n")
    cat("--Height     : plot height if you do not want to be set the max                [default = max of colMeans]\n")
    cat("--Control    : name of control to compare all samples to                       [default = max of colMeans]
                         if set makes a FC plot over this sample.\n")
    cat("--GeneList   : text file gene ids (rownames) to subsample by                   [default = none]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir  <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern    <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    outName    <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols       <- sub( '--cols=', '',args[grep('--cols=',args)])
    upStream   <- sub( '--upStream=', '',args[grep('--upStream=',args)])
    downStream <- sub( '--downStream=', '',args[grep('--downStream=',args)])
    Height     <- sub( '--Height=', '',args[grep('--Height=',args)])
    Control    <- sub( '--Control=', '',args[grep('--Control=',args)])
    GeneList   <- sub( '--GeneList=', '',args[grep('--GeneList=',args)])
}

if (identical(Control, character(0))){
    FcPlot <- 0
}else{
    FcPlot <- 1
}

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

Dir
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo

#if ( !identical(GeneList, character(0)) ){
#    print(paste("Filter for genes in:", GeneList))
#    genes <- read.table(GeneList, sep="\t", header=FALSE)
#}


## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
    colnames(df) <- paste(sub(".df", "", oname), 1:ncol(df), sep=".")
 #   if (exists("GeneList")){
 #       print(paste("Filter for genes in:", GeneList, foo[i]))
 #       df <- df[rownames(df) %in% paste(genes[,1]),]
 #   }
    assign(oname, df)
}

SAMPLES <- ls(pattern=".df$")
SAMPLES

## make a data frame with matrix column averages
df <- do.call(cbind, lapply(SAMPLES, function(x){
    avg           <- as.data.frame( colMeans(get(x) ))
    colnames(avg) <- sub(".df", "", x)
    avg
  }
))
df$x <- seq(1,nrow(df))
df.long <- melt(df, id.vars = "x")

mround <- function(x,base){ 
    base*ceiling(x/base) 
} 

if( max(df[,grep("^x$",names(df),invert=TRUE)]) > 0.5 ){
    Ymax <- mround(ceiling(max(df[,grep("^x$",names(df),invert=TRUE)])),2)
    Min <- 0
}else{
    Ymax <- mround(max(df[,grep("^x$",names(df),invert=TRUE)]),0.05)
    Min <- mround(min(df[,grep("^x$",names(df),invert=TRUE)]),-0.05)
}

if (identical(Height, character(0))){
    Height <- Ymax
}else{
    Height <- as.numeric(Height)
}

if (identical(cols,character(0))){
    Cols <- viridis(length(foo))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}


pdf(file=sub("$", ".metaGenePlot.pdf", outName),width=8,height=5)
print({
    p <-
    ggplot(df.long, aes(x=x, y=value, color=variable)) +
    geom_line(size=0.8)+
    scale_color_manual("sample", values = Cols )+
    scale_x_continuous(breaks=c(0
                               ,100
                               ,100+(400/3)
                               ,100+(400/3*2)
                               ,100+(400/3*3)
                               ,600
                                )
                      ,labels=c(paste0("-",upStream)
                               ,"TSS"
                               ,"33%"
                               ,"66%"
                               ,"TES"
                               ,paste0("+", downStream)
                                ))+
    ylab("r.p.m.")+
    xlab("")+
    ylim(0,Height)+
    ggtitle(basename(outName))+
    theme(panel.grid.major = element_blank()
         ,panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key=element_blank())+
    theme(text = element_text(size=12),
          axis.text.x = element_text(vjust=1,color="black",size=12),
          axis.text.y = element_text(color="black",size=12),
          plot.title=element_text(size=12))

})
dev.off()


if(FcPlot==1){
    print(paste("Making Fold Change Plot vs Control:",Control))
    exps <- SAMPLES[SAMPLES!=paste0(Control, ".df")]
    fc <- do.call(cbind, lapply(exps, function(x){
        exp         <- sub(".df", "", x)
        con         <- sub(".df", "", SAMPLES[SAMPLES==paste0(Control, ".df")])
        FC            <- log2(df[,exp]/df[,con])
        avg           <- as.data.frame( FC )
        colnames(avg) <- sub(".df", "", x)
        avg
    }
    ))
    fc$x <- seq(1,nrow(fc))
    df.long <- melt(fc, id.vars = "x")
    Max <- mround(max(abs(fc[,grep("x", colnames(fc), invert=TRUE)])), 0.5)
    if ( min(fc[,grep("x", colnames(fc), invert=TRUE)]) > 0 ){
        Min <- 0
    }else{
        Min <- -Max
    }
    pdf(file=sub("$", ".metaGeneFC.pdf", outName),width=8,height=5)
    print({
        p <-
            ggplot(df.long, aes(x=x, y=value, color=variable)) +
            geom_line(size=1)+
            scale_color_manual("sample", values = Cols )+
            scale_x_continuous(breaks=c(0
                                       ,100
                                       ,100+(400/3)
                                       ,100+(400/3*2)
                                       ,100+(400/3*3)
                                       ,600
                                        )
                              ,labels=c(paste0("-",upStream)
                                       ,"TSS"
                                       ,"33%"
                                       ,"66%"
                                       ,"TES"
                                       ,paste0("+", downStream)
                                        ))+
            ylab(paste0("log2(FC over ", Control, ")"))+
            xlab("")+
            ylim(Min,Max)+
            ggtitle(basename(outName))+
            theme(panel.grid.major = element_blank()
                 ,panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key=element_blank())+
            theme(text = element_text(size=12),
                  axis.text.x = element_text(vjust=1,color="black",size=12),
                  axis.text.y = element_text(color="black",size=12),
                  plot.title=element_text(size=12))
    })
    dev.off()
}
