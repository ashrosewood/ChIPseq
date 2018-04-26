args <- commandArgs()

help <- function(){
    cat("ChIP_makeMetaPlot.R :
- Plot ecdf of pausing indexes calculated from ChIP_CalcPausingIndex.R and save as pdf\n")
    cat("Usage: \n")
    cat("--Dir     : path to heatmap tables                                          [required]\n")
    cat("--Pattern : grep pattern for heatmap tables to use quotes (ex. PolII.*293T) [optional]\n")    
    cat("--outName : prefix for plot name. Type will be appended to this in script   [required]\n")
    cat("--cols    : need the same number as samples separated by comma              [default = viridis palette]\n")
    cat("--bins    : number of bins each window is (ie 25 or 50)                     [default = 25]\n")
    cat("--Type    : Tss, Tes, Peaks                                                 [default = Tss]\n")
    cat("--Height  : plot height                                                     [default = max of colMeans]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir  <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    outName <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols    <- sub( '--cols=', '',args[grep('--cols=',args)])
    bins    <- sub( '--bins=', '',args[grep('--bins=',args)])
    Type    <- sub( '--Type=', '',args[grep('--Type=',args)])
    Height  <- sub( '--Height=', '',args[grep('--Height=',args)])
}

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

#setwd("/projects/b1025/arw/analysis/kevin/SEC/")
#Dir="tables/heatmaps/MYC_H2171_1021"
#Pattern="H2171_1021.*_Peaks5000.binary"
#outName="plots/meta_plots/H2171_1021_Peaks5000binary"
#cols="plots/plotColors/1021_colors.txt"
#bins=25
#Type="Peaks"

## set defaults
if (identical(bins,character(0))){
   bins <- 25
}else{
    bins <- as.numeric(bins)
}

if (identical(Type,character(0))){
   Type <- "Tss"
}

Dir
Pattern
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo <- foo[grep(Type,foo,invert=FALSE)]
foo

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub("_Tss.*|_Tes.*|_Peak.*",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
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

head(df)

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
    Cols <- rainbow(length(foo))#viridis(length(foo))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}

pdf(file=sub("$", paste0(Type, ".metaPlot.pdf"), outName),width=8,height=5)
print({
    p <-

    ggplot(df.long, aes(x=x, y=value, color=variable)) +
    geom_line(size=1.5)+
    scale_color_manual("sample", values = Cols )+
        scale_x_continuous(breaks=c(0
                                   ,nrow(df)/4
                                   ,nrow(df)/2
                                   ,nrow(df)/4*3
                                   ,nrow(df)
                                    )
                          ,labels=c(paste0("-",nrow(df)/2*bins)
                                   ,paste0("-",nrow(df)/4*bins)
                                   ,Type
                                   ,paste0("+",nrow(df)/4*bins)
                                   ,paste0("+", nrow(df)/2*bins))
                           )+
    ylab("r.p.m.")+
    xlab("")+
        ylim(Min,Height)+
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
