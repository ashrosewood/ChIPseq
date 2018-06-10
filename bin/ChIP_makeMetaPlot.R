args <- commandArgs()

help <- function(){
    cat("ChIP_makeMetaPlot.R :
- Plot ecdf of pausing indexes calculated from ChIP_CalcPausingIndex.R and save as pdf\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                          [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T) [optional]\n")    
    cat("--outName    : prefix for plot name. Type will be appended to this in script   [required]\n")
    cat("--cols       : need the same number as samples separated by comma              [default = viridis palette]\n")
    cat("--bins       : number of bins each window is (ie 25 or 50)                     [default = 25]\n")
    cat("--Type       : Tss, Tes, Peaks                                                 [default = Tss]\n")
    cat("--Height     : plot height                                                     [default = max of colMeans]\n")
    cat("--upStream   : number of bp upstream of Tss in matrix                          [required]\n")
    cat("--downStream : number of bp downstream of Tss in matrix                        [required]\n")
    cat("--Label      : Alternate plot label if not Type (example: Pause_site)          [default = Type]
                         Type is used to grep the files. If file name match use Type.                   \n")

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
    upStream   <- as.numeric( sub( '--upStream=', '',args[grep('--upStream=',args)]) )
    downStream <- as.numeric( sub( '--downStream=', '',args[grep('--downStream=',args)]) )
    Label      <- sub( '--Label=', '',args[grep('--Label=',args)])
}

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

## set defaults
if (identical(bins,character(0))){
   bins <- 25
}else{
    bins <- as.numeric(bins)
}

if (identical(Type,character(0))){
   Type <- "Tss"
}

if (identical(Label,character(0))){
   Label <- Type
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
SAMPLES <- SAMPLES[grep("^anti.df$|^sense.df$", SAMPLES, invert=TRUE)]
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
    Cols <- viridis(length(foo))#rainbow(length(foo))#
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE,  comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[sub(".df", "", SAMPLES), "color"])
}

Labels <- c(paste0("-",upStream)
           ,Label
           ,paste0("+", downStream))


if(bins==0){
    Breaks <- c(0
               ,upStream
               ,upStream + downStream)
}else{
    Breaks <- c(0
               ,upStream/bins
               ,(upStream + downStream)/bins)   
}

pdf(file=sub("$", paste0(Type, ".metaPlot.pdf"), outName),width=8,height=5)
print({
    p <-
    ggplot(df.long, aes(x=x, y=value, color=variable)) +
    geom_line(size=1)+
    scale_color_manual("sample", values = Cols )+
    scale_x_continuous(breaks = Breaks
                      ,labels = Labels
                       )+
    ylab("signal")+
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
