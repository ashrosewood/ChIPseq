args <- commandArgs()

help <- function(){
    cat("ChIP_plotEcdf.R :
- Plot ecdf of pausing indexes calculated from ChIP_CalcPausingIndex.R and save as pdf\n")
    cat("Usage: \n")
    cat("--table   : table with tss and body averages (pausingIndexAverageCoverages.txt)   [required]\n")    
    cat("--outName : defaults to table prefix in current dir (do not put extention)        [default = basename(bigWigFile) ]\n")
    cat("--cols    : need the same number as samples separated by comma                    [default = viridis palette]\n")
    cat("--Xmax    : xlim max for ecdf plot                                                [default = max log2(PI)]\n")
    cat("--Xmin    : xlim min for ecdf plot                                                [default = min log2(PI or PRR)]\n")
    cat("--PRR     : plot the promoter release ratio (body/promoter) instead of pi (0/1)   [default = 0 plot pausing index]\n")
    cat("--Control : Sample control prefix (no _tss or _body) name to perform ks test with [default = no p-values calculated]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    tab     <- sub( '--table=', '', args[grep('--table=', args)] )
    outName <- sub( '--outName=', '',args[grep('--outName=',args)])
    cols    <- sub( '--cols=', '',args[grep('--cols=',args)])
    Xmax    <- sub( '--Xmax=', '',args[grep('--Xmax=',args)])
    Xmin    <- sub( '--Xmin=', '',args[grep('--Xmin=',args)])
    PRR     <- sub( '--PRR=', '',args[grep('--PRR=',args)])
    Control <- sub( '--Control=', '',args[grep('--Control=',args)])
}

## debug
#setwd("/projects/b1025/arw/analysis/kevin/SEC/")
#tab <- "tables/pausing_index/Pol2_293T_DMSO_817_rep1.pausingIndexAverageCoverages.txt"
#outName <- "plots/ecdf/Pol2_293T_DMSO_817_rep1.pausingIndex"
#Xmax <- 10
#cols <- "plots/plotColors/293T_rep1_817_colors.txt"
#PRR <- 0
#Control <- "PolII_DMSO_293T_817_rep1"

if (identical(outName,character(0))){
   outName <- gsub("AverageCoverages|.txt", "", basename(df))
}

if (identical(PRR,character(0))){
    PRR <- 0
}else{
    PRR <- as.numeric(PRR)
}

if (identical(Control,character(0))){
    pVals <- 0
}else{
    pVals <- 1
}

library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)

###############################
## calculate pausing index
###############################

df           <- read.delim(tab)
rownames(df) <- paste(df$gene_id)
dim(df)

Tss          <- df[ , grep("tss$", names(df)) ]
Body         <- df[ , grep("body$", names(df)) ]

## be sure they are in the same order
Order        <- sub( "_tss", "", names(Tss) )
Tss.or       <- Tss[ , paste(Order,"tss", sep="_") ]
Body.or      <- Body[ , paste(Order,"body", sep="_") ]

stopifnot( sub("_tss", "", colnames(Tss.or)) == sub("_body", "", colnames(Body.or)))
stopifnot( rownames(Tss.or) == rownames(Body.or) )

## pausing index
PI <- do.call(cbind, lapply(Order, function(x){
    pi           <- log2( Tss.or[, paste0(x, "_tss")] / Body.or[, paste0(x, "_body")] )
    pi           <- as.data.frame(pi)
    rownames(pi) <- rownames(Tss.or)
    colnames(pi) <- sub("$", "_pi", x)
    pi
  }
))

Prr <- do.call(cbind, lapply(Order, function(x){
    pi           <- log2( (Body.or[, paste0(x, "_body")]) / (Tss.or[, paste0(x, "_tss")]) )
    pi           <- as.data.frame(pi)
    rownames(pi) <- rownames(Tss.or)
    colnames(pi) <- sub("$", "", x)
    pi
  }
))

###############################
## make ecdf plot
###############################
if (identical(cols,character(0))){
   Cols <- viridis(length(Order))
}else{
    df.col           <- read.table(cols,sep="\t", header=TRUE, comment.char = "")
    rownames(df.col) <- paste(df.col$sample)
    Cols             <- paste(df.col[Order, "color"])
}

if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

if(PRR==0){
    df.long <- melt(PI[ is.finite(rowSums(PI)), ] )
    if (identical(Xmax,character(0))){
        Xmax <- round(max(df.long$value))
    }else{
        Xmax <- as.numeric(Xmax)
    }
    if (identical(Xmin,character(0))){
        Xmin <- round(min(df.long$value))
    }else{
        Xmin <- as.numeric(Xmin)
    }
    pdf(file=sub("$", ".pdf", outName),width=8,height=5)
    print({
        p <-ggplot(df.long, aes(x=value, colour=variable)) +
            stat_ecdf(size=1.5) +
            scale_color_manual("sample", values = Cols )+
            ylab("fraction")+
            xlab("log2(pausing index)")+
            xlim(Xmin,Xmax)+
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
    if(pVals > 0){
        exps <- names(PI)[grep(Control, names(PI), invert=TRUE)]
        pvals.df <- do.call(cbind, lapply(1:length(exps), function(i){
            print(paste("ks.test Pval", exps[i], "vs Control:", Control))
            p <- ks.test(PI[, exps[i] ]
                        ,PI[, paste0(Control, "_pi") ]
                        ,alternative="two.sided")
            p.df <- data.frame(t(data.frame(p$p.value, p$statistic)))
            names(p.df) <- paste0(exps[i], "_pval")
            p.df
        }))
        write.table(pvals.df, sub("$", ".PI.pvalues.txt", outName), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) 
    }
    PI$gene_id <- rownames(PI)
    write.table(PI, sub("$", ".PI.txt", outName), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) 
}else{
    df.long <- melt(Prr[ is.finite(rowSums(Prr)), ] )
    if (identical(Xmax,character(0))){
        Xmax <- round(max(df.long$value))
    }else{
        Xmax <- as.numeric(Xmax)
    }
    if (identical(Xmin,character(0))){
        Xmin <- round(min(df.long$value))
    }else{
        Xmin <- as.numeric(Xmin)
    }
    pdf(file=sub("$", ".pdf", outName),width=8,height=5)
    print({
        p <-
            ggplot(df.long, aes(x=value, colour=variable)) +
            stat_ecdf(size=1.5) +
            scale_color_manual("sample", values = Cols )+
            ylab("fraction")+
            xlab("log2(PRR)")+
            xlim(Xmin,Xmax)+
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
    if(pVals > 0){
        exps <- names(Prr)[grep(Control, names(Prr), invert=TRUE)]
        pvals.df <- do.call(cbind, lapply(1:length(exps), function(i){
            print(paste("ks.test Pval", exps[i], "vs Control:", Control))
            p <- ks.test(Prr[, exps[i] ]
                        ,Prr[, Control ]
                        ,alternative="two.sided")
            p.df <- data.frame(t(data.frame(p$p.value, p$statistic)))
            names(p.df) <- paste0(exps[i], "_pval")
            p.df
        }))
        write.table(pvals.df, sub("$", ".PRR.pvalues.txt", outName), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) 
    }
    dev.off()
    Prr$gene_id <- rownames(Prr)
    write.table(Prr, sub("$", ".PRR.txt", outName), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) 
}
