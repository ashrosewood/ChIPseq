args <- commandArgs()

help <- function(){
    cat("ChIPseq_differentialPeaksEdgeR.R :
- NOTE: this program is for when you have 2 or more replicates for the condition and control.
- You MUST also provide an option of 0 to CalcNormFactors in order to turn this off.
- Each row is the name of the comparison output and the column names are the sample with replicate names.
- In the comparison file 1 means control, -1 means experiment, 0 means do not consider for this comparison\n")
    cat("Usage: \n")
    cat("--countFile       : counts file (first 6 columns: seqnames, start, end, width, strand, name) [required]\n")
    cat("--comparisonFile  : file that defines what the control and experiment are in each comparison [required]\n")
    cat("--outputDirectory : Path to write output files                                               [required]\n")
    cat("--runMDS          : Generate mds plot for the top500 regions (0/1)                           [default = 1]\n")
    cat("--libSizeFile     : Library sizes to be used rather than column sums for normalization       [default = uses column sums]
                            Must have the same column names as the samles in the countFile.\n")
    cat("--CalcNormFactors : Normalize based on library size (0/1)                                    [default = 1]
                            Turn off if you already have normalzied counts (ex. spike in normalized).\n")
    cat("--adjp            : FDR adjusted p-value cutoff for significant flag columns                 [default = 0.05]\n")
    cat("--FC              : Foldchange cutoff for significant flag columns                           [default = 1.5]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    countFile       <- sub('--countFile=', '', args[grep('--countFile=', args)])
    comparisonFile  <- sub('--comparisonFile=', '', args[grep('--comparisonFile=',args)])
    outputDirectory <- sub('--outputDirectory=', '',args[grep('--outputDirectory=',args)])
    runMDS          <- sub('--runMDS=', '',args[grep('--runMDS=',args)])    
    CalcNormFactors <- sub('--CalcNormFactors=', '', args[grep('--CalcNormFactors=', args)])
    libSizeFile     <- sub('--libSizeFile=', '', args[grep('--libSizeFile=',args)])
    adjp            <- sub('--adjp=', '', args[grep('--adjp=',args)])
    FC              <- sub('--FC=', '', args[grep('--FC=',args)])
}

## libsize
if (identical(libSizeFile,character(0))){
    normLibSizes <- 0
}else{
    normLibSizes <- 1
}

print(c("normLibSize:", normLibSizes))

## make outputDirectory if needed
if(!(file.exists(outputDirectory))) {
      dir.create(outputDirectory,FALSE,TRUE)
}

# If runMDS isn't specified, default is to run it.
if (identical(runMDS,character(0))){
   runMDS <- 1
}

# If CalcNormFactors isn't specified, default is to normalize counts by library size.
if (identical(CalcNormFactors,character(0))){
    CalcNormFactors<-1
}else{
    CalcNormFactors <- as.numeric(CalcNormFactors)
}
print(c("CalcNormFactors:", CalcNormFactors))
      
# If comparison file isn't specified, default is to skip it.
runComp<-1
if (identical(comparisonFile,character(0))){
   runComp<-0
}

# If output directory isn't specified, default is current directory
if (identical(outputDirectory,character(0))){
   outputDirectory<-""
} else { 
  outputDirectory<-paste(outputDirectory,"/",sep="")
}

if (identical(adjp,character(0))){
   adjp <- 0.05
}else{
   adjp <- as.numeric(adjp)
}

if (identical(FC,character(0))){
   FC <- 1.5
}else{
   FC <- as.numeric(FC)
}


library(edgeR)
library(GenomicRanges)
library(biomaRt)

##----------load counts --------#
##load count data
print("Loading count data")
print(countFile)

if(grepl('rda',countFile)){
   counts <- get(load(file=countFile))
}
if(grepl('txt',countFile)){
    counts <- read.delim(file=countFile,header=TRUE,sep="\t")
}
rownames(counts) <- counts$peak
head(counts)
dim(counts)

print("Peek at gnModel")
gnModel <- counts[,1:6]
head(gnModel)

##----------MDS plot--------------------#

allData <-data.frame(counts[,7:length(colnames(counts))])	

head(allData)

print (runMDS)
if (runMDS==1){
    print("Generating MDS plot.")
    newnames <- counts[,7:length(colnames(counts))]
    ##newnames <- gsub("\\.\\d+$","",newnames)
    dge <- DGEList(allData)      #filter low counts & calc lib.sizes
    if(normLibSizes>0){
        libSize              <- read.table(file=libSizeFile,sep="\t",as.is=T,header=T)  
        dge$samples$lib.size <- as.numeric(libSize[,names(allData)])
    }
    if(CalcNormFactors > 0){
        dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (composition effect)
    }
    dge <- estimateCommonDisp(dge,verbose=T)
    dge <- estimateTagwiseDisp(dge)
    p <- plotMDS(dge, method="logFC", top=500)
    df <- data.frame(x=p$x, y=p$y, name=names(p$x))
    pdfname <-paste(outputDirectory,"allData.MDS.pdf",sep="")
    pdf(pdfname)
    plot(df$x, df$y, pch=19, col="blue", main="Top 500 Peak Counts MDS"
        ,cex=1.2
        ,cex.lab=1.2
        ,cex.main=1.2
        ,cex.axis=1.2
        ,ylim=c(min(df$y)-0.25,max(df$y)+0.25)
        ,xlim=c(min(df$x)-0.5,max(df$x)+0.5)
        ,ylab="Leading logFC dim 2"
        ,xlab="Leading logFC dim 1"
         )
    text(df$x, df$y, labels=df$name, cex=0.9, pos=3)
    dev.off()
}

##-----------gene expression--------------#
runEdgeR <- function(data,comparison){
    head(data)
    print(comparison)
    group1 <- c()
    group2 <- c()
    ## for this comparison extract controls and experiments counts 
    for(i in 1:length(comparison)){
        value <- as.integer(comparison[i])
        if (value == 1){
            group1<-c(group1,(colnames(data)[i]))
        }						
        if (value == -1){				
            group2<-c(group2,(colnames(data)[i]))
        }	  		  
    }
    group1  <- data.frame(subset(data,select=group1))
    group2  <- data.frame(subset(data,select=group2))
    subdata <-cbind(group1,group2)
    setup   <- c(rep("group1",each=length(group1)),rep("group2",each=length(group2)))
    ## print groups
    print(head(group1))
    print(head(group2))
    print(head(subdata))
    print(setup)
    ## make design for edgeR
    grp              <- factor(setup)
    grp              <- relevel(grp,ref="group1")
    design           <- model.matrix(~grp)
    colnames(design) <-levels(grp)
    print(design)
    ##
    ## run edgeR with glm's 
    ##
    dge <- DGEList(subdata, group=grp)#DGEList(subdata[rowSums(cpm(subdata)>1)>=2,], group=grp) #filter low counts & calc lib.sizes
    if(normLibSizes>0){
        dge$samples$lib.size <- as.numeric(libSize[,names(subdata)])
    }
    if(CalcNormFactors > 0){
        dge <- calcNormFactors(dge)  #calcs normalization factors based on lib.size (RNA composition effect)
    }
    dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
    dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
    dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter 
    fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
    lrt <- glmLRT(fit,coef="group2")
    ## annotate edgeR table 
    foo <- rownames(lrt$table)                #List of differentially expressed genes
    stopifnot(rownames(lrt$table)==rownames(gnModel[foo,]))
    stopifnot(rownames(dge$counts)==rownames(gnModel[foo,]))
    ## calculate RPKM using exonLength
    stopifnot(rownames(lrt$table)==rownames(dge$counts))
    ## combine gnModel, RPKM, edgeR table, and calculate the adjusted pval
    df            <- cbind(as.data.frame( gnModel[foo,]),dge$counts
                          ,signif(lrt$table, 4)
                          ,adj.p=signif(p.adjust(lrt$table$PValue,method="BH"), 4)
                           )
    names(df)     <- sub("log", "log2", names(df))
    totals        <- data.frame(t(lrt$samples[,"lib.size"]))
    names(totals) <- sub("$", "_total", rownames(lrt$samples))
    df            <- merge(df, totals) ## add total reads to the df
    ## add flags
    cat("adding flags\n")    
    up <- df$adj.p < adjp & df$log2FC > log2(FC)
    flag <- rep(0, nrow(df))
    flag[up] <- 1
    df$up <- flag
    dn <- df$adj.p < adjp & df$log2FC < -log2(FC)
    flag1 <- rep(0, nrow(df))    
    flag1[dn] <- 1
    df$dn <- flag1
    df
}

if (runComp == 1) {
   print(comparisonFile)
   comparisons <- read.table(file=comparisonFile,sep="\t",as.is=T,header=T,row.names=1)
   if (normLibSizes>0){
       libSize <- read.table(file=libSizeFile,sep="\t",as.is=T,header=T)    
   }
   print("Comparisons file:")
   print(colnames(comparisons))
                                        #colnames(comparisons) <- gsub("\\_","\\.",colnames(comparisons))
   allDataSampleOrder <- colnames(allData)
   print("SampleOrder:")
   print(allDataSampleOrder)
   comparisons <- comparisons[,allDataSampleOrder]
   print(allDataSampleOrder)
   print(colnames(comparisons))
   for(c in 1:length(rownames(comparisons))){
       comparison <- (paste(c(comparisons[c,1:length(colnames(comparisons))])))
       print("Running EdgeR comparison")
       print(paste(c(comparisons[c,1:length(colnames(comparisons))])))
       print(rownames(comparisons)[c])
       comp <-runEdgeR(data=allData,comparison)
       label1<- paste(outputDirectory,paste(rownames(comparisons)[c],"edgeR.df.rda",sep="."),sep="")
       label2<- paste(outputDirectory,paste(rownames(comparisons)[c],"edgeR.txt",sep="."),sep="")
       save(comp,file=label1)
       write.table(comp,file=label2,sep="\t",col.names=TRUE, quote=F, row.names=FALSE)
   }
}

