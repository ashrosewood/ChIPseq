args <- commandArgs()

help <- function(){
    cat("ChIP_makeMetaGeneHeatMap.R :
- Use kmeans clustering to cluster metagene heatmaps cdts.
- This script is meant for multiple files that need to be compared to a control.
- They will be separated by cluster than sorted by rowSums in each cluster
- A fold change cdt is also made.\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                                           [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T)                  [optional]\n")    
    cat("--clusNum    : number of clusters for heatmap                                                   [default = 4]\n")
    cat("--Control    : number of bp upstream of Tss in matrix                                           [required]\n")
    cat("--Combine    : make a cdt file with all samples the control in first then the rest sorted (0/1) [default = 1]\n")
    cat("--Separate   : make the cdt files for each sample separate (0/1)                                [default = 0]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir        <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern    <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    clusNum    <- sub( '--clusNum=', '',args[grep('--clusNum=',args)])
    Control    <- sub( '--Control=', '',args[grep('--Control=',args)])
    Combine    <- sub( '--Combine=', '',args[grep('--Combine=',args)])
    Separate   <- sub( '--Separate=', '',args[grep('--Separate=',args)])
}

if (identical(Combine,character(0))){
   Combine <- 1
}else{
   Combine <- as.numeric(Combine)
}

if (identical(Separate,character(0))){
   Separate <- 0
}else{
   Separate <- as.numeric(Separate)
}

if (identical(clusNum,character(0))){
   clusNum <- 4
}else{
   clusNum <- as.numeric(clusNum)
}

#Dir <- "tables/metaGene"
#Pattern <- "H3K4me3_111816"
#Control <- "cps25_WT_LOG_H3K4me3_111816"
#FC <- 1
#clusNum <- 3
#outName <- "tables/metaGene/H3K4me3_111816.metaGene.cdt"
#setwd("../../yohhei/cps25/")

## select rda files for pattern of interest
Dir
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
    colnames(df) <- paste(sub(".df", "", oname), 1:ncol(df), sep=".")
    assign(oname, df)
}

SAMPLES <- ls(pattern=".df$")
SAMPLES

## sanity check must be the same genes in the same order
stopifnot(rownames(get(SAMPLES[1]))==rownames(get(SAMPLES[2])))

## put the control first and set experiments
exps <- SAMPLES[SAMPLES!=paste0(Control, ".df")]

## report the order in the data frame
df.all <- get(SAMPLES[SAMPLES==paste0(Control, ".df")])
print(paste("control sample 1", paste0(Control, ".df")))
for( i in 1:length(exps) ){
    print(paste("sample", i+1, exps[i]))
    df.all <- cbind(df.all, get(exps[i]))
}


## coverage heatmap
print( paste("cluster by kmeans for", clusNum, "clsuters") )
## kmeans parameters
go.paras       <- list(knc=clusNum, max.iter=20, nrs=30)
## force kmeans to use the same random seed everytime
set.seed(20)    
km             <- kmeans(df.all, centers=go.paras$knc, 
                         iter.max=go.paras$max.iter, nstart=go.paras$nrs)
## get cluster numbers and sort by occupancy within each cluster
clus          <- data.frame(km$cluster)
clus$totalCov <- rowSums(df.all)
clus.or <- clus[with(clus, order(km.cluster, -totalCov)), ]

## report the number in each cluster
for ( x in 1:clusNum){
    print(paste("cluster", x, sum(clus.or$km.cluster==x), "genes")) 
}    

## make the cdt files for each sample
if( Separate==1 ){
    for ( i in 1:length(SAMPLES) ){
        sub.df <- get(SAMPLES[i])
        ## add arbitrary divider
        sep.cdt <- data.frame()
        for ( x in 1:clusNum){
            print(paste("Sample",SAMPLES[i],"cluster",  x))
            cdt <- cbind( UID=rownames(clus.or[clus.or$km.cluster==x,])
                        ,NAME=rownames(clus.or[clus.or$km.cluster==x,])
                        ,sub.df[rownames(clus.or[clus.or$km.cluster==x,]),]
                         )
            fake <- data.frame(cbind(UID=paste0("clus",x)
                                    ,NAME=paste0("clus",x)
                                    ,t(data.frame(num=rep(ceiling(max(df.all)),ncol(sub.df))))
                                     ))
            colnames(fake) <- colnames(cdt)
            rownames(fake) <- paste0("clus",x)
            sep.cdt <- rbind(sep.cdt, cdt, fake)
        }
        ## write the cdt for each sample
        write.table(sep.cdt
                   ,file=paste(Dir, paste(sub(".df", "", SAMPLES[i]), "kmeans", clusNum, "cdt", sep="."), sep="/")
                   ,sep="\t",row.names=FALSE)
    }
    for (i in 1:length(exps)){
        con.df                  <- get(SAMPLES[SAMPLES==paste0(Control, ".df")])
        sub.df                  <- get(exps[i])
        fc.cdt                  <-  data.frame()
        for (x in 1:clusNum){
            print(paste("Sample FC",exps[i],"cluster",  x))
            or                  <- rownames(clus.or[clus.or$km.cluster==x,])
            FC                  <- as.matrix(sub.df[or,] / con.df[or,])
            FC[is.nan(FC)]      <- 1
            FC[is.infinite(FC)] <- max(FC[!FC %in% "Inf"])
            FC                  <- log2(FC)
            FC[is.infinite(FC)] <- min(FC[!FC %in% "-Inf"])
            cdt                 <- cbind( UID=rownames(clus.or[clus.or$km.cluster==x,])
                                        ,NAME=rownames(clus.or[clus.or$km.cluster==x,])
                                        ,FC
                                         )
            fake                <- data.frame(cbind(UID=paste0("clus",x)
                                                   ,NAME=paste0("clus",x)
                                                   ,t(data.frame(x=rep(0,ncol(FC))))
                                                    ))
            colnames(fake)      <- colnames(cdt)
            rownames(fake)      <- paste0("clus",x)
            fc.cdt              <- rbind(fc.cdt, cdt, fake)    
        }
        write.table(fc.cdt
                   ,file=paste(Dir, paste(sub(".df", "", SAMPLES[i]), "kmeans", clusNum, "log2FC.cdt", sep="."), sep="/")
                   ,sep="\t",row.names=FALSE)
    }
}


## make combined heatmap with control first
if( Combine==1){
    all.cdt <- data.frame()
    ## add arbitrary divider for the combined heatmap
    for ( x in 1:clusNum){
        print(paste("Combined heat cluster",  x))
        cdt <- cbind(UID=rownames(clus.or[clus.or$km.cluster==x,])
                    ,NAME=rownames(clus.or[clus.or$km.cluster==x,])
                    ,df.all[rownames(clus.or[clus.or$km.cluster==x,]),]
                     )
        fake <- data.frame(cbind(UID=paste0("clus",x)
                                ,NAME=paste0("clus",x)
                                ,t(data.frame(num=rep(ceiling(max(df.all)),ncol(df.all))))
                                 ))
        colnames(fake) <- colnames(cdt)
        rownames(fake) <- paste0("clus",x)
        all.cdt <- rbind(all.cdt, cdt, fake)   
    }
    ## write the cdt for each sample
    write.table(all.cdt
               ,file=paste(Dir, paste(sub("$", ".combined", Control), "kmeans", clusNum, "cdt", sep="."), sep="/")
               ,sep="\t",row.names=FALSE)
    ## write combined FC in the order of exps
    fc.cdt                  <- data.frame()
    for( x in 1:clusNum){
        print(paste("Combined FC heat cluster",  x))
        or                  <- rownames(clus.or[clus.or$km.cluster==x,])
        fc.mat              <- do.call(cbind, lapply(1:length(exps), function(j){
            con.df              <- get(SAMPLES[SAMPLES==paste0(Control, ".df")])
            sub.df              <- get(exps[j])
            FC                  <- as.matrix(sub.df[or,] / con.df[or,])
            FC[is.nan(FC)]      <- 1
            FC[is.infinite(FC)] <- max(FC[!FC %in% "Inf"])
            FC                  <- log2(FC)
            FC[is.infinite(FC)] <- min(FC[!FC %in% "-Inf"])
            FC
        }))
        ## 
        cdt <- cbind(UID=or
                    ,NAME=or
                    ,fc.mat
                     )
        fake <- data.frame(cbind(UID=paste0("clus",x)
                                ,NAME=paste0("clus",x)
                                ,t(data.frame(num=rep(0,ncol(fc.mat))))
                                 ))
        colnames(fake) <- colnames(cdt)
        rownames(fake) <- paste0("clus",x)
        fc.cdt <- rbind(fc.cdt, cdt, fake)
    }
        write.table(fc.cdt
                   ,file=paste(Dir, paste(sub("$", ".combined", Control), "kmeans", clusNum, "log2FC.cdt", sep="."), sep="/")
                   ,sep="\t",row.names=FALSE)    
}
