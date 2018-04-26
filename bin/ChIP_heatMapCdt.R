args <- commandArgs()

help <- function(){
    cat("ChIP_heatMapCdt.R :
- Heatmap cdt sorted by occupancy in the control.
- This script is meant for multiple files that need to be compared to a control.
- A fold change cdt is also made in the same order.\n")
    cat("Usage: \n")
    cat("--Dir        : path to heatmap tables                                                           [required]\n")
    cat("--Pattern    : grep pattern for heatmap tables to use quotes (ex. PolII.*293T)                  [optional]\n")    
    cat("--Control    : number of bp upstream of Tss in matrix                                           [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir        <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern    <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    Control    <- sub( '--Control=', '',args[grep('--Control=',args)])
}

## select rda files for pattern of interest
Dir
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
foo

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda|_Peaks.*rda|_Tss.*rda",".df",basename(foo[i]))
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



fc <- do.call(cbind, lapply(1:length(exps), function(i){
    print(paste("FC sample", i, exps[i]))
    exp         <- as.matrix(get(exps[i]))
    con         <- get(SAMPLES[SAMPLES==paste0(Control, ".df")])
    pseudo      <- min(c(min(exp[!exp == 0]), min(con[!con == 0])))
    exp         <- exp + pseudo
    con         <- con + pseudo
    FC          <- log2(exp/con)
    FC
}))

or            <- order(-rowSums(get(SAMPLES[SAMPLES==paste0(Control, ".df")])))

write.table(data.frame(HeatOrder=rownames(fc[or,]))
           ,file=paste(Dir, paste(sub("$", ".combined", Control), "ocpSortRowIds.txt", sep="."), sep="/")
           ,sep="\t",row.names=FALSE, quote=FALSE)


cdt      <- cbind( UID=or,NAME=or, df.all[or,] ) 
fc.cdt   <- cbind( UID=or,NAME=or, fc[or,]  )

## write tables
write.table(cdt,file=paste(Dir, paste0(Control, ".combined.ocpSort.cov.cdt"), sep="/"), sep="\t",row.names=FALSE)
write.table(fc.cdt, file=paste(Dir, paste0(Control, ".combined.ocpSort.log2FC.cdt"),  sep="/"), sep="\t",row.names=FALSE)
