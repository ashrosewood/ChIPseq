args <- commandArgs()

help <- function(){
    cat("ChIP_kmeansHeatMap.R :
- Combine heat map matricies (log2FC or average coverage) into one table for heatmaps.
- This table can then be passed to ChIP_kmeansHeatMap.R or ChIP_heatMapCdt.R depending if you want to cluster.\n")
    cat("Usage: \n")
    cat("--Dir         : path to heatmap tables                                           [required]\n")
    cat("--Pattern     : grep pattern for heatmap tables to use quotes (ex. PolII.*293T)  [default = all .rda files in a directoy]\n")    
    cat("--sampleOrder : Order to combine matrices (separated by a pipe |)                [default = default string order]
                             do not include the suffix in the files (.metaGene.rda|_Peaks.*rda|_Tss.*rda|_Tes.*rda)\n")  
    cat("--outName     : prefix to your out file names (No .extention)                    [required]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    Dir         <- sub( '--Dir=', '', args[grep('--Dir=', args)])
    Pattern     <- sub( '--Pattern=', '', args[grep('--Pattern=', args)])
    sampleOrder <- sub( '--sampleOrder=', '',args[grep('--sampleOrder=',args)])
    outName     <- sub( '--outName=', '',args[grep('--outName=',args)])
}

## show params
cat("Dir:", Dir, sep="\n" )
cat("Pattern:", Pattern, sep="\n" )
cat("outName:", outName, sep="\n" )

if (identical(sampleOrder,character(0))){
    reOrder <- 0
    cat("sampleOrder:", sampleOrder, sep="\n" )
}else{
    reOrder <- 1
    cat("sampleOrder:", "not provided", sep="\n" )
}

## make output directory if does not exist
if(!(file.exists( dirname(outName) ))) {
    print(paste("mkdir", dirname(outName)))
    dir.create(dirname(outName),FALSE,TRUE)  
}

## select rda files for pattern of interest
foo <- list.files(Dir, pattern=".rda", full.names=TRUE)
foo <- foo[grep(Pattern,foo,invert=FALSE)]
cat("Tables:", foo, sep="\n")

## load files and assign to file name
for (i in 1:length(foo))
{
    oname = gsub(".metaGene.rda|_Peaks.*rda|_Tss.*rda|_Tes.*rda",".df",basename(foo[i]))
    oname <- gsub("-","_",oname)
    df <- get(load(foo[i]))
    colnames(df) <- paste(sub(".df", "", oname), 1:ncol(df), sep=".")
    assign(oname, df)
}

SAMPLES <- ls(pattern=".df$")
SAMPLES

## sanity check must be the same genes in the same order
stopifnot(rownames(get(SAMPLES[1]))==rownames(get(SAMPLES[2])))

if(reOrder==1){
    ## check that the sample order is the same as the number of SAMPLES
    stopifnot( length(grep(sampleOrder, SAMPLES)) == length(SAMPLES) )
    newOrder     <- strsplit( sampleOrder, "\\|" )[[1]]
    cat( "matrix order:", newOrder, sep="\n" )
    regX        <-  function(x){
                        grep(x,SAMPLES)
    }
    lapply(newOrder, function(x){ regX(x) }) 

    SAMPLES.or   <- SAMPLES[do.call(rbind, lapply(newOrder, function(x){ regX(x) }) )[,1]]
    df           <- do.call(cbind, lapply(SAMPLES.or,
                                          function(x){
                                              get(x)
                                          }
                                          )
                            )

    save(df, file=sub("$", ".rda", outName))
}else{
    df           <- do.call(cbind, lapply(SAMPLES,
                                          function(x){
                                              get(x)
                                          }
                                          )
                            )
    save(df, file=sub("$", ".rda", outName))    
}

print("done")
