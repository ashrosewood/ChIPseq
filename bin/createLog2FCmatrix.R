args <- commandArgs()

help <- function(){
    cat("createLog2FCmatrix.R :
- Create log2 fold change (exp/con) data frame for an experiment and control data frame in RData or rda format of equal dimensions.
- A pseudo count will be added using the minimal non zero value in both datasets.\n")
    cat("Usage: \n")
    cat("--conFile  : Control matrix file path (FC denominator)                [required]\n")
    cat("--expFile  : Experiment matrix file path (FC numerator)               [required]\n")    
    cat("--outName  : Prefix for your output file (extension excluded)         [default = expFile.log2FC.rda]\n")    
    cat("\n")
    q()
}

## Save values of each argument

if(length(args)==0 || !is.na(charmatch("-help",args))){
    help()
} else {
    conFile     <- sub('--conFile=', '', args[grep('--conFile=', args)])
    expFile     <- sub('--expFile=', '', args[grep('--expFile=', args)])
    outName     <- sub('--outName=', '', args[grep('--outName=', args)])
}

print(paste("conFile:", conFile))
print(paste("expFile:", expFile))

if (identical(outName,character(0))){
   outName <- sub(".rda$|.RData$", ".log2FC.rda", expFile)
}

con        <- get(load(conFile))
exp        <- get(load(expFile))

pseudo     <- min( c(min(exp[exp!=0]) , min(con[con!=0])) )
print(paste("pseudo used:", signif(pseudo)))

if(dim(con)==dim(exp)){
    print("data frames are the same dimention")    
    log2FC <- log2( (exp + pseudo) / (con + pseudo) )
    save(log2FC, file=outName)
}else{
    print("data frames are NOT the same dimention\n log2FC dataframes will not be made.")
}

print("done")
