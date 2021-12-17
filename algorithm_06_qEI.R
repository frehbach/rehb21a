print("libPaths:")
print(.libPaths())

library(SPOT)
require(reticulate)
library(dplyr)

require(DiceKriging)
require(DiceOptim)
require(mnormt)

files.sources = dir("helperScriptsFunctions/")
files.sources <- paste0("helperScriptsFunctions/",files.sources[endsWith(files.sources, ".R")])
invisible(sapply(files.sources, source))

args = commandArgs(trailingOnly=TRUE)

print("Args:")
print(args)

retries = 10
while(tryCatch({
    use_python("/usr/bin/python2", required = T)
    py_config()
    return(0)
}, error = function(e) {
    return(1)
})){
    retries <- retries - 1
    if(retries <= 0){
        break
    }
}

### RUN Parameters #########################################

### Recieve Setup
### 
seed <- as.numeric(args[1])
set.seed(seed)
funID <- as.numeric(args[2])
algoID <- as.numeric(args[3])

nDim <- as.numeric(args[4])
budget <- as.numeric(args[5])
batchSize <- as.numeric(args[6])
experimentPath <- args[7]
args <- args[-7]

maxIters <- budget
maxEvals <- maxIters * batchSize

############################
############################

solver <- function(fun,lower,upper,solverParameterList){
    ########target function wrapper for SPOT
    tfun <- function(x){
        as.numeric(apply(x,1,fun))
    }

    ## Create a design with nDim*maximum batchSize
    designSize <- length(lower) * 2 * batchSize
    X <- designLHD(,lower,upper,control = list(size = designSize))
    y <- tfun(X)
    
    optimizerControl = list(funEvals = as.integer(log(length(lower)) * 1000), 
                            populationSize = 5 * length(lower))
    
    while(nrow(X) < maxEvals){
        print(paste("nrow(x):",nrow(X)))
        model <- km(~1, design=X, response=y,
                      covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5),
                      nugget = 0.000001,nugget.estim = T, iso = F)
        getQEI <- function(x){
            res <- -qEI(matrix(x,nrow = batchSize),model)
            if(is.nan(res)) res <- 0
            return(res)
        }
        result <- optimDE(,fun = getQEI, rep(lower,batchSize),rep(upper,batchSize),optimizerControl)$xbest
        
        newX <- matrix(result,nrow = batchSize)
        newY <- tfun(newX)
        
        X <- rbind(X, newX)
        y <- c(y,newY)
    }
}


wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("06QEI",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)

