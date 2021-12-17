print("libPaths:")
print(.libPaths())

files.sources = dir("helperScriptsFunctions/")
files.sources <- paste0("helperScriptsFunctions/",files.sources[endsWith(files.sources, ".R")])
invisible(sapply(files.sources, source))

library(SPOT)
require(reticulate)

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

## Definition of a simple 1+1-ES
do1plus1EsIter <- function(histList, f){
    sr=0.2 #succes rate limit
    a=1.2 #step size multiplier
    g=10 #length of memory list (for success rate)
    
    if(is.null(histList$sigma)){
        sigma=0.1 #initial step size
    }else{
        sigma=histList$sigma
    }
    
    x <- histList$x
    y <- histList$y
    memory <- histList$memory
    
    n <- length(x) #number of variables
    
    
    x1 <- x+sigma*rnorm(n,0,1)
    y1 <- f(x1)
    
    if(y1 < y){ #success
        x <- x1
        y <- y1
        if(length(memory)<g){
            memory <- c(memory,1)
        }else{
            memory <- c(memory[-1],1)
        }
    }else{ #fail
        if(length(memory)<g){
            memory <- c(memory,0)
        }else{
            memory <- c(memory[-1],0)
        }
    }
    
    if(length(memory)>=g){
        successrate <- sum(memory)/g
        print(paste("Sigma adaptation:",successrate))
        if(successrate>sr)
            sigma <- sigma * a
        else
            sigma <- sigma / a
        print(paste("New Sigma:",sigma))
    }	
    histList <- list(x=x,y=y,sigma=sigma,memory=memory)
    return(histList)
}

################# 1+1es
solver <- function(fun,lower,upper,solverParameterList){
    print(lower)
    
    for(j in 1:batchSize){
        print("new ES started!")
        xStart <- runif(length(lower), min = lower, max = upper)
        yStart <- fun(xStart)
        
        histList <- list(x=xStart,y=yStart)
        for(i in 1:(maxIters-1)){
            histList <- do1plus1EsIter(histList, fun)
        }
    }
}

wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("04Local",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)