library(SPOT)
require(reticulate)
library(dplyr)

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
        matrix(apply(x,1,fun), ncol=1)
    }
    
    designSize <- length(lower) * 2 * batchSize
    X <- designLHD(,lower,upper,control = list(size = designSize))
    y <- tfun(X)
    
    optimizerControl = list(funEvals = as.integer(log(length(lower)) * as.integer(1000/batchSize)), 
                            populationSize = 5 * length(lower))
    
    modelControl = list(target = c("y", "s"),
                        useLambda = T, 
                        lambdaUpper = -4,
                        algTheta = optimDE,
                        budgetAlgTheta = 250)
    
    while(nrow(X) < maxEvals){
        print(paste("nrow(x):",nrow(X)))
        model <- buildKriging(X,y,modelControl)
        model$target <- c("y","s")
        
        evalIPI = function(x,tVal){
            pred <- predict(model, x)
            yVal = pred$y
            yMin = min(y)
            s = pred$s
            return(-(0.5*pnorm((yMin-yVal)/(1.05-tVal))
                     +pnorm((-(s-tVal)^2)/(0.05))))
        }
        
        tVals <- NULL
        for(i in 1:batchSize)
            tVals = c(tVals,i*(1/(batchSize+1)))
        
        newX = NULL
        for(tVal in tVals)
        {
            newX = rbind(newX, optimDE(,fun = function(x){evalIPI(x,tVal = tVal)}, lower,upper,optimizerControl)$xbest)
        }
        newY <- tfun(newX)
        
        X <- rbind(X, newX)
        y <- rbind(y,newY)
    }
}

wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("07IPI",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)
