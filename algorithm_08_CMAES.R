require(reticulate)

files.sources = dir("helperScriptsFunctions/")
files.sources <- paste0("helperScriptsFunctions/",files.sources[endsWith(files.sources, ".R")])
invisible(sapply(files.sources, source))

args = commandArgs(trailingOnly=TRUE)

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
# 1) SEED
# 2) Function ID (1 = OBK, 2 = KOELN)
# 3) Algo ID

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

################# 
solver <- function(fun,lower,upper){
    popSize <- batchSize * 4
    
    cma<-import("cma")
    initial <- runif(length(lower),lower,upper)
    
    properties = list('bounds'= list(lower, upper), 'popsize'= popSize) 
    maxiter = maxIters/4 #-1
    es<-cma$CMAEvolutionStrategy(initial, 0.5, properties) 
    while(es$countiter<=maxiter){
        x<-es$ask()
        xp<-r_to_py(x)
        es$tell(xp,sapply(x,fun))
        gc()
    }
    return(es$result$fbest)
}

wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("08CMAES",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)
