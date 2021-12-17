files.sources = dir("helperScriptsFunctions/")
files.sources <- paste0("helperScriptsFunctions/",files.sources[endsWith(files.sources, ".R")])
invisible(sapply(files.sources, source))


## Run Random Search on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 1,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_RandomSearch")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run 1+1ES on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 4,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_ES")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run CMAES on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 8,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_CMAES")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run Multi-Local EI on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 13,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_MultiLocalEI")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run moiMBO on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 5,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_MOIMBO")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run qEI on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 6,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_QEI")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

## Run IPI on BBOB Funcs
paramList = list("seed" = 1:15, 
                 "functionID" = 1:24, 
                 "algoID" = 7,
                 "nDim" = c(2,5),
                 "budget" = 150,#50,
                 "batchSize" = 2,#c(2,4,8,16),
                 "expResults/exp1_IPI")
print(paste("This will produce", nrow(expand.grid(paramList)), "experiments!"))
runAlgorithmByID(paramList)

