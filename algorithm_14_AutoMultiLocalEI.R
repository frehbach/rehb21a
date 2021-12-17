#### Algorithm depends on function curve!!!!
#### 
curveDf <- data.frame("batchSize" = c(16,8,4,2,1), "s" = c(1,1.7,3,5.3,7.8))

curveDf$budgetTakenPerIter <- rev(curveDf$s)

experimentPath <- "exAuto2804"

getEvalCost <- function(currentBatchSize){
    return(curveDf$budgetTakenPerIter[curveDf$batchSize == currentBatchSize])
}

print("libPaths:")
print(.libPaths())

library(SPOT)
source("objFun.R")
require(reticulate)
library(dplyr)

require(DiceKriging)
require(DiceOptim)
require(mnormt)
library(flacco)
library(checkmate)

args = commandArgs(trailingOnly=TRUE)
print("Args:")
print(args)

retries = 10
while(tryCatch({
    use_python("/usr/bin/python", required = T)
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
### Recieve Setup ###
### 
seed <- as.numeric(args[1])
set.seed(seed)
funID <- as.numeric(args[2])
algoID <- as.numeric(args[3])
nDim <- as.numeric(args[4])
budget <- as.numeric(args[5])
maxBatchSize <- as.numeric(args[6])

## In this algo the batchsize is the max possible batch Size!!

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
    histList <- list(x=x,y=y,sigma=sigma,memory=memory, x1 = x1, y1=y1)
    return(histList)
}

featureModel <- readRDS("autoRF.rds")

predictCores <- function(features, maxBatchSize){
    features <- features[featureModel$trainCols]
    features[is.na(features)] <- 0
    minVal <- -10e12
    maxVal <- 10e12
    features <- data.frame(t(sapply(features, function(y) min(max(y,minVal),maxVal))))
    p <- predict(featureModel$rf, as.matrix(features))
    p <- as.numeric(as.character(p))
    return(min(p, maxBatchSize))
}

optim1plus1ES <- function(x,f,control=list(),...){
    con<-list(
        sr=0.2, #succes rate limit
        sigma0=1, #initial step size
        a=1.2, #step size multiplier
        g=8,#length of memory list (for success rate)
        budget=100) #budget of function evaluations
    con[names(control)] <- control
    control<-con
    
    sigma <- control$sigma0
    a <- control$a
    budget <- control$budget
    g <- control$g
    sr <- control$sr
    
    n <- length(x) #number of variables
    
    y <- f(x)
    xhist <- x
    yhist <- y
    
    memory <- NULL
    
    for(i in 1:(budget-1)){
        x1 <- x+sigma*rnorm(n,0,1)
        y1 <- f(x1)
        xhist <- append(xhist,x1)
        yhist <- c(yhist,y1)
        
        if(y1 < y){ #success
            x <- x1
            y <- y1
            if(length(memory)<g)
                memory <- c(memory,1)	
            else
                memory <- c(memory[-1],1)	
        }else{ #fail
            if(length(memory)<g)
                memory <- c(memory,0)	
            else
                memory <- c(memory[-1],0)	
        }
        
        if(i>g){
            successrate <- sum(memory)/g
            if(successrate>sr)
                sigma <- sigma * a
            else
                sigma <- sigma / a			
        }		
    }
    return(list(x=x,y=y,sigma=sigma,xhist=xhist,yhist=yhist))
}

getPVPointDACE <- function(x, y, lower, upper){
    optimizerControl = list(funEvals = as.integer(log(length(lower)) * 1000), 
                            populationSize = 5 * length(lower))

    model <- km(~1, design=x, response=y,
                covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5),
                nugget = 0.000001,nugget.estim = T, iso = F)
    
    getPV <- function(xNew){
        krig <- predict(object = model, newdata = xNew, type = "UK", 
                        se.compute = FALSE, cov.compute = TRUE, checkNames = FALSE)
        krig$mean
    }
    
    optimDE(,fun = getPV, lower, upper,optimizerControl)$xbest
}

createEIPoints <- function(xAll, yAll, lower, upper, pointsEI){
    optimizerControl = list(funEvals = as.integer(log(length(lower)) * 1000 * pointsEI), ## 1000!!!!
                            populationSize = 5 * length(lower))
    
    if(pointsEI > 1){
        ## EI should be used to create multiple points in this iteration
        print(paste("Creating", pointsEI, "Points with q-EI"))
        #browser()
        model <- km(~1, design=xAll, response=yAll,
                    covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5),
                    nugget = 0.000001,nugget.estim = T, iso = F)
        getQEI <- function(x){
            res <- -qEI(matrix(x,nrow = pointsEI), model)
            if(is.nan(res)) res <- 0
            return(res)
        }
        result <- optimDE(,fun = getQEI, rep(lower,pointsEI),rep(upper,pointsEI),optimizerControl)$xbest
        newX <- matrix(result,nrow = pointsEI)
    }else{
        print("singleCore SBO (PV) will be combined with ES")
        print("Creating a point with SBO (PV)")
        
        newX <- matrix(getPVPointDACE(xAll, yAll, lower, upper),ncol = length(lower))
    }
    
    return(newX)
}

createESPoints <- function(lastESAmount, currentAmntES, xEI, yEI, esList, fun){
    print(paste("Creating", currentAmntES, "Points with ES"))
    
    ## TODO Method for restarting old bad ES is missing
    
    ## If necessary create new ES
    if(currentAmntES > length(esList)){
        esToStart <- currentAmntES - length(esList)
        
        xEI <- xEI[order(yEI),]
        yEI <- yEI[order(yEI)]
        
        for(i in 1:esToStart){
            xStart <- xEI[i,]
            yStart <- yEI[i]
            
            print(paste("Starting a new ES at",xStart," Value:", yStart, collapse=" "))
            esList[[length(esList) + 1]] <- list(x=xStart,y=yStart)
        }
    }
    
    ## Sort ES run only the best ones
    esYResults <- unlist(lapply(esList, function(l){return(l$y)}))
    esList[order(esYResults)]
    
    ## Run all ES
    for(i in 1:currentAmntES){
        esList[[i]] <- do1plus1EsIter(esList[[i]], fun)
    }
    return(esList)
}

getFeatures <- function(X, y){
    X <- as.matrix(X)
    fObj <- createFeatureObject(X = X, y = y)
    ctrl <- list(allow_cellmapping = FALSE, show_progress = F, 
                 subset = c("basic","ela_meta",
                            "ic"))
    features <- data.frame(calculateFeatures(fObj, control = ctrl))
    return(features)
}

getESAmnt <- function(algoBudget, remainingBudget, batchSize){
    progress <- (algoBudget-remainingBudget)/algoBudget
    progressSwitches <- seq(from = 0, to = 1, by = 1/batchSize)
    return(max(sum(progress > progressSwitches)-1,0))
}

################# 1+1es
solver <- function(fun,lower,upper,solverParameterList){
    tfun <- function(x){
        apply(x,1,fun)
    }
    
    print(paste("Lower bounds:",paste(lower, collapse = " ")))
    print(paste("maxBatchSize:",maxBatchSize))
    
    initialDesignSize <- length(lower) * 2 * maxBatchSize
    initialDesign <- designLHD(x = NULL,
                               lower = lower, 
                               upper = upper, 
                               control = list(
                                   size = initialDesignSize,
                                   retries=1000))
    xAll <- initialDesign
    yAll <- tfun(xAll)
    xMBO <- NULL
    
    print(paste("Total budget:",budget))
    remainingBudget <- budget - 2 * length(lower) * getEvalCost(maxBatchSize)
    print(paste("Budget after initDesign:",remainingBudget))
    algorithmBudget <- remainingBudget
    
    if(maxBatchSize == 1){
        print("BatchSize is 1, using singleCore SBO (PV)")
        for(i in 1:remainingBudget){
            print("Creating a point with SBO (PV)")
            print(yAll)
            
            newX <- matrix(getPVPointDACE(xAll, yAll, lower, upper),ncol = length(lower))
            newY <- tfun(newX)
            xAll <- rbind(xAll, newX)
            yAll <- c(yAll,newY)
        }
    }else{
        ### Start Multi-Local EI ####################
        print("MaxBatchSize is > 1, using multi-localEI")
        xEI <- xAll
        yEI <- yAll
        lastESAmount <- 0
        esList <- list()
        batchSizeDF <- NULL
        while(remainingBudget >= 1){
            ## Calculate features and predict best batch size based on them
            features <- getFeatures(xAll, yAll)
            batchSize <- predictCores(features, maxBatchSize)
            print(paste("Prediction says batchSize=",batchSize))
            batchSizeDF <- rbind(batchSizeDF, data.frame("time" = budget - remainingBudget, 
                                                         "batchSize" = batchSize))
            
            ## If not enough budget is left, batchSize might need to be reduced:
            batchSize <- min(batchSize, 
                             max(curveDf$batchSize[curveDf$budgetTakenPerIter<= remainingBudget]))
            
            ## Reduce remaining budget by the respective amount
            remainingBudget <- remainingBudget - getEvalCost(batchSize)
            print(paste("Budget was reduced by",getEvalCost(batchSize),"to",remainingBudget))

            ## get the desired amount of ES cores
            currentAmntES <- getESAmnt(algorithmBudget, remainingBudget, batchSize)
            amntEI <- batchSize - currentAmntES
            
            print(paste("the algorithm will use", currentAmntES, "ES evals and", amntEI, "EI evals"))
            
            ## run ESs
            newESX <- NULL
            newESy <- NULL
            if(currentAmntES > 0){
                esList <- createESPoints(lastESAmount, currentAmntES, xEI, yEI, esList, fun)
                for(i in 1:currentAmntES){
                    newESX <- rbind(newESX, esList[[currentAmntES]]$x1)
                    newESy <- c(newESy, esList[[currentAmntES]]$y1)
                }
            }
            
            ## create EI Cores
            newXEI <- createEIPoints(xAll, yAll, lower, upper, amntEI)
            newYEI <- tfun(newXEI)
            
            ## add EI results to result list
            xEI <- rbind(xEI, newXEI)
            yEI <- c(yEI, newYEI)
            
            print("new EI points:")
            print(newXEI)
            print(newYEI)
            print("new ES points:")
            print(newESX)
            print(newESy)
            
            
            ## add all results to overall list
            xAll <- rbind(xAll, newESX, newXEI)
            yAll <- c(yAll, newESy, newYEI)
        }
        csvName <- paste0(experimentPath,"/",paste("14AutoMultiLocalEI",paste(args,collapse="_"),sep="_"),".csv")
        write.csv(batchSizeDF, csvName)
    }
}

wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("14AutoMultiLocalEI",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)