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

################# 1+1es
solver <- function(fun,lower,upper,solverParameterList){
    tfun <- function(x){
        apply(x,1,fun)
    }
    
    print(paste("Lower bounds:",paste(lower, collapse = " ")))
    print(paste("BatchSize:",batchSize))
    
    initialDesignSize <- length(lower) * 2 * batchSize
    initialDesign <- designLHD(x = NULL,
                               lower = lower, 
                               upper = upper, 
                               control = list(
                                   size = initialDesignSize))
    xAll <- initialDesign
    yAll <- tfun(xAll)
    
    print(paste("Total Amount of iterations:",maxIters))
    itersToDo <- maxIters-2*length(lower)
    print(paste("Iterations after initDesign:",itersToDo))
    
    if(batchSize == 1){
        print("BatchSize is 1, using singleCore SBO (PV)")
        for(i in 1:itersToDo){
            print("Creating a point with SBO (PV)")
            print(yAll)
            
            newX <- matrix(getPVPointDACE(xAll, yAll, lower, upper),ncol = length(lower))
            newY <- tfun(newX)
            xAll <- rbind(xAll, newX)
            yAll <- c(yAll,newY)
        }
    }else{
        print("BatchSize is > 1, using multi-localEI")
        
        xEI <- xAll
        yEI <- yAll
        
        switchAfter <- round(itersToDo / batchSize,0)
        print(paste("Method will be switched after:",switchAfter,"iterations"))
        
        esList <- list()
        
        subIters <- 0
        switchPointsReached <- 0
        for(i in 1:itersToDo){
            ## This loop will do the remaining amount of iterations in the algorithm
            if(subIters == switchAfter){
                print("Switch point reached!")
                switchPointsReached <- switchPointsReached + 1
                subIters <- 0
            }
            
            pointsEI <- batchSize - switchPointsReached
            pointsES <- batchSize - pointsEI
            
            optimizerControl = list(funEvals = as.integer(log(length(lower)) * 1000), 
                                    populationSize = 5 * length(lower))
            
            if(pointsEI > 1){
                ## EI should be used to create multiple points in this iteration
                print(paste("Creating", pointsEI, "Points with q-EI"))
                
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
                newY <- tfun(newX)
                print(newX)
                
                xEI <- rbind(xEI, newX)
                yEI <- c(yEI,newY)
            }else{
                print("singleCore SBO (PV) will be combined with ES")
                print("Creating a point with SBO (PV)")
                
                newX <- matrix(getPVPointDACE(xAll, yAll, lower, upper),ncol = length(lower))
                newY <- tfun(newX)
            }
            
            if(pointsES > 0){
                print(paste("Creating", pointsES, "Points with ES"))
                oldESLength <- length(esList)
                if(pointsES > oldESLength){
                    ## A new ES has to be started
                    bestEIIndex <- which.min(yEI)
                    
                    xStart <- xEI[bestEIIndex,]
                    yStart <- yEI[bestEIIndex]
                    
                    print("Starting a new ES at")
                    print(xStart)
                    
                    esList[[pointsES]] <- list(x=xStart,y=yStart)
                    esList[[pointsES]] <- do1plus1EsIter(esList[[pointsES]], fun)
                }
                if(oldESLength > 0){
                    for(j in 1:oldESLength){
                        esList[[j]] <- do1plus1EsIter(esList[[j]], fun)
                    }
                }
                if(length(esList) > 0){
                    for(j in 1:length(esList)){
                        xAll <- rbind(xAll, esList[[j]]$x1)
                        yAll <- c(yAll,esList[[j]]$y1)
                    } 
                }
            }
            
            xAll <- rbind(xAll, newX)
            yAll <- c(yAll,newY)
            
            subIters <- subIters + 1
        }
    }
}

wrapped <- getBBOBWrappedFun(functionID = funID, 
                             algoName = paste("13MultiLocalEI",paste(args,collapse="_"),sep="_"), 
                             experimentPath = experimentPath,
                             nDim = nDim,
                             iid = seed)

start_time<-Sys.time()
solver(wrapped$fun, wrapped$lower, wrapped$upper)
end_time<-Sys.time()
print("Time taken: \n")
print(end_time-start_time)