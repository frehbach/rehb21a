##
## Main File for analysing the experiment data and creating plots.
## Running this file should produce .pdf Files with the plots similar to the papers supplementary material.
## Running this Code on the data generated with (FULL_EXPERIMENTS <- F) results in some warnings, these can
## be ignored.
##

library(ggplot2)
library(plotly)
library(parallel)

setBBOBFunToZero <- function(df){
    getZeroed <- function(cList){
        functionID <- as.numeric(cList[1])
        nDim <- as.numeric(cList[2])
        seed <- as.numeric(cList[3])
        opt <- smoof::getGlobalOptimum(smoof::makeBBOBFunction(dimensions = as.numeric(as.character(nDim)), 
                                                               fid = as.numeric(as.character(functionID)), 
                                                               iid = as.numeric(as.character(seed))))$value
        subDf <- df[df$functionID == functionID &
                       df$nDim == nDim &
                       df$seed == seed,]
        subDf$y <- subDf$y - opt
        subDf
    }
    
    configs <- expand.grid(list(unique(df$functionID), unique(df$nDim), unique(df$seed)))
    configs <- split(configs, seq(nrow(configs)))
    df <- bind_rows(mclapply(configs, getZeroed,mc.cores = 10))
    
    # for(functionID in unique(df$functionID)){
    #     for(nDim in unique(df$nDim)){
    #         for(seed in unique(df$seed)){
    #             opt <- smoof::getGlobalOptimum(smoof::makeBBOBFunction(dimensions = as.numeric(as.character(nDim)), 
    #                                                                    fid = as.numeric(as.character(functionID)), 
    #                                                                    iid = as.numeric(as.character(seed))))$value
    #             print(opt)
    #             df[df$functionID == functionID &
    #                    df$nDim == nDim &
    #                    df$seed == seed,]$y <- df[df$functionID == functionID &
    #                                                  df$nDim == nDim &
    #                                                  df$seed == seed,]$y - opt
    #         }
    #     }
    # }
    return(df)
}


createAllPlots <- function(plotData, strName, plotType = 1, curve = NULL){
    if(plotType == 1 || plotType == 5){
        pdf(strName, width = 15,height = 4)
    }else if(plotType == 6){
        pdf(strName, width = 4.5,height = 6)
    }else if(plotType == 3){
        pdf(strName, width = 8.5,height = 6)
    }else{
        pdf(strName, width = 10,height = 6)
    }
    
    if(plotType == 3){
        plotData$batchSize <- as.numeric(as.character(plotData$batchSize))
        plotData$budget <- as.numeric(as.character(plotData$budget))
        plotData <- plotData[plotData$iteration == (plotData$budget * plotData$batchSize),]
        print(tilePlotFinalCores(plotData))
    }else if(plotType == 4){
        plotData$batchSize <- as.numeric(as.character(plotData$batchSize))
        plotData$budget <- as.numeric(as.character(plotData$budget))
        plotData <- plotData[plotData$iteration == (plotData$budget * plotData$batchSize),]
        print(tilePlotFinalCoresDomination(plotData))
    }else if(plotType == 5){
        plotData$batchSize <- as.numeric(as.character(plotData$batchSize))
        plotData$budget <- as.numeric(as.character(plotData$budget))
        
        print(convergencePlotBBOBOneLine(plotData, title = ""))
    }else if(plotType == 6){
        plotData$batchSize <- as.numeric(as.character(plotData$batchSize))
        plotData$budget <- as.numeric(as.character(plotData$budget))
        #plotData <- plotData[plotData$iteration == (plotData$budget * plotData$batchSize),]
        print(tilePlotBestBatchSize(plotData, curve))
    }else{
        for(i in 1:24){
            at <- filter(plotData,functionID %in% i)
            if(nrow(at) == 0){
                next
            }
            
            at$batchSize <- as.numeric(as.character(at$batchSize))
            at$budget <- as.numeric(as.character(at$budget))
            
            batchSizes <- unique(at$batchSize)
            for(b in batchSizes){
                subData <- at %>% filter(batchSize == b)
                if(plotType == 1){
                    subData <- subData[subData$iteration <= (subData$budget * subData$batchSize),]
                    print(convergencePlotBBOB(subData, title = paste("FunID:",i, "- BatchSize:",b)))
                }
                if(plotType == 2){
                    subData <- subData[subData$iteration == (subData$budget * subData$batchSize),]
                    print(bbobBoxPlot(subData, title = paste("FunID:",i, "- BatchSize:",b)))
                }
            }
        }
    }
    dev.off()
}

createListOfPlots <- function(plotData, strName, plotType = 1, curve = NULL){
    plotList <- list()
    for(i in 1:24){
        at <- filter(plotData,functionID %in% i)
        if(nrow(at) == 0){
            next
        }
        if(plotType == 1){
            plotList[[length(plotList)+1]] <- batchSizeConvergencePlot(at, title = paste("FunID:",i), curve)  
        }
    }
    plotList
}

loadResultFolder <- function(dirToRun){
    readFoldersAndWriteRDS(dirs = c(dirToRun), 
                           variablesInAlgoName = c("seed","functionID","algoID","nDim","budget","batchSize"),
                           algorithmIDsToReplace = c(1,2,3,4,5,6,7,8,9,10,12,13,14),
                           algorithmNamesToInsert = c("RandomSearch",
                                                      "SPOT",
                                                      "DIRECT",
                                                      "1+1ES",
                                                      "moimbo",
                                                      "Q-EI",
                                                      "IPI",
                                                      "CMAES",
                                                      "DE",
                                                      "BOBYQA",
                                                      "NEWUOA",
                                                      "MultiLocalEI",
                                                      "AutoMultiLocalEI"))
}
