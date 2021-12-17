###
### This file collects functions that can be used to load and analyze the results of BBOB benchmark results
###

library(gridExtra)
#library(plyr)
library(dplyr)
require(ggplot2)

negate_match_df <- function (y, x, on = NULL) 
{
    if (is.null(on)) {
        on <- intersect(names(x), names(y))
        #message("Matching on: ", paste(on, collapse = ", "))
    }
    keys <- join.keys(x, y, on)
    x[!(keys$x %in% keys$y), , drop = FALSE]
}

## Give 1 or multiple lists of settings. Checks if all these settings are present in a given data frame.
## If not it returns the configs which should be rerun
checkForMissingExperiments <- function(listsOfSettings, data){
    for(i in 1:ncol(data)){
        if(typeof(data[,i])=="factor"){
            data[,i] <- as.character(data[,i])
        }
    }
    allMatches <- NULL
    for(l in listsOfSettings){
        requiredSettings <- expand.grid(l)
        for(i in 1:ncol(requiredSettings)){
            requiredSettings[,i] <- as.character(requiredSettings[,i])
        }
        allMatches <- bind_rows(allMatches,negate_match_df(data, requiredSettings))
    }
    allMatches
}


##
## Read one or multiple experiment data folders and save the results into an rds file with the same name
readFoldersAndWriteRDS <- function(dirs, variablesInAlgoName, 
                                   algorithmIDsToReplace,
                                   algorithmNamesToInsert){
    for(dir in dirs){
        df <- readBBOB(dir)
        df <- applyFrNamingScheme(df)
        df <- df %>% select(-"seed")
        df <- appendVariablesByAlgoName(df, variablesInAlgoName)
        df$algName <- replaceConfigVector(df$algoID, algorithmIDsToReplace, algorithmNamesToInsert)
        df <- addVariablesToAlgoName(df, "batchSize")
        df$nDim <- as.numeric(as.character(df$nDim))
        df <- unique(df)
        saveRDS(df,paste0(dir,".rds"))
    }
}

##
## In order to figure out the parameters with which an algorithm was run on bbob, the algorithm name is parsed.
## This assumes that the algorithm name was appended with all run parameters before the algorithm was run.
## Each parameter is seperated with a '_'. Keeping the order of the parameters is required!
cutParametersFromString <- function(strName, 
                                    paramNameList = c("seed","budget","functionID","infillID","nDim")){
  require(stringr)
  paramValues <- str_split(strName,"_")[[1]][-1]
  names(paramValues) <- paramNameList
  return(as.data.frame(t(paramValues), stringsAsFactors=F))
}

## 
## Append a data.frame with additional columns regarding the run parameters of an algorithm.
## The lines are filled based on the the parameters that can be parsed from the algorithm name via 'cutParametersFromString'
appendVariablesByAlgoName <- function(df, paramNameList = c("seed","budget","functionID","infillID","nDim")){
  require(dplyr)
  print("Starting to append variables")
  cutParams <- function(str){cutParametersFromString(str, paramNameList)}
  additionalVars <- bind_rows(pbapply::pblapply(df$algName,cutParams))
  additionalVars[,which(names(additionalVars) %in% names(df))] <- NULL
  return(cbind(df,additionalVars))
}

## 
## Exchange numeric parameters with character strings
## for example replace infillID=1 with infillID=BP
replaceConfigVector <- function(vec, numericLevels, strLevels){
  for(i in 1:length(numericLevels)){
    nLev <- numericLevels[i]
    strLev <- strLevels[i]
    if(length(which(vec == nLev)) > 0){
      vec[which(vec == nLev)] <- strLev
    }
  }
  return(vec)
}

cleanAlgoName <- function(at){
  splits <- stringr::str_split(at$algName,stringr::fixed("_"))
  splits <- sapply(splits,function(x) x[1])
  at$algName <- splits
  return(at)
}

addVariablesToAlgoName <- function(at, varNames){
  for(v in varNames){
    at$algName <- paste(at$algName,at[[v]],sep="_")
  }
  return(at)
}

## Create a convergence plot
convergencePlotBBOB <- function(at, title = "", logScale = TRUE, scaleTo0 = FALSE){
  require(dplyr)

    if(scaleTo0){
        at$y <- at$y - min(at$y) + 0.00000001
    }
    
  #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
    at <- at %>% group_by(iteration, nDim, algName) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
    
  if(any(as.numeric(as.character(at$functionID)) > 24)){
      at$nDim[at$nDim=="8"] <- "Snake"
      at$nDim[at$nDim=="12"] <- "Gecko"
      at$nDim[at$nDim=="16"] <- "Spider"
  }
  
    p <- ggplot(at, aes(x=iteration,y=med))  + 
        geom_line(aes(color = algName)) + 
        geom_ribbon(aes(fill = algName, ymin=lower, ymax=upper), alpha=0.3) + 
        #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
        facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
        theme(text=element_text(size=16)) + 
        labs(x="evaluation", y = paste0(title,"\ny"))# +
    if(logScale){
        p <- p + scale_y_continuous(trans = 'log10')
    }
    return(p)
}

## Create a convergence plot
convergencePlotBBOBOneLine <- function(at, title = "", logScale = TRUE, scaleTo0 = FALSE){
  require(dplyr)
  
  if(scaleTo0){
    at$y <- at$y - min(at$y) + 0.00000001
  }
  
  #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
  at <- at %>% group_by(iteration, nDim, algName, functionID, batchSize) %>%  mutate(upper = summary(y)[[5]]) %>%
    mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
  
  algNames <- unlist(strsplit(at$algName,"_",T))
  at$algName <- algNames[seq(1,length(algNames),2)]
  
  at$algName[at$algName == "MultiLocalEI"] <- "ML-EI"
  at$algName[at$algName == "RandomSearch"] <- "Random"
  
  funNames <- c("1" = "Sphere", "2"="Ellipsoid", "15"="Rastrigin")
  
  at$functionID <- as.numeric(as.character(at$functionID))
  at$functionID <- paste("BBOB Function ID:",at$functionID, funNames[as.character(at$functionID)])
  at$batchSize <- paste("Batchsize:",at$batchSize)
  
  at$functionID <- factor(at$functionID, levels = c("BBOB Function ID: 1 Sphere", 
                                                    "BBOB Function ID: 2 Ellipsoid", 
                                                    "BBOB Function ID: 15 Rastrigin"))
  at <- at[order(at$functionID),]
  
  p <- ggplot(at, aes(x=iteration,y=med))  + 
    geom_line(aes(color = algName)) + 
    geom_ribbon(aes(fill = algName, ymin=lower, ymax=upper), alpha=0.3) + 
    #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
    facet_wrap(facets = vars(functionID, batchSize), nrow = 1, scales = "free") + 
    theme(text=element_text(size=16)) + 
    #scale_color_brewer(palette="Set1") + 
    labs(x="Evaluation", y = paste0(title,"Objective Function Value\n(best so far)"),
         color = "Algorithm", fill = "Algorithm")# +
  if(logScale){
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

batchSizeConvergencePlot <- function(at, title = "", curveDf = NULL, splitNCore = F){
  #browser()
    logScale = TRUE
    
    require(dplyr)

    at$y <- at$y - min(at$y) + 1e-9
    
    at$iteration <- trunc((at$iteration-1) / as.numeric(at$batchSize)) + 1
    
    #browser()
    at$unitTime <- as.numeric(at$iteration)
    at$batchSize <- as.numeric(at$batchSize)
    for(bs in curveDf$batchSize){
        at$unitTime[at$batchSize == bs] <- at$unitTime[at$batchSize == bs] / curveDf$s[curveDf$batchSize == bs]
    }
    
    at <- at[at$unitTime <= 50, ]
    
    at <- at %>% group_by(unitTime, nDim, algName) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
    
    p <- ggplot(at, aes(x=unitTime,y=med))  + 
        geom_line(aes(color = algName)) + 
        geom_ribbon(aes(fill = algName, ymin=lower, ymax=upper), alpha=0.3) + 
        #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
        theme(text=element_text(size=16)) + 
        labs(x="unitTime", y = paste0(title,"\ny")) + 
        scale_y_continuous(trans = 'log10')
    if(!splitNCore){
        p <- p + facet_wrap(facets = vars(nDim), nrow = 1, scales = "free")
    }else{
        p <- p + facet_wrap(facets = vars(nDim, batchSize), nrow = 1, scales = "free")
    }
    return(p)
}

## Create a convergence plot
convergencePlotBBOBParallel <- function(at, title = ""){
    require(dplyr)

    #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
    at <- at %>% group_by(iteration, nDim, algName, nCores) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))

    at$iteration <- trunc((at$iteration-1) / as.numeric(at$nCores)) + 1
    
    plotList <- list()
    for(n in sort(as.numeric(unique(at$nCores)))){
        localAt <- filter(at, nCores == as.character(n))
        
        plotList[[length(plotList)+1]] <- ggplot(localAt, aes(x=iteration,y=med))  + 
            geom_line(aes(color = algName)) + 
            geom_ribbon(aes(fill = algName, ymin=lower, ymax=upper), alpha=0.3) + 
            #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
            facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
            theme(text=element_text(size=16)) + 
            labs(x="iteration", y = paste0(title,"\ny")) + scale_y_continuous(trans = 'log10') + ggtitle(paste0("nCores = ",n))
    }
    plotList
    # +
    #ggtitle(title)
}

bbobBoxPlot <- function(at, title = "", iter = NULL){
  if(!is.null(iter)){
    at <- at %>% filter(iteration == iter)
  }
    
    
    ranks <- NULL
    for(dim in unique(at$nDim)){
        d <- filter(at, nDim == dim)
        
        ### only for robot functions!!!! Remove this later!!!
        factor <- 1
        if(unique(as.numeric(as.character(at$functionID))) > 24){
            factor <- -1
        }
        lRank <- rankDataByPostHocKruskalConover(d$y*factor, as.factor(d$algName))
        
        ranks <- rbind(ranks, data.frame("algName" = names(lRank), "rank" = lRank, "nDim" = dim))
    }
    rownames(ranks) <- NULL
    for(i in 1:nrow(ranks)){
        ranks$up[i] <- summary(at$y[at$algName == ranks[i,]$algName & at$nDim == ranks[i,]$nDim])[[5]]
    }
    
    if(any(as.numeric(as.character(at$functionID)) > 24)){
        at$nDim[at$nDim=="8"] <- "Snake"
        at$nDim[at$nDim=="12"] <- "Gecko"
        at$nDim[at$nDim=="16"] <- "Spider"
        ranks$nDim[ranks$nDim=="8"] <- "Snake"
        ranks$nDim[ranks$nDim=="12"] <- "Gecko"
        ranks$nDim[ranks$nDim=="16"] <- "Spider"
    }
    
    ggplot(at, aes(x=algName,y=y))  + geom_boxplot() + facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
        theme(text=element_text(size=16)) + scale_y_continuous(trans = 'log10') +
        labs(x="Algorithm", y = paste0(title, "\nEvaluation: ", iter, "\ny")) + theme(axis.text.x = element_text(angle = 90)) +
        geom_text(
            data    = ranks,
            mapping = aes(y = up, label = rank),
            hjust   = -0.1,
            vjust   = -1,
            color = "red"
        )# + ggtitle()
}

errorBarPlot <- function(at, title = ""){
    require(dplyr)
    
    #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
    at <- at %>% group_by(iteration, nDim, algName) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
    
    ## calculate positions for error bars
    desPos <- seq(21,101,30)
    ind <- unique(at$iteration)
    ind <- ind[-which(ind==20)]
    ind <- ind[-which(ind==40)]
    chosenP <- NULL
    for(p in desPos){
        chosenP <- c(chosenP,ind[which.min(abs(ind-p))])
    }

    errorBarData <- at %>% filter(iteration %in% chosenP)
    
    ggplot(at, aes(x=iteration,y=med))  + 
        geom_line(aes(color = algName)) + 
        #geom_ribbon(aes(fill = algName, ymin=lower, ymax=upper), alpha=0.3)+#+, linetype=2) + 
        geom_errorbar(data = errorBarData, aes(ymin=lower, ymax=upper, color = algName), width=15, position = "dodge") +
        #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
        facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
        theme(text=element_text(size=16)) + 
        labs(x="iteration", y = "y") + scale_y_continuous(trans = 'log10') +
        ggtitle(title)
}

rankWilcox <- function(at,...){
  y1 <- filter(at, infillID == "PM")$y
  y2 <- filter(at, infillID == "EI")$y
  
  res <- suppressWarnings(wilcox.test(x = y1, y = y2, alternative = "two.sided", paired = F, conf.int = T))
  if(is.na(res$p.value)){
    return(data.frame("dominant" = NA, "p" = res$p.value, "leading" = "PM"))
  }
  if(res$p.value > 0.05){
    return(data.frame("dominant" = NA, "p" = res$p.value, "leading" = ifelse(res$estimate < 0,"PM","EI")))
  }
  if(res$estimate < 0){
    return(data.frame("dominant" = "PM", "p" = res$p.value, "leading" = "PM"))
  }
  return(data.frame("dominant" = "EI", "p" = res$p.value, "leading" = "EI"))
}

rankTTest <- function(at,...){
  y1 <- filter(at, infillID == "PM")$y
  y2 <- filter(at, infillID == "EI")$y
  
  res <- suppressWarnings(t.test(x = y1, y = y2, alternative = "two.sided", paired = F, conf.int = T))
  if(is.na(res$p.value)){
    return(data.frame("dominant" = NA, "p" = res$p.value, "leading" = "PM"))
  }
  if(res$p.value > 0.05){
    return(data.frame("dominant" = NA, "p" = res$p.value, "leading" = ifelse(res$estimate[1] < res$estimate[2],"PM","EI")))
  }
  if(res$estimate[1] < res$estimate[2]){
    return(data.frame("dominant" = "PM", "p" = res$p.value, "leading" = "PM"))
  }
  return(data.frame("dominant" = "EI", "p" = res$p.value, "leading" = "EI"))
}

rankEIvsBP <- function(at, iter, dim, returnSum = T, useWilcox = T){
  ## Filter the data for the given iteration and dimension
  at <- filter(at, iteration == iter)
  at <- filter(at, nDim == dim)
  
  ## Group by functions
  at <- at %>% group_by(functionID)
  
  if(nrow(at) == 0){
    warning("Missing Data in Rank Test")
    return(NULL)
  }
  
  ## Calculate ranks (1-x) and subtract 1 to know amount of dominations instead of rank
  if(useWilcox){
    rank <- group_modify(at, rankWilcox) 
  }else{
    rank <- group_modify(at, rankTTest) 
  }
  
  
  if(!returnSum){
    rank$iter <- iter
    rank$nDim <- dim
    return(rank)
  }
  
  ## Aggregate sum over all functions
  sumRank <- rank %>% group_by() %>% 
    mutate(sumEI = length(which(rank$dominant == "EI"))) %>% 
    mutate(sumBP = length(which(rank$dominant == "PM"))) 
  
  ## Remove unnecessary rows and columns from aggregated df
  sumRank <- sumRank %>% filter(functionID == min(as.numeric(as.character(sumRank$functionID))))
  sumRank <- sumRank[,c(3:4)]
  
  ## Add columns with run-information
  sumRank$iter <- iter
  sumRank$nDim <- dim
  sumRank
}

## Rank EI vs BP at a given iteration and dimension
## 
rankEIvsBPMultiComparison <- function(at, iter, dim){
  appRankData <- function(at,...){
    df <- suppressWarnings(rankDataByPostHocKruskalConover(
      at$y,as.factor(at$infillID)))
    df <- t(data.frame(df))
    return(as.data.frame(df))
  }
  
  ## Filter the data for the given iteration and dimension
  at <- filter(at, iteration == iter)
  at <- filter(at, nDim == dim)
  
  ## Group by functions
  at <- at %>% group_by(functionID)
  
  if(nrow(at) == 0){
    warning("Missing Data in Rank Test")
    return(NULL)
  }
  
  ## Calculate ranks (1-x) and subtract 1 to know amount of dominations instead of rank
  rank <- group_modify(at, appRankData) 
  rank[,2:3] <- rank[,2:3]-1
  
  ## Aggregate sum over all functions
  sumRank <- rank %>% group_by() %>% 
    mutate(sumEI = sum(EI)) %>% 
    mutate(sumBP = sum(PM)) 
  
  ## Remove unnecessary rows and columns from aggregated df
  sumRank <- sumRank %>% filter(functionID == min(as.numeric(as.character(sumRank$functionID))))
  sumRank <- sumRank[,c(4:5)]
  
  ## Add columns with run-information
  sumRank$iter <- iter
  sumRank$nDim <- dim
  sumRank
}

##
## Count amount of dominations and plot over given iterations and dimensions
rankDominatedAmountsEIvsBP <- function(at, iters = unique(at$iteration), nDims = unique(at$nDim), doPlot = T){
  ## Iterate over all dimensions and iters, apply ranking function and collect result in single data.frame
  df <- NULL
  for(d in nDims){
    getSingleDF <- function(i){
      df <- rankEIvsBP(at,i,d)
      df
    }
    df <- dplyr::bind_rows(df,lapply(iters,getSingleDF))
  }
  
  names(df) <- c("EI","PM","iteration","nDim")
  
  df <- df[,c(2,1,3,4)]
  
  if(doPlot){
    ## melt the data for plotting
    melted <- reshape2::melt(df, id.vars = c("iteration","nDim"))
    names(melted) <- c("iteration","nDim","InfillCriterion","CountDominated")

    ## generate plot
    p <- ggplot(melted, aes(x=iteration, y = CountDominated)) + geom_line(aes(color = InfillCriterion)) + 
      facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
      theme(text=element_text(size=16)) + ylab("Amount of Dominations") + xlab("Iteration")
    return(p)
  }
  df
}

## Create a convergence plot for multiple cores in one graph
convergencePlotMeasurementParallel <- function(at, title = ""){
    require(dplyr)
    
    #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
    at <- at %>% group_by(iteration, nDim, algName, nCores) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
    
    #at$iteration <- trunc((at$iteration-1) / as.numeric(at$nCores)) + 1
    
    plotList <- list()
    for(alg in sort(unique(at$algName))){
        localAt <- filter(at, algName == alg)
        localAt$nCores <- factor(localAt$nCores, levels = as.character(sort(as.numeric(unique(localAt$nCores)))))
        
        plotList[[length(plotList)+1]] <- ggplot(localAt, aes(x=iteration,y=med))  + 
            geom_line(aes(color = nCores)) + 
            geom_ribbon(aes(fill = nCores, ymin=lower, ymax=upper), alpha=0.3) + 
            #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
            facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
            theme(text=element_text(size=16)) + 
            labs(x="iteration", y = paste0(title,"\ny")) + scale_y_continuous(trans = 'log10') + ggtitle(paste0("Algorithm = ",alg))
    }
    plotList
    # +
    #ggtitle(title)
}

plotMedianParallelPerformance <- function(at, title = ""){
    require(dplyr)
    
    at <- at %>% group_by(iteration, nDim, algName, nCores) %>% mutate(med = summary(y)[[3]])
    at$iteration <- trunc((at$iteration-1) / as.numeric(at$nCores)) + 1
    at <- at %>% select(-y) %>% unique()
    
    addMin <- function(df,...){
        dfMinCores <- df[which(df$nCores == min(as.numeric(df$nCores))),]
        bestMedian <- min(dfMinCores$med)
        df %>% mutate(referencePerformance = bestMedian)
    }
    at <- at %>% group_by(nDim, algName) %>% group_modify(addMin) 
    
    whenReferenceReached <- function(df,...){
        dfReachedReference <- df[which(df$med <= df$referencePerformance),]
        minIter <- min(dfReachedReference$iteration)
        df %>% mutate(referenceIteration = minIter)
    }
    at <- at %>% group_by(nDim, algName, nCores) %>% group_modify(whenReferenceReached) 
    
    at <- at[which(at$iteration == max(at$iteration)),]
   
        
    at <- at %>%
        group_by(nDim, algName, nCores) %>%
        mutate(relativePerformance = referencePerformance/med)
    at$referenceIteration[at$referenceIteration>100] <- 110
    
    
    
    #at$nCores <- factor(at$nCores, levels = as.character(sort(as.numeric(unique(at$nCores)))))
    at$nCores <- as.numeric(at$nCores)
    at <- at %>% select(nCores, referenceIteration, algName, nDim)
    
    ggplot(at, aes(x=nCores,y=referenceIteration))  + 
        geom_line(aes(color = algName, group = algName)) + 
        geom_point(aes(color = algName, group = algName)) + 
        facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
        theme(text=element_text(size=16)) + 
        labs(x="points per iteration", y = paste0(title,"\ny")) #+ scale_y_continuous(trans = 'log10')
    #}
    #plotList
    # +
    #ggtitle(title)
}

convergencePlotParallel <- function(at, title = ""){
    require(dplyr)
    
    #at <- at %>%  mutate(isSPOT = str_starts(algName,"SPOT"))
    at <- at %>% group_by(iteration, nDim, algName, nCores) %>%  mutate(upper = summary(y)[[5]]) %>%
        mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
    
    at$iteration <- trunc((at$iteration-1) / as.numeric(at$nCores)) + 1
    
    at <- at %>% filter(nCores %in% c(2,4,8,16))
    
    plotList <- list()
    for(alg in sort(unique(at$algName))){
        localAt <- filter(at, algName == alg)
        localAt$nCores <- factor(localAt$nCores, levels = as.character(sort(as.numeric(unique(localAt$nCores)))))
        
        plotList[[length(plotList)+1]] <- ggplot(localAt, aes(x=iteration,y=med))  + 
            geom_line(aes(color = nCores)) + 
            geom_ribbon(aes(fill = nCores, ymin=lower, ymax=upper), alpha=0.3) + 
            #facet_grid(rows = vars(),cols = vars(nDim), scales = "free_y") + 
            facet_wrap(facets = vars(nDim), nrow = 1, scales = "free") + 
            theme(text=element_text(size=16)) + 
            labs(x="iteration", y = paste0(title,"\ny")) + scale_y_continuous(trans = 'log10') + ggtitle(paste0("Algorithm = ",alg))
    }
    plotList
    # +
    #ggtitle(title)
}

calculateParallelEfficiency <- function(at, testIter){
    at <- at %>% group_by(iteration, nDim, algName, nCores, functionID) %>% mutate(med = median(y))
    at <- at %>% select(-y,-instanceID,-budget,-assignBudgetPerCore) %>% unique()
    
    addMin <- function(df,...){
        dfMinCores <- df[which(df$nCores == min(as.numeric(df$nCores))),]
        bestMedian <- min(dfMinCores$med)
        df %>% mutate(referencePerformance = bestMedian)
    }
    at <- at %>% group_by(nDim, algName, functionID) %>% group_modify(addMin) 
    
    whenReferenceReached <- function(df,...){
        dfReachedReference <- df[which(df$med <= df$referencePerformance),]
        df$referenceIteration <- suppressWarnings(min(dfReachedReference$iteration))
        return(df)
        #minIter <- suppressWarnings(min(dfReachedReference$iteration))
        #df %>% mutate(referenceIteration = minIter)
    }
    at <- at %>% group_by(nDim, algName, nCores, functionID) %>% group_modify(whenReferenceReached) 
    at <- at[which(at$iteration == max(at$iteration)),]
    ## Scale to max iteration to keep same scale everywhere
    
    #at$referenceIteration <- at$referenceIteration/at$iteration*100
    
    ## scale initial value up to always be 100
    scaleValuesToInitial <- function(df,...){
        if(df$algoID == 1){
            browser()
        }
        refAt2 <- df[which(df$nCores == 2),]$referenceIteration
        if(length(refAt2) == 0)return(data.frame())
        idealIter <- ceiling(df$iteration/as.numeric(df$nCores)*2)
        df$referenceIteration <- idealIter / df$referenceIteration
        #df$referenceIteration <- df$referenceIteration/refAt2*100
        #df$referenceIteration <- 1 / (df$referenceIteration*as.numeric(df$nCores)/2) * 100
        return(df)
    }
    
    at <- at %>% group_by(nDim, algName, functionID) %>% group_modify(scaleValuesToInitial)
    at$referenceIteration[at$referenceIteration>1.00] <- 1.00
    return(at)
}

tilePlotParallelEfficiency <- function(at){
    ## One plot per algorithm should be generated
    
    ## Transform evaluations into iterations based on am,ount of parallel cores
    at$iteration <- trunc((at$iteration-1) / as.numeric(at$nCores)) + 1
    at <- at %>% filter(nCores %in% c(2,4,8,16,32,64))
    iters <- unique(at$iteration)
    nDims <- unique(at$nDim)
    
    ## Es muss für jede Kombination von Dimension,Funktion,iteration,nCores eine Farbe erzeugt werden    
    df <- NULL
    for(iter in c(10,25,50,79,100)){
        #for(fun in at$functionID){
            localAt <- at[at$iteration <= iter,]
            newDf <- calculateParallelEfficiency(localAt,iter)
            df <- dplyr::bind_rows(df,newDf)   
        #}
    }
    
    df$functionID <- as.numeric(as.character(df$functionID))
    df$functionID <- factor(df$functionID, levels = c(24:1))
    df$nCores <- factor(df$nCores, levels = sort(unique(as.numeric(df$nCores)),decreasing = F))
    df$iter <- as.factor(df$iteration)
    
    pdf("tilePlots.pdf", width = 15,height = 9)
    for(alg in unique(df$algName)){
        localDf <- df %>% filter(algName == alg)
        print(ggplot(localDf, aes(nCores, functionID)) +
            geom_tile(aes(fill = localDf$referenceIteration), colour = "grey10", width = 0.9, height = 0.9, size = 0.1) +
            facet_grid(rows = vars(nDim),cols = vars(iter)) +
            scale_fill_gradientn(colours = hcl.colors(4)) + 
            theme(text=element_text(size=16)) + 
            ylab("BBOB Function ID") + xlab("nCores") + labs("fill" = "Scaling") +
            ggtitle(alg))
    }
    dev.off()
}


tilePlotFinalCores <- function(at){
  require(dplyr)
  
  at$y <- at$y - min(at$y) + 0.00000001
  at$iteration <- as.numeric(at$iteration)
  at$batchSize <- as.numeric(at$batchSize)
  
  algNames <- unlist(strsplit(at$algName,"_",T))
  at$algName <- algNames[seq(1,length(algNames),2)]
  
  at <- at %>% group_by(iteration, nDim, algName, functionID, batchSize) %>%  mutate(upper = summary(y)[[5]]) %>%
    mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
  
  getBestAlgo <- function(df,...){
    bestMedInd <- which.min(df$med)
    df$bestAlgo <- df[bestMedInd,]$algName
    return(df)
  }
  
  at <- at %>% group_by(iteration, nDim, functionID, batchSize) %>% group_modify(getBestAlgo)
  at$functionID <- factor(at$functionID, levels = 1:24)
  at$batchSize <- factor(at$batchSize)
  
  at <- at[,which(names(at) %in% c("nDim", "functionID", "batchSize", "bestAlgo"))]
  #at$bestAlgo <- factor(at$bestAlgo, levels = c(unique(at$bestAlgo), "Random"))
  
  cbp1 <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00")
  
  at$nDim <- paste("Dimensions:", at$nDim)
  
  p <- ggplot(at, aes(batchSize, functionID)) +
    geom_tile(aes(fill = at$bestAlgo), colour = "grey10", width = 0.9, height = 0.9, size = 0.1) +
    facet_grid(rows = vars(nDim),cols = vars()) +
    scale_y_discrete(breaks=c(1,5,10,15,20,24)) + 
    #scale_x_discrete(breaks=iters[seq(1, length(iters), 3)]) + 
    scale_fill_manual(values = cbp1) +
    theme(text=element_text(size=16)) + 
    ylab("BBOB Function ID") + xlab("Batch Size") + labs("fill" = "Algorithm")
  return(p)
}

tilePlotBestBatchSize <- function(at, curveDF){
  require(dplyr)
  
  at$y <- at$y - min(at$y) + 0.00000001
  at$iteration <- as.numeric(at$iteration)/as.numeric(at$batchSize)
  at$batchSize <- as.numeric(at$batchSize)
  
  algNames <- unlist(strsplit(at$algName,"_",T))
  at$algName <- algNames[seq(1,length(algNames),2)]
  
  at$unitTime <- as.numeric(at$iteration)
  at$batchSize <- as.numeric(at$batchSize)
  for(bs in curveDf$batchSize){
    at$unitTime[at$batchSize == bs] <- at$unitTime[at$batchSize == bs] / curveDf$s[curveDf$batchSize == bs]
  }
  
  at <- at[at$unitTime <= 10, ]
  subAt <- NULL
  for(b in unique(at$batchSize)){
    print(b)
    maxUnitTime <- max(at[at$batchSize ==b,]$unitTime)
    subAt <- rbind(subAt, at[at$unitTime == maxUnitTime,])
    print(nrow(subAt))
  }
  
  at <- at %>% group_by(iteration, nDim, algName, functionID, batchSize) %>%  mutate(upper = summary(y)[[5]]) %>%
    mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
  
  getBestBatchSize <- function(df,...){
    bestMedInd <- which.min(df$med)
    df$bestBatchSize <- df[bestMedInd,]$batchSize
    return(df)
  }
  
  at <- at %>% group_by(nDim, functionID, algName) %>% group_modify(getBestBatchSize)
  at$functionID <- factor(at$functionID, levels = 1:24)
  at$batchSize <- factor(at$batchSize)
  
  at <- at[,which(names(at) %in% c("nDim", "functionID", "algName", "bestBatchSize"))]
  
  at$bestBatchSize <- factor(as.character(at$bestBatchSize), levels = c("2", "4", "8", "16"))
  
  at <- unique(at)
  
  at$algName[at$algName == "MultiLocalEI"] <- "ML-EI"
  at$algName[at$algName == "RandomSearch"] <- "Random"
  
  
  cbp1 <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00")
  
  at$nDim <- paste("Dimensions:", at$nDim)
  
  p <- ggplot(at, aes(algName, functionID)) +
    geom_tile(aes(fill = bestBatchSize), colour = "grey10", width = 0.9, height = 0.9, size = 0.1) +
    facet_grid(rows = vars(nDim),cols = vars()) +
    scale_y_discrete(breaks=c(1,5,10,15,20,24)) + 
    scale_fill_manual(values = cbp1) +
    theme(text=element_text(size=16), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("BBOB Function ID") + xlab("Algorithm") + labs("fill" = "Batch Size")
  return(p)
}

tilePlotFinalCoresDomination <- function(at){
  require(dplyr)
  
  getFunRanks <- function(at){
    ranks <- NULL
    for(fun in unique(at$functionID)){
      for(dim in unique(at$nDim)){
        for(batchS in unique(at$batchSize)){
          d <- filter(at, nDim == dim, batchSize == batchS, functionID == fun)
          
          lRank <- rankDataByPostHocKruskalConover(d$y, as.factor(d$algName))
          
          ranks <- rbind(ranks, data.frame("algName" = names(lRank), "rank" = lRank, 
                                           "nDim" = dim, 
                                           "batchSize" = batchS,
                                           "funID" = fun))
        }
      }
    }
    rownames(ranks) <- NULL
    return(ranks)
  }
  allRanks <- getFunRanks(at)
  
  algNames <- unlist(strsplit(as.character(allRanks$algName),"_",T))
  allRanks$algName <- algNames[seq(1,length(algNames),2)]
  allRanks$algName[allRanks$algName == "MultiLocalEI"] <- "ML-EI"
  allRanks$algName[allRanks$algName == "RandomSearch"] <- "Random"
  
  allRanks$funID <- factor(allRanks$funID, levels = 1:24)
  
  allRanks$rank[allRanks$rank > 1] <- 2
  
  cbp1 <- c("#D3D3D3", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  allRanks$rank[allRanks$rank == 2] = "NotDominated"
  allRanks$rank[allRanks$rank == 1] = allRanks$algName[allRanks$rank == 1]
  
  allRanks$rank = factor(allRanks$rank, levels = c("NotDominated", unique(allRanks$algName)))
  
  allRanks$nDim <- paste("Dimensions:", allRanks$nDim)
  allRanks$batchSize <- paste("Batch Size:", allRanks$batchSize)
  allRanks$batchSize <- factor(allRanks$batchSize, levels = c("Batch Size: 2", 
                                                              "Batch Size: 4", 
                                                              "Batch Size: 8",
                                                              "Batch Size: 16"))
  #at$functionID <- factor(at$functionID, levels = c("Fun. ID: 1", "Fun. ID: 2", "Fun. ID: 15"))
  #at <- at[order(at$functionID),]
  
  p = ggplot(allRanks, aes(algName, funID)) +
    geom_tile(aes(fill = rank), colour = "grey10", width = 0.9, height = 0.9, size = 0.1) +
    facet_grid(rows = vars(nDim),cols = vars(batchSize)) +
    scale_y_discrete(breaks=c(1,5,10,15,20,24)) + 
    scale_fill_manual(values = cbp1) +
    theme(text=element_text(size=16), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none",
          strip.text.y.left = element_text(angle = 0)) +
    ylab("BBOB Function ID") + xlab("Algorithm") + labs("fill" = "Algorithm")
  return(p)
}


tilePlotBatches <- function(at, curveDf){
  require(dplyr)
  
  at$y <- at$y - min(at$y) + 0.00000001
  
  at$iteration <- trunc((at$iteration-1) / as.numeric(at$batchSize)) + 1
  
  at$unitTime <- as.numeric(at$iteration)
  at$batchSize <- as.numeric(at$batchSize)
  for(bs in curveDf$batchSize){
    at$unitTime[at$batchSize == bs] <- at$unitTime[at$batchSize == bs] / curveDf$s[curveDf$batchSize == bs]
  }
  
  at <- at[at$unitTime <= 50, ]
  
  at <- at %>% group_by(unitTime, nDim, algName, functionID) %>%  mutate(upper = summary(y)[[5]]) %>%
    mutate(lower = summary(y)[[2]]) %>% mutate(med = summary(y)[[3]]) %>% mutate(mmin = mean(y))
  
  at$unitTime <- round(at$unitTime)
  
  getBestAlgo <- function(df,...){
    bestMedInd <- which.min(df$med)
    df$bestAlgo <- df[bestMedInd,]$algName
    return(df)
  }
  
  at <- at %>% group_by(unitTime, nDim, functionID) %>% group_modify(getBestAlgo)
  at$functionID <- factor(at$functionID, levels = 1:24)
  
  p <- ggplot(at, aes(unitTime, functionID)) +
    geom_tile(aes(fill = at$bestAlgo), colour = "grey10", width = 0.9, height = 0.9, size = 0.1) +
    facet_grid(rows = vars(nDim),cols = vars()) +
    scale_y_discrete(breaks=c(1,5,10,15,20,24)) + 
    #scale_x_discrete(breaks=iters[seq(1, length(iters), 3)]) + 
    #scale_fill_manual(values = c("red","blue","grey90")) + 
    theme(text=element_text(size=16)) + 
    ylab("BBOB Function ID") + xlab("Unit Time") + labs("fill" = "Algorithm")
  return(p)
}