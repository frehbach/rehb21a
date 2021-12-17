#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
getAmountDominations = function(df){
    nDoms <- rep(0,length(unique(c(as.character(df[,1]),as.character(df[,2])))))
    names(nDoms) <- unique(c(as.character(df[,1]),as.character(df[,2])))
    for(i in 1:nrow(df)){
        if(df[i,3] == 2){
            nDoms[[as.character(df[i,1])]] <- nDoms[[as.character(df[i,1])]] + 1
        }
        if(df[i,3] == 1){
            nDoms[[as.character(df[i,2])]] <- nDoms[[as.character(df[i,2])]] + 1
        }
    }
    nDoms
}

#' Title
#'
#' @param melted
#'
#' @return
#' @export
#'
#' @examples
getRanksPValRanking = function(melted){
    finalRanks = rep(0,length(unique(c(as.character(melted[,1]),as.character(melted[,2])))))
    names(finalRanks) <- unique(c(as.character(melted[,1]),as.character(melted[,2])))
    while(min(finalRanks) == 0){
        if(nrow(melted) == 0){
            finalRanks[which(finalRanks == 0)] <- max(finalRanks) + 1
        }else{
            dom <- getAmountDominations(melted)
            best <- which(dom == 0)
            finalRanks[which(finalRanks == 0)][best] <- max(finalRanks) + 1

            melted <- melted[-which((melted$Var1 %in% names(best)) | (melted$Var2 %in% names(best))), ,drop=F]
        }
    }
    finalRanks
}

#' rankDataByPostHocKruskalConover
#'
#' Use posthoc.kruskal.conover.test to statistically rank the given methods.
#'
#' @param results vector of endresults
#' @param methods vector of strings, method that was used to create the results
#'
#' @return list with ranks
#' @import PMCMR
#' @export
rankDataByPostHocKruskalConover = function(results, methods){
    ## Do statistical test
    kct <- PMCMR::posthoc.kruskal.conover.test(results,methods)
    pvals <- kct$p.value

    ## Classify p - values as No_Effect (-1) or Effect (0)
    pvals[which(is.na(pvals))] <- -1
    pvals[which(pvals > 0.05)] <- -1
    pvals[which(pvals > -1)] <- 0

    ## rank the results
    ranks <- rank(results)

    ## aggregate by methods
    aggrRanks <- aggregate(ranks~methods,FUN = mean)
    ranks <- aggrRanks$ranks
    names(ranks) <- aggrRanks$methods

    #Melt into data.frame
    melted <- reshape2::melt(pvals)

    ## For those rows which have an effect see if something dominates (1) or is being dominated (2)
    rowsWithEffect <- which(melted$value == 0)
    for(i in rowsWithEffect){
        if(ranks[[as.character(melted$Var1[i])]] > ranks[[as.character(melted$Var2[i])]]){
            melted$value[i] <- 2
        }else{
            melted$value[i] <- 1
        }
    }

    #Calculate the ranks from the melted domination data.frame
    ranksOfResults <- getRanksPValRanking(melted)
    ranksOfResults
}

#' Title
#'
#' @param df
#' @param atIter
#'
#' @return
#' @export
#'
#' @examples
getBests <- function(df, atIter = Inf) {
    names <- unique(df$algoName)
    final <- NULL
    for(n in names){
        for(r in unique(df$replication[which(df$algoName == n)])){
            subDF <- df[which(df$algoName == n &
                                  df$replication == r),]
            subDF <- subDF[which(subDF$iteration < atIter),]
            minimum <- min(subDF$iterBest)
            final <- rbind(final, subDF[which(subDF$iterBest == minimum)[1],])
        }
    }
    final
}


#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
getRankDF <- function(df) {
    problems <- unique(df$problemName)
    algos <- unique(df$algoName)

    rankTable = NULL
    for(f in problems){
        subData <- df[which(df$problemName == f),]
        subData <- getBests(subData)
        subData$algoName <- as.factor(subData$algoName)
        ranks = rankDataByPostHocKruskalConover(subData$iterBest,subData$algoName)
        rankTable <- rbind(rankTable, ranks)
    }
    rownames(rankTable) <- problems
    rankTable
}

#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
rankAtIters <- function(df, atIters) {
    getRanksAtIter <- function(maxIter){
        getRankDF(df[which(df$iteration < maxIter),])
    }
    lapply(atIters,getRanksAtIter)
}
