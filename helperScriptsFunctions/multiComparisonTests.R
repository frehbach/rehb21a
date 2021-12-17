require(PMCMR)

getRanks <- function(pVals, ranks){
    initialNames <- names(ranks)
    finalRanks <- rep(0,length(ranks))
    names(finalRanks) <- initialNames
    while(min(finalRanks) == 0)
    {
        localRanks <- ranks[which(finalRanks == 0)]
        best <- names(localRanks)[which.min(localRanks)]
        
        nonDominated <- names(which(pVals[which(rownames(pVals)==best),]==1))
        nonDominated <- nonDominated[finalRanks[nonDominated] == 0]
        bestNames <- c(best, nonDominated)
        ind <- which(initialNames %in% bestNames)
        finalRanks[ind] <- max(finalRanks) + 1
    }
    finalRanks
}

rankDataByPostHocKruskalConover = function(results, methods)
{
    res <- kruskal.test(results,methods)
    if(is.nan(res$p.value)){
        ranks <- rep(1,length(levels(methods)))
        names(ranks) <- levels(methods)
        return(ranks)
    }
    if(res$p.value > 0.05){
        ranks <- rep(1,length(levels(methods)))
        names(ranks) <- levels(methods)
        return(ranks)
    }
    kct = posthoc.kruskal.conover.test(results,methods)
    pvals = kct$p.value
    ranks = rank(results)
    temp = aggregate(ranks~methods,FUN = median)
    pvals <- cbind(pvals, rep(NA,nrow(pvals)))
    pvals <- rbind(c(NA,pvals[,1]),pvals)
    rownames(pvals)[1] <- colnames(pvals)[1]
    for(i in 1:nrow(pvals)){
        for(j in 1:ncol(pvals)){
            if(is.na(pvals[i,j]))
                pvals[i,j] = pvals[j,i]
        }
    }
    for(i in 1:nrow(pvals)){
        for(j in 1:ncol(pvals)){
            if(is.na(pvals[i,j]))
                pvals[i,j] = 1
            if(pvals[i,j] > 0.05)
            {
                pvals[i,j] = 1
            }else{
                pvals[i,j] = 0
            }
        }
    }
    ranks <- temp$ranks
    names(ranks) <- levels(methods)
    finalRanks = getRanks(pvals, ranks)
    names(finalRanks) <- levels(methods)
    finalRanks
}
