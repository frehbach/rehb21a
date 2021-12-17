library(callr)

getBaseObjFun <- function(functionID, nDim, iid){
    #fileName <- paste0("simulatedFunctions/cobbsSim_",functionID,"_",nDim,"_",iid,".rds")
    fileName <- paste0("simulatedFunctions/cobbsSim_",functionID,"_",nDim,"_",1,".rds") ## Trying to run it only on the very first simulation
    funC <- readRDS(fileName)$simulation[[1]]
    function(x){
        funC(matrix(x, nrow=1))
    }
}

getBBOBWrappedFun <- function(functionID, algoName, experimentPath, nDim = 2, iid = 1){
    print(paste("Generating objective function: functionID:", functionID, " ,nDim:", nDim, " ,iid:",iid))
    if(functionID < 25){
        fun <- smoof::makeBBOBFunction(dimensions = nDim, fid = functionID, iid = iid)
        
        rList <- list(
            "fun" = wrapToBBOBFunction(fun = fun, 
                                              functionID = functionID, 
                                              nDim = nDim, 
                                              algoID = algoName,
                                              instanceID = iid, 
                                              experimentPath = experimentPath),
            "lower" = as.numeric(smoof::getLowerBoxConstraints(fun)),
            "upper" = as.numeric(smoof::getUpperBoxConstraints(fun))
        )
        return(rList)
    }
    fun <- smoof::makeBBOBFunction(dimensions = nDim, fid = functionID-24, iid = iid)
    list(
        "fun" = wrapToBBOBFunction(fun = getBaseObjFun(functionID-24, nDim, iid), 
                                          functionID = functionID, 
                                          nDim = nDim, 
                                          algoID = algoName,
                                          instanceID = iid, 
                                          experimentPath = experimentPath),
        "lower" = as.numeric(smoof::getLowerBoxConstraints(fun)),
        "upper" = as.numeric(smoof::getUpperBoxConstraints(fun))
    )
}
