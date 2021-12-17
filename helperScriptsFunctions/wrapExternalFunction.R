
initFilesAndFolders <- function(functionID, nDim, algoID, instanceID, experimentPath){
  if(stringr::str_ends(experimentPath,stringr::fixed("/"))){
    experimentPath <- stringr::str_trunc(experimentPath, stringr::str_length(experimentPath)-1, "right", ellipsis = "")
  }
  path <- paste0(experimentPath, "/", getSubFolderName(algoID), "/")
  dir.create(path, showWarnings = FALSE, recursive = T)

  lines <- c(
    paste0("suite = 'bbob', funcId = ", functionID,
           ", DIM = " ,nDim , ", Precision = 1.000e-08, algId = '", algoID,
           "', coco_version = '2.3.1', logger = 'FR-bbob', data_format = 'FR-bbob'"),
    "%",
    paste0("data_f1/bbobexp_f1_DIM2_i1.dat, ", instanceID, ":100|5.3e+01")
  )

  infoPath <- paste0(path, "bbobexp_f", functionID, "_i", instanceID, ".info")
  fileConn<-file(infoPath)
  writeLines(lines, fileConn)
  close(fileConn)

  dataPath <- paste0(path, "data_f", functionID)
  dir.create(dataPath, showWarnings = FALSE, recursive = T)

  tdatLine <- "% f evaluations | g evaluations | best noise-free fitness - Fopt + sum g_i+ | measured fitness | best measured fitness or single-digit g-values | x1 | x2...\n"
  tdatPath <- paste0(dataPath, "/", getTDatFileName(functionID, nDim, instanceID))
  cat(tdatLine, file = tdatPath)
  return(tdatPath)
}

getSubFolderName <- function(algoID){
  return(paste0(algoID,"_on_bbob_budget0100xD"))
}

getTDatFileName <- function(functionID, nDim, instanceID){
  return(paste0("bbobexp_f", functionID, "_DIM", nDim, "_i", instanceID, ".tdat"))
}

## % f evaluations | g evaluations | best noise-free fitness - Fopt (6.671000000000e+01) + sum g_i+ | measured fitness | best measured fitness or single-digit g-values | x1 | x2...
## iter | 0 | best(fun - knownOptimum) | fun | best(fun) | x1 | x2 ...

writeBBOBIter <- function(filePath, iter, bestMinKnown, y, bestY, x, printX, saveIndividuals){
  if(printX){
    line <- paste(iter, 0, bestMinKnown, y, bestY, paste0(x, collapse = " "),"\n")
  }else{
    line <- paste(iter, 0, bestMinKnown, y, bestY,"\n")
  }
  cat(line, file = filePath, append = T)

  if(saveIndividuals){
    ind <- list("x" = x, "y" = y)
    subPath <- stringr::str_remove(filePath,fixed(".tdat"))
    saveRDS(ind, file = paste0(subPath, "_eval_", iter,".rds"))
  }
}

wrapToBBOBFunction <- function(fun, functionID, nDim, algoID, instanceID = 1, experimentPath = "exdata", knownFunctionOptimum = 0,
                               printX = T, saveIndividuals = F){
  ## 'Function-Globals' for the wrapped function
  tdatPath <- initFilesAndFolders(functionID, nDim, algoID, instanceID, experimentPath)
  iteration <- 0
  bestY <- Inf

  return(
    function(x){
      y <- fun(x)
      iteration <<- iteration + 1
      if(y < bestY){
        bestY <<- y
      }
      writeBBOBIter(tdatPath, iteration, bestY - knownFunctionOptimum, y, bestY, x, printX, saveIndividuals)
      return(y)
    }
  )
}

