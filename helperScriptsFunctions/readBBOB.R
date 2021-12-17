readSingleResultFile <- function(filePath){
  lines <- readLines(filePath)
  headers <- cumsum(stringr::str_starts(lines,"%"))
  instances <- split(lines,headers)

  parseTable <- function(str){
    try({
      tbl <- read.table(text = str, sep = " ", skip = 1)
      colnames(tbl) <- c("fEvaluations","gEvaluations","bestNoiseFreeFitness","measuredFitness","bestMeasuredFitness",paste0("x", 1:(ncol(tbl)-5)))
      return(tbl)
    })
    print(paste0("Reading of file ", filePath, " failed, incomplete file!"))
    return(data.frame())
  }
  return(lapply(instances,parseTable))
}

readDataDir <- function(fileDir){
  files <- dir(fileDir)
  fileDir <- paste0(fileDir,"/",files[which(stringr::str_ends(files,"tdat"))])

  return(lapply(fileDir, readSingleResultFile))
}

getVarByCommaSplit <- function(str,varName){
  subStr <- stringr::str_split(str,stringr::fixed(paste0(varName," = ")))[[1]][2]
  return(stringr::str_replace_all(stringr::str_split(subStr,",")[[1]][1],pattern = "'",""))
}

readInfoFile <- function(filePath){
  lines <- suppressWarnings(readLines(filePath))

  getSingleEntryList <- function(lines){
    funID <- getVarByCommaSplit(lines[1],"funcId")
    algName <- getVarByCommaSplit(lines[1],"algId")
    nDim <- getVarByCommaSplit(lines[1],"DIM")

    instances <- stringr::str_split(lines[3],",")[[1]]
    instances <- instances[-1]
    instances <- unlist(stringr::str_split(instances,":"))[seq(1,length(instances)*2,2)]
    instances <- stringr::str_remove_all(instances," ")
    return(list("funID" = funID,"algName" = algName,"nDim" = nDim, "instances" = instances))
  }

  lineSets <- split(lines, ceiling(seq_along(lines)/3))
  return(lapply(lineSets, getSingleEntryList))
}

readBatch <- function(dataDir){
  require(dplyr)
  allFiles <- paste0(dataDir,"/",dir(dataDir))
  infoInd <- stringr::str_ends(allFiles,stringr::fixed(".info"))
  infoFiles <- allFiles[infoInd]
  dataFolders <- allFiles[!infoInd]

  if(length(infoFiles) < 1){
    return(NULL)
  }

  readFileAtIndex <- function(i){
    infoFile <- readInfoFile(infoFiles[i])
    tbls <- readDataDir(dataFolders[i])

    if(length(infoFile) == 0){
      warning(paste("An infoFile is corrupted!",infoFiles[i]))
      return(NULL)
    }

    readSingleEntry <- function(entry){
      infoFile <- infoFile[[entry]]
      tbls <- tbls[[entry]]

      for(j in 1:length(tbls)){
        tbls[[j]] <- tbls[[j]] %>% mutate(funID = infoFile$funID, algName = infoFile$algName,
                                          nDim = infoFile$nDim,instance = infoFile$instances[j])
      }
      return(dplyr::bind_rows(tbls))
    }
    return(dplyr::bind_rows(lapply(1:length(infoFile),readSingleEntry)))
  }

  return(dplyr::bind_rows(lapply(1:length(infoFiles),readFileAtIndex)))
}

readBBOB <- function(dataDir, filterString = NULL){
  foundFolders <- dir(dataDir)
  if(!is.null(filterString)){
    foundFolders <- foundFolders[stringr::str_detect(foundFolders, stringr::fixed(filterString))]
  }
  batchFolders <- paste0(dataDir,"/",foundFolders)
  print("Starting to read BBOB data.")
  print(paste("Found a total of", length(batchFolders), "batches to read"))
  resultDF <- dplyr::bind_rows(pbapply::pblapply(batchFolders,readBatch))
  return(resultDF)
}

applyFrNamingScheme <- function(df, keepX = FALSE, keepY = FALSE){
  newDf <- data.frame("seed" = df$instance, "functionID" = df$funID, "nDim" = df$nDim, "iteration" = df$fEvaluations,
                      "y" = df$bestNoiseFreeFitness, "algName" = df$algName)
  if(keepX){
      newDf <- cbind(newDf, df[,startsWith(names(df),"x")])
  }
  if(keepY){
    newDf <- cbind(newDf, df[,startsWith(names(df),"measuredFitnes")])
  }
  newDf
}
