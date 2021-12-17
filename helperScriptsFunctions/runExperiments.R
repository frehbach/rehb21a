###
### Main File for starting the Experiments.
### There is a small test set of experiments, with reduced budget, repeats, etc. 
### Use this test set (FULL_EXPERIMENTS <- F) to see that all the codes etc. are running correctl on your system. 
### Once that is working you can try to run the full experiment set (FULL_EXPERIMENTS <- T)
###

require(dplyr)

### RUN Parameters #########################################
# 1) SEED
# 2) Function ID (1-24 BBOB, >24 = Simulated functions)
# 3) Algo ID
# 4) nDim
# 5) base budget
# 6) batchSize

##
## Function for starting a specific set of experiments
runExperiments <- function(idList, startFile, filterDir = "exAuto2804"){
    configs <- expand.grid(idList)
    print(paste("Starting Experiments on:", startFile))
    dir.create("slurmOut", showWarnings = FALSE)
    
    print(paste0("Filtering Experiments. Before filtering: ", nrow(configs)))
    if(file.exists(paste0(filterDir,".rds"))){
        
        df <- readRDS(paste0(filterDir,".rds"))
        df <- unique(df %>% select(one_of(names(configs))))
        for(i in 1:ncol(configs)){
            df[,i] <- as.numeric(as.character(df[,i]))
            configs[,i] <- as.numeric(as.character(configs[,i]))
        }
        configs <- anti_join(configs, df)
    }
    print(paste0("Finished filtering, experiments remaining: ", nrow(configs)))
    startSingle <- function(r){
        print(r)
        system(paste0("cp helperScriptsFunctions/runSingleCore.slurm delete.slurm"))
        write(paste("/opt/software/R/R-current/bin/Rscript", startFile ,paste(r, collapse=" ")),file="delete.slurm",append=TRUE)
        system("sbatch delete.slurm")
        system("rm delete.slurm") 
    }
    if(nrow(configs) > 0){
        apply(configs, 1, startSingle) 
    }
    
    print("done!")
}

runAlgorithmByID <- function(paramList){
    for(id in paramList$algoID){
        localParamList <- paramList
        localParamList$algoID <- id
        rFileByID <- list(
            "1" = "algorithm_01_randomSearch.R",
            "2" = "algorithm_02_SPOT.R",
            "3" = "algorithm_03_SPOTDIRECT.R",
            "4" = "algorithm_04_OnlyLocal.R",
            "5" = "algorithm_05_mlrMBO.R",
            "6" = "algorithm_06_qEI.R",
            "7" = "algorithm_07_IPI.R",
            "8" = "algorithm_08_CMAES.R",
            "9" = "algorithm_09_DE.R",
            "10" = "algorithm_10_BOBYQA.sR",
            "12" = "algorithm_12_NEWUOA.R",
            "13" = "algorithm_13_MultiLocalEI.R",
            "14" = "algorithm_14_AutoMultiLocalEI.R"
        )
        runExperiments(localParamList,rFileByID[[as.character(id)]])
    }
}

