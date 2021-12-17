library(ggplot2)
library(plotly)

curveDf <- data.frame("batchSize" = c(16,8,4,2,1), "s" = c(1,1.7,3,5.3,7.8))

files.sources = dir("helperScriptsFunctions/")
files.sources <- paste0("helperScriptsFunctions/",files.sources[endsWith(files.sources, ".R")])
invisible(sapply(files.sources, source))

loadResultFolder("expResults/exp1_RandomSearch")
loadResultFolder("expResults/exp1_ES")
loadResultFolder("expResults/exp1_CMAES")
loadResultFolder("expResults/exp1_IPI")
loadResultFolder("expResults/exp1_MOIMBO")
loadResultFolder("expResults/exp1_MultiLocalEI")
loadResultFolder("expResults/exp1_QEI")

df <- readRDS("expResults/exp1_CMAES.rds")
df <- rbind(df, readRDS("expResults/exp1_ES.rds"))
df <- rbind(df, readRDS("expResults/exp1_IPI.rds"))
df <- rbind(df, readRDS("expResults/exp1_MOIMBO.rds"))
df <- rbind(df, readRDS("expResults/exp1_MultiLocalEI.rds"))
df <- rbind(df, readRDS("expResults/exp1_QEI.rds"))
df <- rbind(df, readRDS("expResults/exp1_RandomSearch.rds"))
df <- setBBOBFunToZero(df)

saveRDS(df, "expResults/exp1_total.rds")
df <- readRDS("expResults/exp1_total.rds")

createAllPlots(df, "imgResults/exp1_convergencePlots.pdf", plotType = 1)
createAllPlots(df, "imgResults/exp1_boxPlots.pdf", plotType = 2)
createAllPlots(df, "imgResults/exp1_tilePlots.pdf", plotType = 3)
createAllPlots(df, "imgResults/exp1_tilePlotsDomination.pdf", plotType = 4)
createAllPlots(df, "imgResults/exp1_bestBatchSize.pdf", plotType = 6, curve = curveDf)

fig3df <- rbind(
    filter(df, functionID == 1, batchSize == 2, nDim == 5),
    filter(df, functionID == 2, batchSize == 4, nDim == 5),
    filter(df, functionID == 15, batchSize == 8, nDim == 5)
)
createAllPlots(fig3df, "imgResults/fig3_convergencePlots.pdf", plotType = 5)

source("02_plotParallelEfficiency.R")

