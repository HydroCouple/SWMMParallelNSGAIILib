###################################################################################################
#
#                               START OF Calibration PROGRAM 
#
###################################################################################################

#import the required libraries
library(qdapRegex)
library(parallel)
library(doParallel)
library(foreach)

#clear out any residual variable values
rm(list=ls())

#specify the path to the recompiled SWMM executable and the SWMM input template file
SWMMexe <- '/uufs/chpc.utah.edu/common/home/u6013151/Downloads/swmm51012_engine_2/source5_1_012/swmm5'
SWMMTemplateFile <- '/uufs/chpc.utah.edu/common/home/u6013151/SWMM/NewApproach_Continuous/Models/NWFC/Continuous_NWFC_2016.inp'

#source the program that has many of the SWMM processing functions
source("/uufs/chpc.utah.edu/common/home/u6013151/R/RSWMM/SWMMFunctions_Continuous_ParallelTest.R")

#load the modified "nsga2" and "apply" functions. This is so I could parallelize the 
#algorithm.
load("/uufs/chpc.utah.edu/common/home/u6013151/R/RSWMM/newFunctions.RData", envir = .GlobalEnv)

#initialize
numSub <- 60 #number of subcatchments in the model
#find the location in the SWMM input file where the subcatchments are listed
model <- strsplit(as.character(readLines(SWMMTemplateFile)), split = "\\s+")
ind.1 <- which(sapply(model, "[",1)=="[SUBCATCHMENTS]")[1]
subStart <- ind.1+3
#create a vector of the siteIDs. this program only considers 1800 North which is 
#siteID=2
siteID <<- c(2)
#variables that tell the program which values to get from the SWMM simulation. 
#more details are in the SWMMFunctions_Continuous_ParallelTest.R file
iType <- c(1)
vIndex <- c(4)
vIndex.tss <- c(6)
vIndex.tp <- c(7)

#get the file that define the bounds of each parameter
paramBoundsFile <- '/uufs/chpc.utah.edu/common/home/u6013151/SWMM/NewApproach_Continuous/Calibration Files/parameterBounds.txt'
#initialize the 
flow.calCSV <<- c()
nodeID <<- c()
tss.calCSV <<- c()
tp.calCSV <<- c()
tdp.calCSV <<- c()
i <- 1
if (length(which(siteID==2))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_1800N_2016.csv"
  nodeID[i] <- "FID2689"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_1800N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_1800N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_1800N.csv"
  i <- i+1
}
if (length(which(siteID==3))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_300N.csv"
  nodeID[i] <- "FID4202"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_300N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_300N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_300N.csv"
  i <- i+1
}
if (length(which(siteID==4))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_1250N.csv"
  nodeID[i] <- "FID3364"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_1250N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_1250N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_1250N.csv"
  
  i <- i+1
}
if (length(which(siteID==5))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_800N.csv"
  nodeID[i] <- "FID4372"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_800N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_800N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_800N.csv"
  
  i <- i+1
}
if (length(which(siteID==6))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_1300N.csv"
  nodeID[i] <- "FID4498"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_1300N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_1300N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_1300N.csv"
  
  i <- i+1
}
if (length(which(siteID==7))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_1000N.csv"
  nodeID[i] <- "FID4806"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_1000N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_1000N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_1000N.csv"
  
  i <- i+1
}
if (length(which(siteID==8))>0){
  flow.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/FlowCal_1400N.csv"
  nodeID[i] <- "FID5865"
  tss.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tssCal_1400N.csv"
  tp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tpCal_1400N.csv"
  tdp.calCSV[i] <- "/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/tdpCal_1400N.csv"
  
  i <- i+1
}

#get the geospatial characteristics of the subcatchments. This comes from output from
#ArcGIS processes.
subData <- read.table("/uufs/chpc.utah.edu/common/home/u6013151/SWMM/CalibrationFiles/Subcatchments.txt", sep = " ", header = F)

colnames(subData) <- c("x","y","FID","TA","Slope","Street","Res","Com","Imp","Build","RG2015","RG2016_Spr", "RG2016_Fal", "MC")

#get the parameter ranges from the text file
parametersTable <- getParameterBounds(paramBoundsFile)


#Initialize the mco list and give it attributes for each of the inputs to the nsga2 function
#that will be called later.
mcoOpt <- {}
mcoOpt$functionCallToEvalForASWMMTimeSeries <- 'getSWMMTimeSeriesData(headObj,iType = iType[k], nodeID[k], vIndex = vIndex[k])'
mcoOpt$functionCallToEvalForASWMMTimeSeriesTSS <- 'getSWMMTimeSeriesData(headObj,iType = iType[k], nodeID[k], vIndex = vIndex.tss[k])'
mcoOpt$functionCallToEvalForASWMMTimeSeriesTP <- 'getSWMMTimeSeriesData(headObj,iType = iType[k], nodeID[k], vIndex = vIndex.tp[k])'
mcoOpt$lower <- c(as.vector(parametersTable["Minimum"]))$Minimum
mcoOpt$upper <- c(as.vector(parametersTable["Maximum"]))$Maximum
mcoOpt$SWMMexe <- SWMMexe
mcoOpt$performanceStats <- c("rootMeanSquaredError","absoluteVolumeDifference", "rootMeanSquaredErrorTSS", "rootMeanSquaredErrorTP")
mcoOpt$baseOutputName <- "/scratch/general/lustre/u6013151/Temp/NWFC_"
mcoOpt$SWMMTemplateFile <- SWMMTemplateFile
mcoOpt$numGenerations <-100 
mcoOpt$popsize <- 100
#mcoOpt$oWeights <- c(10000,1,100)
optimizationHistory <- data.frame()

#create blank input, output, and report SWMM files
for (i in 1:mcoOpt$popsize){
  x.out <- paste(mcoOpt$baseOutputName,toString(i), ".out",sep = "")
  x.inp <- paste(mcoOpt$baseOutputName,toString(i), ".inp",sep = "")
  x.rpt <- paste(mcoOpt$baseOutputName,toString(i), ".rpt",sep = "")
  file.create(x.out)
  file.create(x.inp)
  file.create(x.rpt)
}

#detect the number of cores and the current machine
nCores <- detectCores()
nCores
#create a cluster that uses those cores
cl <- makeCluster(nCores, type = "FORK")

#export all variable values to that cluster
clusterExport(cl, list(ls()))
#register the cluster
registerDoParallel(cl)
#initialize a counter that gets the number of generations, objective function values for 
#each generation, and the model parameters for each value. they are then written to an output
#file
numGen <<- 0
opResults <<- data.frame()
paramResults <<- data.frame()

#load the mco library
library(mco)
#run the nsga2 algorithm and store the results in the "out" variable. 
out <- nsga2(objectiveFunction, 
             idim=length(mcoOpt$lower), 
             odim = length(nodeID)*length(mcoOpt$performanceStats),
             baseOutputName=mcoOpt$baseOutputName,
             SWMMTemplateFile=mcoOpt$SWMMTemplateFile,
             SWMMexe=mcoOpt$SWMMexe,
             numSub=numSub,
             functionCallToEvalForASWMMTimeSeries=mcoOpt$functionCallToEvalForASWMMTimeSeries,
             functionCallToEvalForASWMMTimeSeriesTSS=mcoOpt$functionCallToEvalForASWMMTimeSeriesTSS,
             functionCallToEvalForASWMMTimeSeriesTP=mcoOpt$functionCallToEvalForASWMMTimeSeriesTP,
             performanceStat=mcoOpt$performanceStats,
             oWeights = mcoOpt$oWeights,
             generations = mcoOpt$numGenerations,
             popsize = mcoOpt$popsize,
             lower.bounds = mcoOpt$lower,
             upper.bounds = mcoOpt$upper,
             constraints = NULL)

#write all the output files
colnames(opResults) <- c("generations",rep(mcoOpt$performanceStats, length(nodeID)))
colnames(paramResults) <- c("generations", parametersTable$Name)
write.csv(out$par, file = "/scratch/general/lustre/u6013151/Results/param.csv")
write.csv(out$value, file = "/scratch/general/lustre/u6013151/Results/val.csv")
write.csv(out$pareto.optimal, file = "/scratch/general/lustre/u6013151/Results/pareto.csv")
write.csv(mcoOpt$oWeights, file = "/scratch/general/lustre/u6013151/Results/oWeights.csv")
write.csv(opResults, file = "/scratch/general/lustre/u6013151/Results/valHistory.csv")
write.csv(paramResults, file = "/scratch/general/lustre/u6013151/Results/paramHistory.csv")
stopCluster(cl)
closeAllConnections()


