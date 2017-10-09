runSWMM<- function(inpFile,rptFile,outFile,SWMMexe='swmm5.exe',verbose=T){
  #runs swmm on the inpFile, rptFile, outFile provided.  Does not create directories
  inpFile <- paste('"',inpFile,'"',sep = "")
  rptFile <- paste('"', rptFile,'"',sep="")
  outFile <- paste('"', outFile,'"',sep="")
  
  command <- paste(SWMMexe,inpFile,rptFile,outFile, sep = " ")
  
  #run the initial simulation
  system(command)
  flush.console()
}
checkSWMMForErrors<-function(outFile){
  #Checks a SWMM output file for errors
  f=file(outFile,"rb")
  seek(f,-6*4,"end")
  output={}
  output$position.objectID=readBin(f, integer(), n = 1, size = 4)
  output$position.objectProperties=readBin(f, integer(), n = 1, size = 4)
  output$position.computedResults=readBin(f, integer(), n = 1, size = 4)
  output$numReportingPeriods=readBin(f, integer(), n = 1, size = 4)
  output$errorStatus=readBin(f, integer(), n = 1, size = 4)
  return(output$errorStatus)
}
openSWMMOut<-function(outFile,verbose=T){
  #gets the header of a binary output file
  RECORDSIZE=4
  output={}
  f=file(outFile, "rb")
  
  output$header=readBin(f, integer(), n = 7, size = 4)
  output$numSubc=output$header[4]
  output$numNode=output$header[5]
  output$numLink=output$header[6]
  output$numPoll=output$header[7]
  output$unitCode=output$header[3]
  seek(f,-6*4,"end")
  output$position.objectID=readBin(f, integer(), n = 1, size = 4)
  output$position.objectProperties=readBin(f, integer(), n = 1, size = 4)
  output$position.computedResults=readBin(f, integer(), n = 1, size = 4)
  output$numReportingPeriods=readBin(f, integer(), n = 1, size = 4)
  output$errorStatus=readBin(f, integer(), n = 1, size = 4)
  if(output$errorStatus>0){
    print(paste("SWMM Errored out with code",output$errorStatus))
    return(output)
  }
  #Get Object IDs
  seek(f,output$position.objectID,"start");
  #For all subcatchments
  output$subcNames={}
  if(output$numSubc>0){
    for(i in 1:output$numSubc){
      lengthName=readBin(f, integer(), n = 1, size = 4)
      output$subcNames[i]=readChar(f, lengthName, useBytes = FALSE)
    }
  }
  #For all nodes
  output$nodeNames={};
  if(output$numNode>0){
    for(i in 1:output$numNode){
      lengthName=readBin(f, integer(), n = 1, size = 4)
      output$nodeNames[i]=readChar(f, lengthName, useBytes = FALSE)
    }
  }
  #For all links
  output$linkNames={};
  if(output$numLink>0){
    for(i in 1:output$numLink){
      lengthName=readBin(f, integer(), n = 1, size = 4)
      #       print(lengthName)
      output$linkNames[i]=readChar(f, lengthName, useBytes = FALSE)
    }
  }
  #For all pollutants
  output$pollNames={};
  output$pollUnits={};
  if(output$numPoll>0){
    #A bug was fixed in this section on 1/10/2012 at 1:34 pm
    for(i in 1:output$numPoll){
      
      
      lengthName=readBin(f, integer(), n = 1, size = 4)
      output$pollNames[i]=readChar(f, lengthName, useBytes = FALSE)
      
    }
    for(i in 1:output$numPoll){
      
      
      unitCode=readBin(f, integer(), n = 1, size = 4)
      #print(unitCode)
      if(unitCode==0){
        output$pollUnits[i]="mg/L"
      }else if(unitCode==1){
        output$pollUnits[i]="ug/L"
      }else if(unitCode==2){
        output$pollUnits[i]="counts/L"
      }
    }
    #end 1/10/2012 bug fix
  }
  seek(f,output$position.objectProperties,"start")
  
  
  #Subcatchments
  
  output$numSubcPropSaved= readBin(f, integer(), n = 1, size = 4)
  output$codesSubcPropSaved=readBin(f, integer(), n = 1, size = 4)
  if (output$codesSubcPropSaved==1){
    output$subcArea=readBin(f,what="double",n=output$numSubc,size=4);
    #print(output$subcArea)
  }
  
  #Nodes
  
  output$numNodePropSaved=readBin(f, integer(), n = 1, size = 4)
  output$codesNodePropSaved=readBin(f, integer(),n=output$numNodePropSaved,size=4);
  
  if(output$numNode>0){
    temp=readBin(f,what="double",n=output$numNodePropSaved*output$numNode,size=4)
    codestemp=temp[seq(from=1,by=3,to=length(temp))]
    output$nodeType={}
    count=0
    for( i in codestemp){
      count=count+1
      if(i==0){
        output$nodeType[count]="Junction"
      }else if(i==1){
        output$nodeType[count]="Outfall"
      }else if(i==2){
        output$nodeType[count]="Storage"
      }else if(i==3){
        output$nodeType[count]="Divider"
      }
    }
    output$nodeInvert=temp[seq(from=2,by=3,to=length(temp))]
    output$nodeMaxDepth=temp[seq(from=3,by=3,to=length(temp))]
  }
  
  
  #Links
  output$numLinkPropSaved=readBin(f, integer(), n = 1, size = 4)
  output$codesLinkPropSaved=readBin(f, integer(), n =output$numLinkPropSaved, size = 4)
  
  
  if(output$numLink>0){
    
    temp=readBin(f,what="double",n=output$numLinkPropSaved*output$numLink,size=4)
    codestemp=temp[seq(from=1,by=5,to=length(temp))]
    output$linkType={}
    count=0
    for( i in codestemp){
      count=count+1
      if(i==0){
        output$linkType[count]="Conduit"
      }else if(i==1){
        output$linkType[count]="Pump"
      }else if(i==2){
        output$linkType[count]="Orifice"
      }else if(i==3){
        output$linkType[count]="Weir"
      }else if(i==4){
        output$linkType[count]="Outlet"
      }
    }
    
    output$linkUpstreamInvertOffset=temp[seq(from=2,by=5,to=length(temp))]
    output$linkDownstreamInvertOffset=temp[seq(from=3,by=5,to=length(temp))]
    output$linkMaxDepth=temp[seq(from=4,by=5,to=length(temp))]
    output$linkLength=temp[seq(from=5,by=5,to=length(temp))]
  }
  
  
  output$outFileHandle=f
  output$numSubcVars=readBin(f, integer(), n = 1, size = 4)
  
  output$subcVarCodes=readBin(f,integer(),n=  output$numSubcVars,size=4)
  output$numNodeVars=readBin(f, integer(), n = 1, size = 4)
  
  output$nodeVarCodes=readBin(f,integer(), n=  output$numNodeVars,size=4)
  output$numLinkVars=readBin(f, integer(), n = 1, size = 4)
  output$linkVarCodes=readBin(f,integer(), n=  output$numLinkVars,size=4)
  output$numSysVars=readBin(f,integer(),n=1,size=4)
  output$sysVarCodes==readBin(f,integer(), n=output$numSysVars,size=4)
  
  
  
  output$bytesPerPeriod= 2*RECORDSIZE +
    (output$numSubc*(output$numSubcVars) +
       output$numNode*(output$numNodeVars) +
       output$numLink*(output$numLinkVars) +
       output$numSysVars)*RECORDSIZE;
  
  
  return(output)
}
getSWMMResult<-function(headObj,iType,iIndex,vIndex,period){
  #gets a single SWMM result.  See the documentation in SWMM intefacing for detail
  #iType is either 0,1,2,3 for subcatch, node, link or sys variable
  #vIndex is the index of the variable for subcatch, node, link or sys object
  #iIndex is the position of the subcatch, node, or link among the other subcatch, links, nodes
  #I would recommend using getSWMMTimeSeries so that you don't have to know iIndex.  That function is a wrapper
  #around this one and looks up results based on names not indicies.
  
  #1/24/2012 edit: Fixed bug in returning results for models with pollutants
  RECORDSIZE=4
  SUBCATCH=0
  NODE     = 1;
  LINK     = 2;
  SYS      = 3;
  
  f=headObj$outFileHandle
  StartPos=headObj$position.computedResults
  off = StartPos + period*(headObj$bytesPerPeriod) + 2*RECORDSIZE;
  if ( iType == SUBCATCH )
  {
    off = off+ RECORDSIZE*(iIndex*(headObj$numSubcVars) + vIndex);
  }
  else if (iType == NODE)
  {
    off = off+ RECORDSIZE*(headObj$numSubc*(headObj$numSubcVars) +
                             iIndex*(headObj$numNodeVars) + vIndex);
  }
  else if (iType == LINK)
  {
    off = off+ RECORDSIZE*(headObj$numSubc*(headObj$numSubcVars) +
                             headObj$numNode*(headObj$numNodeVars) +
                             iIndex*(headObj$numLinkVars) + vIndex);
  }
  else if (iType == SYS)
  {
    off = off+ RECORDSIZE*(headObj$numSubc*(headObj$numSubcVars) +
                             headObj$numNode*(headObj$numNodeVars) +
                             headObj$numLink*(headObj$numLinkVars) + vIndex);
    
  }
  
  
  seek(f,off,"start")
  output=readBin(f,what="double",size=4,n=1)
  return(output)
}
getSWMMTimes<-function(headObj){
  #gets the time stamps of the SWMM results in binary file
  f=headObj$outFileHandle
  seek(f,headObj$position.computedResults,"start")
  
  headObj$SWMMTimes<-array(NaN,headObj$numReportingPeriods)
  if(headObj$numReportingPeriods>0){
    for(i in 1:headObj$numReportingPeriods){
      
      headObj$SWMMTimes[i]<-readBin(f,what="double",size=8,n=1)
      #      if(i<100){
      #       print(headObj$SWMMTimes[i])
      #    }
      seek(f,headObj$bytesPerPeriod-8,"current")
    }
  }else{
    stop("No time steps listed in SWMM output file.")
    
  }
  #Convert SWMM times to R POSIXlt datetimes
  headObj$SWMMTimes<-headObj$SWMMTimes*86400.0+as.POSIXct(strptime("12/30/1899", format="%m/%d/%Y",tz="MST"))#edit 2/10/2012 to force GMT time zone rather than locale specific
  return(headObj)
}
getSWMMTimeSeriesData<-function(headObj,iType,nameInOutputFile,vIndex){
  #headObj should be an object obtained by calling openSWMMOut
  #iType should be 0 for Subcatchments
  #                1 for nodes
  #                2 for links
  #                3 for system variables
  #nameInOutputFile should be the exact name in the output file of a subcatchment,
  # link, or node, or leave this as an empty string if searching for system results
  #vIndex should be selected from this lists below,depending on whether iType is a subcatchment,
  # link or node.
  ############################################################
  #######BEGIN vIndex choices##################################
  ############################################################
  #Number of subcatchment variables (currently 6 + number of pollutants).
  # Code number of each subcatchment variable:
  # 0 for rainfall (in/hr or mm/hr),
  # 1 for snow depth (in or mm),
  # 2 for evaporation + infiltration losses (in/hr or mm/hr),
  # 3 for runoff rate (flow units),
  # 4 for groundwater outflow rate (flow units),
  # 5 for groundwater water table elevation (ft or m),
  # 6 for runoff concentration of first pollutant,
  # 5 + N for runoff concentration of N-th pollutant.
  #
  #Number of node variables (currently 6 + number of pollutants)
  #Code number of each node variable:
  # 0 for depth of water above invert (ft or m),
  # 1 for hydraulic head (ft or m),
  # 2 for volume of stored + ponded water (ft3 or m3),
  # 3 for lateral inflow (flow units),
  # 4 for total inflow (lateral + upstream) (flow units),
  # 5 for flow lost to flooding (flow units),
  # 6 for concentration of first pollutant,
  # 5 + N for concentration of N-th pollutant.
  #
  #Number of link variables (currently 5 + number of pollutants)
  #Code number of each link variable:
  #
  # 0 for flow rate (flow units),
  # 1 for flow depth (ft or m),
  # 2 for flow velocity (ft/s or m/s),
  # 3 for Froude number,
  # 4 for capacity (fraction of conduit filled),
  # 5 for concentration of first pollutant,
  # 4 + N for concentration of N-th pollutant.
  #
  #Number of system-wide variables (currently 14)
  #Code number of each system-wide variable:
  # 0 for air temperature (deg. F or deg. C),
  # 1 for rainfall (in/hr or mm/hr),
  # 2 for snow depth (in or mm),
  # 3 for evaporation + infiltration loss rate (in/hr or mm/hr),
  # 4 for runoff flow (flow units),
  # 5 for dry weather inflow (flow units),
  # 6 for groundwater inflow (flow units),
  # 7 for RDII inflow (flow units),
  # 8 for user supplied direct inflow (flow units),
  # 9 for total lateral inflow (sum of variables 4 to 8) (flow units),
  # 10 for flow lost to flooding (flow units),
  # 11 for flow leaving through outfalls (flow units),
  # 12 for volume of stored water (ft3 or m3),
  # 13 for evaporation rate (in/day or mm/day)
  #
  
  if(iType==0){
    iIndex=(0:(headObj$numSubc-1))[headObj$subcNames==nameInOutputFile]
  }else if(iType==1){
    iIndex=(0:(headObj$numNode-1))[headObj$nodeNames==nameInOutputFile]  
  }else if(iType==2){
    iIndex=(0:(headObj$numLink-1))[headObj$linkNames==nameInOutputFile]
  }else if(iType==3){
    
    iIndex=0
  }
  output=array(NaN,headObj$numReportingPeriods)
  for(period in 0:(-1+headObj$numReportingPeriods)){
    #browser()
    output[period+1]=getSWMMResult(headObj=headObj,iType=iType,iIndex=iIndex,vIndex=vIndex,period=period)
  }
  return(output)
}
getCalDataFromCSV<-function(CSVFile,dateFormat="%m/%d/%y %H:%M"){
  #Calibration Data should be in a CSV with datetimes in the first column,
  #and data in the second column
  #The text file is assumed to have a one line header
  #Call this function with the correct dateFormat for your datetimes
  #the dateFormat is passed to strptime, so look for formatting information there
  #for example, dates like this 1/1/07 12:00, can be read with the default dateFormat
  
  temp=read.csv(file=CSVFile, header = TRUE, sep = ",", quote="\"", dec=".",
                fill = TRUE, comment.char="",stringsAsFactors = FALSE)
  cData<-{}
  cData$times<-as.POSIXct(strptime(temp[,1], format=dateFormat,tz="MST")) #edit 2/10/2012 to force GMT time zone rather than locale specific
  cData$obs<-temp[,2]
  return(cData)
  
}
interpCalDataToSWMMTimes<-function(headObj, cData){
  #this linearly interpolates observations for calibration onto the time stamps of the SWMM model
  cData$interpedObs<-approx(cData$times, cData$obs, headObj$SWMMTimes,
                               method="linear",
                               yleft=NaN, yright=NaN, rule = 1)$y
  return(cData)
}
readTemplateFile<-function(SWMMTemplateFile){
  #reads in the SWMM file that has replacement codes
  SWMMTemplate=readLines(con = SWMMTemplateFile, n = -1L, ok = TRUE, warn = TRUE,
                         encoding = "unknown")
  return(SWMMTemplate)
}
replaceCodesInTemplateFile<-function(SWMMTemplate,parameters,replacementCodes){
  #replaces the codes in an input file with parameters from optimization
  for(i in 1:length(parameters)){
    
    SWMMTemplate=gsub(replacementCodes[i], parameters[i], SWMMTemplate,fixed=TRUE)
  }
  return(SWMMTemplate)
}
writeNewInputFile<-function(SWMMTemplate,filename){
  #writes a new input file after replacement codes have been replaced by parameters
  SWMMTemplate=gsub("//", "////", SWMMTemplate, ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)
  writeLines(SWMMTemplate, con = filename, sep = "\n", useBytes = FALSE)
  
}
performanceStatsAsMinimization<-function(correspondingSWMMSeries,correspondingSWMMSeriesTSS,correspondingSWMMSeriesTP,headObj, includeOnlyNonNAN=TRUE, cData, waterQual=FALSE, cData.TSS=NULL, cData.TP=NULL, k){
  #calDataObj should have the fields times, obs, and interpedObs, which can be
  #obtained by calling getCalDataFromCSV and interpCalDataToSWMMTimes
  #headObj should be a header returned by openSWMMOut and modified by getSWMMTimes
  #correspondingSWMMSeries should be a series of just data points for comparison to
  #calDataObj$interpedObs.  It could be obtained by calling getSWMMTimeSeriesData.  It
  #should have the same length as calDataObj$interpedObs
  #All performance stats given are expressed so that they can be minimized to maximize the fit
  #of the model to the data.  For example, correlation is mulitplied by -1,
  #so that minimizing it improves correlation
  #2/21/2012 edit: added includeOnlyNonNAN argument to restrict the calculation to only non-NAN values
  if(waterQual){
    pos.ind <- which(cData.TSS$obs>0)
    cData.TSS$times <- cData.TSS$times[pos.ind]
    cData.TSS$obs <- cData.TSS$obs[pos.ind]
    pos.ind <- which(cData.TP$obs>0)
    cData.TP$times <- cData.TP$times[pos.ind]
    cData.TP$obs <- cData.TP$obs[pos.ind]
    
    tss.ind <- which((as.double(cData.TSS$times)>=as.double(headObj$SWMMTimes[1]))
                     &(as.double(cData.TSS$times)<=as.double(headObj$SWMMTimes[length(headObj$SWMMTimes)])))
    res.tss <- c()
    swmmTimes <- round.POSIXt(headObj$SWMMTimes,units = "secs")
    cd.tss <- round.POSIXt(cData.TSS$times[tss.ind], units = "secs")
    
    for (i in 1:length(tss.ind)){
      ind <- which(swmmTimes==cd.tss[i])
      res.tss[i] <- correspondingSWMMSeriesTSS[ind]-as.numeric(cData.TSS$obs[tss.ind[i]])
    }
    rootMeanSquaredErrorTSS=sqrt(mean(res.tss^2, na.rm = T))
    tp.ind <- which((as.double(cData.TP$times)>=as.double(headObj$SWMMTimes[1]))
                     &(as.double(cData.TP$times)<=as.double(headObj$SWMMTimes[length(headObj$SWMMTimes)])))
    res.tp <- c()
    cd.tp <- round.POSIXt(cData.TP$times[tp.ind], units = "secs")
    for (i in 1:length(tp.ind)){
      ind <- which(swmmTimes==cd.tp[i])
      res.tp[i] <- correspondingSWMMSeriesTP[ind]-as.numeric(cData.TP$obs[tp.ind[i]])
    }
    rootMeanSquaredErrorTP=sqrt(mean(res.tp^2, na.rm = T))
  }
  
  residualz=cData$interpedObs-correspondingSWMMSeries
  if(includeOnlyNonNAN){
    residualz=residualz[!is.nan(residualz)]
  }
  meanAbsoluteError=mean(residualz)
  
  sumOfSquaredError=sum(residualz^2)
  
  rootMeanSquaredError=sqrt(mean((residualz)^2))
  
  
  
  library(caTools)
  f.ind <-1:216615
  s.ind <-234031:245130
  if (siteID[k]==2){
    obs.vol <- trapz(headObj$SWMMTimes[f.ind],cData$interpedObs[f.ind])+trapz(headObj$SWMMTimes[s.ind],cData$interpedObs[s.ind])
    sim.vol <- trapz(headObj$SWMMTimes[f.ind],correspondingSWMMSeries[f.ind])+trapz(headObj$SWMMTimes[s.ind],correspondingSWMMSeries[s.ind])
    absoluteVolumeDifference=abs(obs.vol-sim.vol)/(obs.vol)*100
  }else{
    absoluteVolumeDifference=abs(trapz(headObj$SWMMTimes,cData$interpedObs)
                               -trapz(headObj$SWMMTimes,correspondingSWMMSeries))/(trapz(headObj$SWMMTimes,cData$interpedObs))*100
  }

    
  
  linearCorrelationTimesMinus1=-1*cor(cData$interpedObs[!is.nan(cData$interpedObs)],correspondingSWMMSeries[!is.nan(cData$interpedObs)])
  mean4thPowerError=mean(residualz^4)
  if (waterQual){
    output=data.frame(meanAbsoluteError,sumOfSquaredError,linearCorrelationTimesMinus1,mean4thPowerError,rootMeanSquaredError,absoluteVolumeDifference, rootMeanSquaredErrorTSS,rootMeanSquaredErrorTP)
  }else{
    output=data.frame(meanAbsoluteError,sumOfSquaredError,linearCorrelationTimesMinus1,mean4thPowerError,rootMeanSquaredError,absoluteVolumeDifference)
  }
  #  browser()
  return(output)
  
}
openOptimizationHistoryFile<-function(historyFilename){
  #optimization history is written to a CSV as it goes
  optimizationHistoryFile=file(historyFilename,"w")
  return(optimizationHistoryFile)
}

getParameterBounds<-function(parameterBoundsFile){
  #parameterBoundsFile should be a full filename and path of a CSV with 4 columns
  # and one header row.  The format is as follows (but do not put the comment marker # in there):
  #Code,Minimum,Maximum,Initial
  # $1$,0,1,.5
  # $2$,1,2,1.5
  #The codes are the dollar sign enclosed codes you used in the template file to denote uncertain parameters
  #Min and Max are the ranges of the parameter search for each parameter
  #Initial is the best guess used to initialize the optimization
  parametersTable=read.csv(file=parameterBoundsFile, header = TRUE, sep = ",", quote="\"", dec=".",
                           fill = TRUE, comment.char="#",stringsAsFactors = FALSE)
  return(parametersTable)
}
objectiveFunction<-function(x,
                            performanceStat,
                            functionCallToEvalForASWMMTimeSeries,
                            functionCallToEvalForASWMMTimeSeriesTSS,
                            functionCallToEvalForASWMMTimeSeriesTP,
                            baseOutputName,
                            SWMMTemplateFile,
                            SWMMexe, numSub, 
                            oWeights=NULL, iteration=1){
  # A number of global variables must be defined before this function is called
  #iteration should be set to 1 in the global namespace
  #parametersTable should be created in the global namespace by calling getParameterBounds
  #SWMMTemplate should be created in the global namespace
  #calDataObj should exist in the global namespace after calling calDataObj=readTemplateFile(...)
  #optimizationHistoryFile should exist in the global namespace after calling
  #optimizationHistoryFile=openOptimizationHistoryFile(filename) on a given filename
  #baseFilename should be something like this C:\\SWMM Models\\output
  #the iteration number and file extensions .inp, .rpt, and .out will be added to it

  if (is.null(oWeights)){
    oWeights <- rep(1,length(performanceStat))
  }else if (length(oWeights)!=length(performanceStat)){
    stop(paste("SWMM returned an error. oWeights are not the same 
               length as the performanceStat."))
  }else{
    oWeights <- oWeights
  }
  
  widthResults <- c()
  timpResults <- c()
  wResults <- c()
  fl <- SWMMTemplateFile
  s <- readLines(fl)
  ws <- rm_white(s)
  b <- strsplit(ws, " ")
  sb <- c()
  taResults <- c()
  impResults <- c()
  for (k in 1:numSub){
    sb[k] <- b[[k-1+subStart]][1]
    sel <- subData[subData$FID==as.numeric(paste(strsplit(sb[k], split = "")[[1]][-1],collapse = "")),]
    taResults[k] <- sel$TA
    impResults[k] <- sel$Imp
    wResults[k] <- (taResults[k]*10000)/(2*sel$MC)
  }
  
  
  ind.kw <- grep("kWidth",parametersTable$Name)
  ind.imp <- grep("PCTImp",parametersTable$Name)
  
  widthResults <- x[ind.kw[1]]/100*wResults+wResults
  timpResults <- x[ind.imp[1]]/100*impResults+impResults
  widthResults[32] <- x[ind.kw[2]]/100*wResults[32]+wResults[32]
  timpResults[32] <- x[ind.imp[2]]/100*impResults[32]+impResults[32]
  widthResults[20] <- x[ind.kw[3]]/100*wResults[20]+wResults[20]
  timpResults[20] <- x[ind.imp[3]]/100*impResults[20]+impResults[20]
  widthResults[15] <- x[ind.kw[4]]/100*wResults[15]+wResults[15]
  timpResults[15] <- x[ind.imp[4]]/100*impResults[15]+impResults[15]
  widthResults[7] <- x[ind.kw[5]]/100*wResults[7]+wResults[7]
  timpResults[7] <- x[ind.imp[5]]/100*impResults[7]+impResults[7]
  
  SWMMTemplate <- readTemplateFile(SWMMTemplateFile)
  
  SWMMTemplateModified=replaceCodesInTemplateFile(SWMMTemplate,x,as.matrix(parametersTable["Code"]))
  for (i in 1:numSub){
    SWMMTemplateModified <- replaceCodesInTemplateFile(SWMMTemplateModified,widthResults[i],replacementCodes = paste("$S",paste(strsplit(sb[i], split = "")[[1]][-1],collapse = ""),"$", collapse = "",sep = "")) 
    SWMMTemplateModified <- replaceCodesInTemplateFile(SWMMTemplateModified,timpResults[i], replacementCodes = paste("$IMP",paste(strsplit(sb[i],split = "")[[1]][-1],collapse = ""),"$",collapse = "",sep = ""))
  }
  
  inpFile=paste(baseOutputName,iteration,'.inp',sep="")
  rptFile=paste(baseOutputName,iteration,'.rpt',sep="")
  outFile=paste(baseOutputName,iteration,'.out',sep="")
  writeNewInputFile(SWMMTemplateModified,inpFile)
  runSWMM(inpFile,rptFile,outFile,SWMMexe,verbose=T)
  if(!file.exists(outFile)){errCode=1}else{errCode=checkSWMMForErrors(outFile)}
  
  if(errCode==0){
    

    ps <- c()
    for (k in 1:length(nodeID)){
      headObj<-openSWMMOut(outFile)#removed the returnObjectProperties argment on 2/7/2011
      headObj<-getSWMMTimes(headObj)
      calData <- getCalDataFromCSV(flow.calCSV[k],dateFormat="%Y-%m-%d %H:%M:%S")
      if (siteID[k]==2){calData$obs <- calData$obs*0.0283168}
      calData <- interpCalDataToSWMMTimes(headObj, calData)
      calData.TSS <- getCalDataFromCSV(tss.calCSV[k],dateFormat="%m/%d/%Y %H:%M")
      calData.TP <- getCalDataFromCSV(tp.calCSV[k],dateFormat="%m/%d/%Y %H:%M")
      
      correspondingSWMMSeries <- eval(parse(text=functionCallToEvalForASWMMTimeSeries))
      ind <- which(headObj$SWMMTimes>=calData$times[1]&headObj$SWMMTimes<=calData$times[length(calData$times)])
      headObj$SWMMTimes <- headObj$SWMMTimes[ind]
      calData$interpedObs <- calData$interpedObs[ind]
      correspondingSWMMSeries <- correspondingSWMMSeries[ind]
      correspondingSWMMSeriesTSS <- eval(parse(text=functionCallToEvalForASWMMTimeSeriesTSS))
      correspondingSWMMSeriesTP <- eval(parse(text=functionCallToEvalForASWMMTimeSeriesTP))
      correspondingSWMMSeriesTSS <- correspondingSWMMSeriesTSS[ind]
      correspondingSWMMSeriesTP <- correspondingSWMMSeriesTP[ind]
      #     print(correspondingSWMMSeries[1:20])
      perfStats=performanceStatsAsMinimization(correspondingSWMMSeries, 
                                               correspondingSWMMSeriesTSS, 
                                               correspondingSWMMSeriesTP, headObj, 
                                               cData = calData,
                                               includeOnlyNonNAN = TRUE,
                                               waterQual=TRUE, 
                                               cData.TSS = calData.TSS, 
                                               cData.TP = calData.TP, k = k)
      perfStatsToUse <-as.numeric(perfStats[performanceStat])
      perfStatsToUse <- perfStatsToUse*oWeights
      ps <- c(ps,perfStatsToUse)
    }
    
    


    # browser()
    close(headObj$outFileHandle)
  }else if(errCode>0){
    summaryRow=unlist(c(iteration,x,perfStatsToUse))
    
    
    stop(paste("SWMM returned an error. Optimization stopping.  See", rptFile))
    return(array(NaN,length(performanceStat)))
    
  }

  #browser()

  return(ps) ;
}

readLID<-function(filename,headObj){
  library(zoo)
  #readLID function added 1/11/2012 after version 1 of RSWMM
  #you need the zoo library to run this: install.packages("zoo")
  #headObj should be obtained by calling both openSWMMOut and getSWMMTimes
  #filename is the LID report filename
  #skips the first 9 header lines and reads an LID file
  #   Elapsed      Total      Total	  Surface	     Soil	   Bottom	  Surface	    Drain	  Surface	    Soil/	  Storage
  #    Time	   Inflow	     Evap	    Infil	     Perc	    Infil	   Runoff	  Outflow	    Depth	    Pave 	    Depth
  #   Hours	    in/hr	    in/hr	    in/hr	    in/hr	    in/hr	    in/hr	    in/hr	   inches	    Moist	   inches
  # -------	 --
  
  r=read.table(filename,skip=9,sep='\t')
  names(r)<-c("ElapsedHours","TotalInflow","TotalEvap","SurfaceInfil","SoilPerc","BottomInfil","SurfaceRunoff","DrainOutflow","SurfaceDepth","SoilMoist","StorageDepth")
  pp=as.numeric(headObj$SWMMTimes-headObj$SWMMTimes[1])/3600.0
  sz=zoo(rep(0,length(headObj$SWMMTimes)),pp)
  rz=zoo(r[,2:dim(r)[2]],r[,1])
  rz=rz[,1:dim(rz)[2]-1]
  tz=merge(sz,rz,fill=0.0)
  first.of.timeStep <- function(tt) floor(tt/(pp[2]-pp[1]))*(pp[2]-pp[1])
  
  # average z over quarters
  # 1. via "yearqtr" index (regular)
  # 2. via "Date" index (not regular)
  #z.qtr1 <- aggregate(z.day, as.yearqtr, mean)
  #z.qtr2 <- aggregate(z.day, first.of.quarter, mean)
  tz=aggregate(tz,first.of.timeStep,mean)
  return(tz)
}






