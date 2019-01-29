##########  LOAD DATA SETS #####
## Note : It takes care of differences in Field countsbetween dataFrames/DataSets, by using the names in the previous frame to take them as subset from the next Dataset
## 


datAllSets <- list()
groupsrcdatListPerDataSet <- list() ##Holds File List Per Data Set - Used to Cross Ref FileIdx 
n <- 0

fieldNames <- "expID"
for (i in dataSetsToProcess )
{
  n <- n + 1
  strDataSetIdentifier <- strsplit(strDataSetDirectories[[i]],"/")
  strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
  strDataFileName <- paste("setn1_Dataset_", strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
  
  #strDataFileName <- paste("setn1_Dataset_",strsplit(strDataSetDirectories[[i]],"/")[[1]][[2]],".RData",sep="") ##To Which To Save After Loading
  
  message(paste("...Loading ",strDataFileName) )
  load(paste(strDatDir,strDataFileName,sep="/"),verbose=TRUE)
  
  #datAllFrames<-readRDS(strDataFileName)
  
  groupsrcdatListPerDataSet[[i]] <- groupsrcdatList 
  stopifnot(datAllFrames$dataSet ==  i) ##Identify DataSet
  ##Add Only Sequencially, so there are no Empty Entries in the List
  if (n > 1) ##Only Combine DataSets on the fields existing On the Prior dataFrame
    datAllSets[[n]] <- datAllFrames[fieldNames]
  else
    datAllSets[[n]] <- datAllFrames
  
  if (n==1)
    fieldNames <- names(datAllFrames)
  
  
}
message(paste("Done Loading All datasets. Now Merging... ",dataSetsToProcess[1]," to ",dataSetsToProcess[length(dataSetsToProcess)]) )

#### Need to Filter Out Empty Entries in the List before binding otherwise R session crashes with mem. violation (Issue #2340)
datAllFrames = do.call(rbind,datAllSets);


## Save the aggregated DAta file ##

strDataSetIdentifierFirst <- strsplit(strDataSetDirectories[[head(dataSetsToProcess,1)]],"/")
strDataSetIdentifierFirst <- strDataSetIdentifierFirst[[1]][[ length(strDataSetIdentifierFirst[[1]]) ]]
strDataSetIdentifier <- strsplit(strDataSetDirectories[[tail(dataSetsToProcess,1)]],"/")
strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
strDataFileName <- paste("setn",NROW(dataSetsToProcess),"_Datasets_", strDataSetIdentifierFirst ,"-",strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
save(datAllFrames,file=paste(strDatDir,strDataFileName,sep="/"))


message(paste("Saved Merged files to:",strDatDir,strDataFileName) )

#datAllFrames <- rbindlist(lapply(datAllSets,alloc.col) ) 
#datAllFrames <- rbindlist(datAllSets)
#datAllFrames <- mapply(c,datAllSets)
##Done LOADING Required DataSets

source("HuntingEventAnalysis_lib.r")

##### LOAD Processed Dat Hunt Stat FROM FILE ######
#strCondTags <- list("LE","LL","NE","NL","DE","DL")
#dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
#datHuntStat <- loadHuntEvents(strCondTags,dataSetsToProcess)

### END OF HUNT STAT LOAD ####
