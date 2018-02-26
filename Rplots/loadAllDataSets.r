##########  LOAD DATA SETS #####
##de


datAllSets <- list()
groupsrcdatListPerDataSet <- list() ##Holds File List Per Data Set - Used to Cross Ref FileIdx 
n <- 0
for (i in dataSetsToProcess )
{
  n <- n + 1
  strDataSetIdentifier <- strsplit(strDataSetDirectories[[i]],"/")
  strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
  strDataFileName <- paste("setn1_Dataset_", strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
  
  #strDataFileName <- paste("setn1_Dataset_",strsplit(strDataSetDirectories[[i]],"/")[[1]][[2]],".RData",sep="") ##To Which To Save After Loading
  
  message(paste("...Loading ",strDataFileName) )
  load(strDataFileName,verbose=TRUE)
  
  groupsrcdatListPerDataSet[[i]] <- groupsrcdatList 
  stopifnot(datAllFrames$dataSet ==  i) ##Identify DataSet
  ##Add Only Sequencially, so there are no Empty Entries in the List
  datAllSets[[n]] <- datAllFrames
}
message(paste("Done Loading All datasets. Now Merging... ",dataSetsToProcess[1]," to ",dataSetsToProcess[length(dataSetsToProcess)]) )

#### Need to Filter Out Empty Entries in the List before binding otherwise R session crashes with mem. violation (Issue #2340)
#datAllFrames = do.call(rbind,datAllSets);
#datAllFrames <- rbindlist(lapply(datAllSets,alloc.col) ) 
datAllFrames <- rbindlist(datAllSets)
#datAllFrames <- mapply(c,datAllSets)
##Done LOADING Required DataSets

source("HuntingEventAnalysis.r")

##### LOAD Processed Dat Hunt Stat FROM FILE ######
strCondTags <- list("LE","LL","NE","NL","DE","DL")
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
datHuntStat <- loadHuntEvents(strCondTags,dataSetsToProcess)

### END OF HUNT STAT LOAD ####
