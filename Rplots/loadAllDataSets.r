##########  LOAD DATA SETS #####

datAllSets <- list()
groupsrcdatListPerDataSet <- list() ##Holds File List Per Data Set - Used to Cross Ref FileIdx 
for (i in dataSetsToProcess )
{
  strDataSetIdentifier <- strsplit(strDataSetDirectories[[i]],"/")
  strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
  strDataFileName <- paste("setn1_Dataset_", strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
  
  #strDataFileName <- paste("setn1_Dataset_",strsplit(strDataSetDirectories[[i]],"/")[[1]][[2]],".RData",sep="") ##To Which To Save After Loading
  
  message(paste("...Loading ",strDataFileName) )
  load(strDataFileName,verbose=TRUE)
  
  groupsrcdatListPerDataSet[[i]] <- groupsrcdatList 
  stopifnot(datAllFrames$dataSet ==  i) ##Identify DataSet
  datAllSets[[i]] <- datAllFrames
}

datAllFrames = do.call(rbind,datAllSets);
##Done LOADING Required DataSets


##### LOAD Processed Dat Hunt Stat FROM FILE ######
strCondTags <- list("LE","LL","NE","NL","DE","DL")
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
datHuntStat <- loadHuntEvents(strCondTags,dataSetsToProcess)

### END OF HUNT STAT LOAD ####
