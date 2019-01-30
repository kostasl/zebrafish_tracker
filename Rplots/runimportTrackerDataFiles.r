source("TrackerDataFilesImport_lib.r")

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
groupsrcdatListPerDataSet <- list()
datAllSets <-list()
n <- 0
#### List Of Data files / and result label assuming organized in Directory Structure ###
for ( idxDataSet in firstDataSet:lastDataSet )
{
  n <- n +1
  d = strDataSetDirectories[[idxDataSet]]
  groupsrcdatList = list()
  strCondR  <- "*.csv"; 
  groupsrcdatList[["LE"]] <- list(getFileSet("LiveFed/Empty/",d),"-LiveFed-Empty")
  
  groupsrcdatList[["LL"]] <- list(getFileSet("LiveFed/Live/",d),"-LiveFed-Live")
  
  groupsrcdatList[["NE"]] <- list(getFileSet("NotFed/Empty/",d),"-NotFed-Empty")
  
  groupsrcdatList[["NL"]] <- list(getFileSet("NotFed/Live/",d),"-NotFed-Live")
  
  groupsrcdatList[["DE"]] <- list(getFileSet("DryFed/Empty/",d),"-DryFed-Empty")
  
  groupsrcdatList[["DL"]] <- list(getFileSet("DryFed/Live/",d),"-DryFed-Live")
  
  
  ##OutPutFIleName
  strDataSetIdentifier <- strsplit(d,"/")
  strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
  strDataFileName <- paste(strDataExportDir,"/setn1_Dataset_", strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
  strDataFileNameRDS <- paste(strDataExportDir,"/setn1_Dataset_", strDataSetIdentifier,".rds",sep="") ##To Which To Save After Loading
  message(paste(" Importing to:",strDataFileName))
  ##RUN IMPORT FUNCTION
  datAllFrames <-importTrackerFilesToFrame(groupsrcdatList,"extractFileNameParams_huntingExp")
  datAllFrames$dataSet <- idxDataSet ##Identify DataSet
  
  datAllSets[[n]] <- datAllFrames
  
  ##CHeck If Exp Ids not found 
  stopifnot(NROW(datAllFrames[which(is.na(datAllFrames$expID)), ]) == 0)
  
  groupsrcdatListPerDataSet[[idxDataSet]] <- groupsrcdatList 
  save(datAllFrames,groupsrcdatList,file=strDataFileName) ##Save With Dataset Idx Identifier
  saveRDS(datAllFrames, file = strDataFileNameRDS)
  
  #idxDataSet = idxDataSet + 1
} ##For Each DataSet Directory
#### END OF IMPORT TRACKER DATA ############

##Save the File Sources and all The Frames Combined - Just In case there are loading Problems Of the Individual RData files from each set
save(groupsrcdatListPerDataSet,file=paste(strDataExportDir,"/groupsrcdatListPerDataSet_Ds-",firstDataSet,"-",lastDataSet,".RData",sep=""))

#datAllFrames <- rbindlist(datAllSets);
datAllFrames = do.call(rbind,datAllSets);
save(datAllFrames,file=paste(strDataExportDir,"datAllFrames_Ds-",firstDataSet,"-",lastDataSet,".RData",sep=""))
