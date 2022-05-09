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
  groupsrcdatList[["LE"]] <- list((getFileSet("LE/",d,strCondR= "tracks_[0-9].csv")),"-TestEmpty")
  groupsrcdatList[["LR"]] <- list((getFileSet("LR/",d,strCondR= "*tracks_[0-9].csv")),"-TestPrey")
  
  ##OutPutFIleName
  strDataSetIdentifier <- strsplit(d,"/")
  strDataSetIdentifier <- strDataSetIdentifier[[1]][[ length(strDataSetIdentifier[[1]]) ]]
  strDataFileName <- paste(strDataExportDir,"/setn1_Dataset_", strDataSetIdentifier,".RData",sep="") ##To Which To Save After Loading
  strDataFileNameRDS <- paste(strDataExportDir,"/setn1_Dataset_", strDataSetIdentifier,".rds",sep="") ##To Which To Save After Loading
  message(paste(" Importing to:",strDataFileName))
  ##RUN IMPORT FUNCTION
  #datAllFrames <- importTrackerFilesToFrame(groupsrcdatList,"extractFileNameParams_HungerExp_camB") ##"extractFileNameParams_huntingExp")
  datAllFrames <- importTrackerFilesToFrame(groupsrcdatList,"extractFileNameParams_OliviaAssay") ##"extractFileNameParams_huntingExp")
  
  datAllFrames$dataSet <- idxDataSet ##Identify DataSet
  
  datAllSets[[n]] <- datAllFrames
  message("Imported frames for GroupID:",unique(datAllFrames$groupID) )
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
