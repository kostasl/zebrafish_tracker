source("TrackerDataFilesImport.r")

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 


#### List Of Data files / and result label assuming organized in Directory Structure ###
for ( idxDataSet in firstDataSet:lastDataSet )
{
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
  strDataFileName <- paste("setn1_Dataset_",strsplit(d,"/")[[1]][[2]],".RData",sep="") ##To Which To Save After Loading
  message(paste(" Importing to:",strDataFileName))
  ##RUN IMPORT FUNCTION
  datAllFrames <-importTrackerFilesToFrame(groupsrcdatList)
  datAllFrames$dataSet <- idxDataSet ##Identify DataSet
  save(datAllFrames,groupsrcdatList,file=strDataFileName) ##Save With Dataset Idx Identifier
  
  #idxDataSet = idxDataSet + 1
} ##For Each DataSet Directory
#### END OF IMPORT TRACKER DATA ############
