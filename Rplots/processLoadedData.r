######################################   START DATA Processing SCRIPT   ###############################
#######################################################################################
source("HuntingEventAnalysis_lib.r")
source("TrajectoryAnalysis.r")

lTrackletStat <- list();
lHuntStat     <- list();
lMotionStat   <- list();


### TRAJECTORIES Indicating Hunting  - With distinct colour for each larva ####
#source("plotTrackScatterAndDensities.r")
##########

strCondTags <- names(groupsrcdatList)
#### Process Files  #####
### Calculates Statistics/Makes Plots 
## Extracts Hunting events and saves them on output csv files / one for each conditiongroup - showing start-end frames and file
for (i in strCondTags)
{
  message(paste("#### ProcessGroup ",i," ###############"))
  subsetDat = groupsrcdatList[[i]]
  strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);
  
  ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
  #vexpID = unique(filtereddatAllFrames$expID)
  ##Select Larvaof this Group
  
  datAllGroupFrames <- datAllFrames[which(datAllFrames$group == i),]
  #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
  vexpID = unique(datAllGroupFrames$expID)
  idxDataSet <- unique(datAllGroupFrames$dataSet)
  
  #  lHuntStat[[i]] = calcHuntStat(datAllGroupFrames,vexpID)
  ##Combine Hunting Events across fish in this Condition In One
  #datHuntEvent = do.call(rbind,lHuntStat[[i]]$vHuntingEventsList )
  
  lTrackletStat[[i]] <- calcRecordingEventSpeed(datAllGroupFrames,vexpID,idxDataSet)
  
  ## Extract Hunting Events From Data
  #lMotionStat[[i]] <- calcMotionStat(datAllGroupFrames,vexpID,dataSetsToProcess)

  #datHuntEvent = detectHuntEvents(datAllGroupFrames,vexpID,dataSetsToProcess)
  #writeHuntEventToFile(datHuntEvent,dataSetsToProcess,groupsrcdatListPerDataSet)
  
  #lHuntStat[[i]] <- calcHuntStat3(datHuntEvent)
  #stopifnot(length(lHuntStat[[i]]$vHLarvaEventCount) > 0)
  
  ##Reconstruct DataSet File List - So As to link fileIdx To Files
  #filelist <- getFileSet("LiveFed/Empty/",strDataSetDirectories[[idxDataSet]])
  #filelist <-  
  
  
} ##Process Stat For Each Condition / Write Hunting Events


## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)#
datMotionStat = do.call(rbind,lMotionStat)

#datHuntStat <- rbindlist(lHuntStat)
#datMotionStat <-rbindlist(lMotionStat)
#Data Exported In One Dir -> strDataExportDir, and read from another - so as to Avoid accidental Overwrites
save(datHuntStat, file=paste(strDataExportDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datHuntStat.RData",sep=""))
save(datMotionStat, file=paste(strDataExportDir,"/","setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datMotionStat.RData",sep=""))


