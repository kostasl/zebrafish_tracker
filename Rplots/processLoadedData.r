######################################   START DATA Processing SCRIPT   ###############################
#######################################################################################
source("HuntingEventAnalysis.r")
source("TrajectoryAnalysis.r")

lHuntStat <- list();
lMotionStat <- list();


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
  
  ## Extract Hunting Events From Data
  datHuntEvent = detectHuntEvents(datAllGroupFrames,vexpID,dataSetsToProcess)
  lMotionStat[[i]] <- calcMotionStat(datAllGroupFrames,vexpID,dataSetsToProcess)
  
  lHuntStat[[i]] <- calcHuntStat3(datHuntEvent)
  
  stopifnot(length(lHuntStat[[i]]$vHLarvaEventCount) > 0)
  ##Reconstruct DataSet File List - So As to link fileIdx To Files
  #filelist <- getFileSet("LiveFed/Empty/",strDataSetDirectories[[idxDataSet]])
  #filelist <-  
  
  
  if (NROW(datHuntEvent) > 0 )
  {
    message("Writing Hunting Data...")
    ## Recover The File names from fileIdxs
    ## Run It for Each DataSet Idx - Not Sure Of Better Way extracting the relevant filelist
    datHuntEvent$filenames = "." ##Create Field
    for (d in idxDataSet)
    {
      
      ##Get Files Used for This DataSet, and this Condition
      filelist <- getVideofilePath(unlist(groupsrcdatListPerDataSet[d][[1]][[i]][[1]]),strVideoFilePath)
      
      ##Set File Name
      
      datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0,]$filenames <- filelist[ datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0,]$fileIdx ]
    }
    
    strDataFileName <- paste("setn",NROW(dataSetsToProcess),"-D",dataSetsToProcess[1],"-",dataSetsToProcess[NROW(dataSetsToProcess)],"-HuntEvents-",i,sep="") ##To Which To Save After Loading
    
    write.csv(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".csv",sep="" ) , row.names=FALSE ) 
    ###Save Hunt Event Data Frame
    
    message(paste(" Exporting to:",strDataFileName))
    ##ExPORT 
    datHuntEvent$groupID = i
    save(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    saveRDS(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".rds",sep="" )) ##So It can be loaded into a custom Named Structure
    
  }else{
    message("No Hunting Event to write!")
  }
  
  
} ##Process Stat For Each Condition / Write Hunting Events


## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)#
datMotionStat = do.call(rbind,lMotionStat)

#datHuntStat <- rbindlist(lHuntStat)
#datMotionStat <-rbindlist(lMotionStat)
#Data Exported In One Dir -> strDataExportDir, and read from another - so as to Avoid accidental Overwrites
save(datHuntStat, file=paste(strDataExportDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datHuntStat.RData",sep=""))
save(datMotionStat, file=paste(strDataExportDir,"/","setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datMotionStat.RData",sep=""))



### Loads Each HuntEvent RData File Found in A given directory and merges the records, which can then be identified by groupID
mergeHuntEventRecords <- function(strSrcDir,strExt = "*.RData")
{

  ldatHunt <- list()
#  groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
#  strCondTags <- names(groupsrcdatList)
  ldatFiles <- getFileSet("",strSrcDir,strExt)
  
  i <- 1
  for (f in ldatFiles)
  {
    message(paste("#### Load HuntEvent File ",f," ###############"))
    load(file=paste(f,"",sep="" )) ##Save With Dataset Idx Identifier
    
    #datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13),labels=vHuntEventLabels )##Set To NoTHuntMode

    ldatHunt[[i]] <- datHuntEvent
    i <- i +1
  }

  datHuntEvents = do.call(rbind,ldatHunt)#

  nDat <- NROW(unique(datHuntEvents$dataSetID) )
  first <- levels(datHuntEvents$dataSetID)[min(unique(as.numeric(datHuntEvents$dataSetID)) )]
  last <- levels(datHuntEvents$dataSetID)[max(unique(as.numeric(datHuntEvents$dataSetID) ) )]

  saveRDS(datHuntEvents,file=paste(strDataExportDir,"/setn",nDat,"-D",first,"-",last,"-HuntEvents-Merged.rds",sep="")  )
}
