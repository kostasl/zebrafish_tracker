######################################   START DATA Processing SCRIPT   ###############################
#######################################################################################
source("HuntingEventAnalysis_lib.r")
source("TrajectoryAnalysis.r")

lTrackletStat <- list();
lHuntStat     <- list();
lHuntEvents   <- list();
lMotionStat   <- list();


### TRAJECTORIES Indicating Hunting  - With distinct colour for each larva ####
#source("plotTrackScatterAndDensities.r")
##########

strRGroupsTags <- unique(datAllFrames$group) ##Rearing Group names(groupsrcdatList)
strTCondTags <- unique(datAllFrames$testCond) ##Test Condition names(groupsrcdatList)
#### Process Files  #####
### Calculates Statistics/Makes Plots 
## Extracts Hunting events and saves them on output csv files / one for each conditiongroup - showing start-end frames and file
i = 0
for (g in strRGroupsTags)
{
  for (c in strTCondTags)
  {
    i = i +1
    message(paste("#### Process Rearing Group ",g," Tested in",c, "###############"))
    subsetDat = groupsrcdatList[[c]]
    #strCond   <- paste(c,subsetDat[2],collapse=NULL);
    
    ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
    #vexpID = unique(filtereddatAllFrames$expID)
    ##Select Larvaof this Group
    
    datAllGroupFrames <- datAllFrames[which(datAllFrames$group == g &  datAllFrames$testCond == c) ,]
    if (NROW(datAllGroupFrames) < 1)
    {
      warning("No Frames in Group, skipping")
      next
    }
    #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
    vexpID = unique(datAllGroupFrames$expID)
    idxDataSet <- unique(datAllGroupFrames$dataSet)
    
    datHuntEvent = detectHuntEvents(datAllGroupFrames,vexpID,c,dataSetsToProcess)
    lHuntEvents[[i]] = datHuntEvent ## Collect TO One Data File
    writeHuntEventToFile(datHuntEvent,dataSetsToProcess,groupsrcdatListPerDataSet)
    
    #lTrackletStat[[i]] <- calcRecordingEventSpeed(datAllGroupFrames,vexpID,idxDataSet)
  
    ## DISABLED as these are not being used Anymore ##  
    ## Extract Hunting Events From Data
    #lMotionStat[[i]] <- calcMotionStat(datAllGroupFrames,vexpID,dataSetsToProcess)
  
    
    lHuntStat[[g]] <- calcHuntStat3(datHuntEvent)
    #lHuntStat[[i]] = calcHuntStat3(datHuntEvent)
    stopifnot(length(lHuntStat[[g]]$vHLarvaEventCount) > 0)
  
    ##Combine Hunting Events across fish in this Condition In One
    #datHuntEvent = do.call(rbind,lHuntStat[[i]]$vHuntingEventsList )
    ##Reconstruct DataSet File List - So As to link fileIdx To Files
    #filelist <- getFileSet("LiveFed/Empty/",strDataSetDirectories[[idxDataSet]])
  }
    
} ##Process Stat For Each Condition / Write Hunting Events
message("~~ Finished Going Through Groups and Test Conditions ~~ ")

datAllHuntEvents = do.call(rbind,lHuntEvents)

## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)#
datMotionStat = do.call(rbind,lMotionStat)
datTrackletStat = data.frame(do.call(rbind,lTrackletStat))
#datHuntStat <- rbindlist(lHuntStat)
#datMotionStat <-rbindlist(lMotionStat)
#Data Exported In One Dir -> strDataExportDir, and read from another - so as to Avoid accidental Overwrites
saveRDS(datAllHuntEvents,file=paste(strDataExportDir,"/HB_allHuntEvents.rds",sep="" )) ##So It can be loaded into a custom Named Structure

save(datTrackletStat,lTrackletStat,file =paste(strDataExportDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datTrackletStat.RData",sep="")) 

save(datHuntStat, file=paste(strDataExportDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datHuntStat.RData",sep=""))
save(datMotionStat, file=paste(strDataExportDir,"/","setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datMotionStat.RData",sep=""))

## Track Lengths ##
#datTrackletStat_filt <- datTrackletStat[datTrackletStat$Length_mm < 200,]
#boxplot(unlist(Length_mm) ~ unlist(groupID),data=datTrackletStat_filt, notch = TRUE)
# 
# ### Examine Activity Between Groups  ##
# vPathLengthPerLarva_LE <- tapply(unlist(lTrackletStat[["LE"]]$Length_mm),unlist(lTrackletStat[["LE"]]$expID),sum)
# vPathLengthPerLarva_DE <- tapply(unlist(lTrackletStat[["DE"]]$Length_mm),unlist(lTrackletStat[["DE"]]$expID),sum)
# vPathLengthPerLarva_NE <- tapply(unlist(lTrackletStat[["NE"]]$Length_mm),unlist(lTrackletStat[["NE"]]$expID),sum)
# vPathLengthPerLarva_LL <- tapply(unlist(lTrackletStat[["LL"]]$Length_mm),unlist(lTrackletStat[["LL"]]$expID),sum)
# vPathLengthPerLarva_DL <- tapply(unlist(lTrackletStat[["DL"]]$Length_mm),unlist(lTrackletStat[["DL"]]$expID),sum)
# vPathLengthPerLarva_NL <- tapply(unlist(lTrackletStat[["NL"]]$Length_mm),unlist(lTrackletStat[["NL"]]$expID),sum)
# 
# ## Empty Test - Path Length per larva
# pdf(file= paste(strPlotExportPath,"/boxplot_totalPathLengthPerGroup.pdf",sep=""),onefile=TRUE ) 
# boxplot(vPathLengthPerLarva_LE,vPathLengthPerLarva_DE,vPathLengthPerLarva_NE,
#         vPathLengthPerLarva_LL,vPathLengthPerLarva_DL,vPathLengthPerLarva_NL,
#         ylim=c(0,10000),notch = TRUE,main="Total Path Length Per Larva  ",
#         ylab="mm",xlab="",
#         show.names=TRUE,names=c("LE","DE","NE","LL","DL","NL"),col=c(colourR[2],colourR[1],colourR[3] ) )
# dev.off()
# 
# 
# 
# ## Live Test - Path Length per larva
# pdf(file= paste(strPlotExportPath,"/boxplotEM_totalPathLengthPerGroup_Live.pdf",sep=""),onefile=TRUE ) 
# boxplot(vPathLengthPerLarva_LL,vPathLengthPerLarva_DL,vPathLengthPerLarva_NL,
#         ylim=c(0,1000),notch = FALSE,main="Total Path Length Per Larva  ",
#         ylab="mm",xlab="",
#         show.names=TRUE,names=c("LL","DL","NL"),col=c(colourR[2],colourR[1],colourR[3] ) )
# dev.off()

