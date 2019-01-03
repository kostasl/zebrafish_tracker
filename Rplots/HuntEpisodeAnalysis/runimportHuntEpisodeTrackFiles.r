###### ## Import THe HUNTEVENT Track Data And Merge Prey And Fish Records ### (Assummes single fish , single/multiple Prey) ##
## Kostas L May 2018
## After Tracking all events and screening Labelling the identified Hunt Events we run the
## labelling process again, now filtering for specific labels in the Hunt Event list - For example the succesfull ones
## This opens the tracker and allows supervised carefull retracking of events of interests tracking BOTH Fish And The Prey / or Preys of interest
## This Script Imports the Fish Tracking data - all Retracked Datasets CSV files are assumed  to exist in a single folder - and Joins the Food tracks with the 
## the fish tracking.
## A specicif hunt event, with fish and food motion,  can then be played back as a cartoon using the function renderHuntEventPlayback
## The list of Analysed Hunt Events Is Given by datTrackedEventsRegister
###


############# Function to Analyse HuntEpisode tracks Combined With food #########
## Import The retracked Hunt/Food Episodes - 
## Find/Record first turn to Prey  ##
## Analyse Distance To Prey Progress, Number of Bouts ... Bout length changes ##


library(signal)

source("TrackerDataFilesImport_lib.r")
source("plotTrackScatterAndDensities.r")
#################IMPORT HuntEvent TRACKER FILES # source Tracker Data Files############################### 
##OutPutFIleName
strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_SetB",".RData",sep="") ##To Which To Save After Loading
message(paste(" Importing to:",strDataFileName))


lHuntEventFOODfileSrc <- list()
lHuntEventTRACKSfileSrc <- list()


  #n <- n +1

  #lHuntEventTRACKSfileSrc[["LE"]] <- list(getFileSet("LiveFed/Empty/",strTrackeroutPath,"tracks"),"-HuntEvent-LiveFed-Empty")
  #lHuntEventFOODfileSrc[["LE"]] <- list(getFileSet("LiveFed/Empty/",strTrackeroutPath,"food"),"-HuntEventFood-LiveFed-Empty")

  lHuntEventTRACKSfileSrc[["LL"]] <- list(getFileSet("LiveFed/Live",strTrackeroutPath,"tracks"),"-HuntEvent-LiveFed-Empty")
  lHuntEventFOODfileSrc[["LL"]] <- list(getFileSet("LiveFed/Live",strTrackeroutPath,"food"),"-HuntEventFood-LiveFed-Empty")
  
  #lHuntEventTRACKSfileSrc[["DL"]] <- list(getFileSet("DryFed/Live",strTrackeroutPath,"tracks"),"-HuntEvent-DryFed-Live")
  #lHuntEventFOODfileSrc[["DL"]] <- list(getFileSet("DryFed/Live",strTrackeroutPath,"food"),"-HuntEventFood-DryFed-Live")
  
  lHuntEventTRACKSfileSrc[["NL"]] <- list(getFileSet("NotFed/Live",strTrackeroutPath,"tracks"),"-HuntEvent-NotFed-Live")
  lHuntEventFOODfileSrc[["NL"]] <- list(getFileSet("NotFed/Live",strTrackeroutPath,"food"),"-HuntEventFood-NotFed-Live")
  
  
  ##RUN IMPORT FUNCTION
  datHuntEventFrames <-importTrackerFilesToFrame(lHuntEventTRACKSfileSrc,"extractFileNameParams_huntingExp")
  ##Clean Out Non (Fish) Tracks -  Dublicates - Use Template Score TO detect Irrelevant Tracks
  datHuntEventFrames <- datHuntEventFrames[datHuntEventFrames$templateScore > 0.70,] 
  ## IMport Food Tracks And Merge With Fish Tracks In Hunt Events # Found in TrackerDataFilesImport_lib
  datHuntEventMergedFrames<- mergeFoodTrackerFilesToFrame(lHuntEventFOODfileSrc,datHuntEventFrames) ##Load Food Tracking Files And Attach Data On Respective Hunt Event Frames.
  #datHuntEventFrames$dataSet <- idxDataSet ##Identify DataSet
  
  ##CHeck If Exp Ids not found 
  stopifnot(NROW(datHuntEventFrames[which(is.na(datHuntEventFrames$expID)), ]) == 0)
  
  #Make an Updated list of ReTracked Hunt Events that have been imported
  datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID,as.character(datHuntEventMergedFrames$group) ) ))
  names(datTrackedEventsRegister) <- c("expID","eventID","trackID","groupID")
  
  save(datHuntEventMergedFrames,datTrackedEventsRegister,lHuntEventTRACKSfileSrc,lHuntEventFOODfileSrc,file=strDataFileName) ##Save With Dataset Idx Identifier
  # #### END OF IMPORT HUNT EVENT TRACKER DATA ############
   
##Save the File Sources and all The Frames Combined - Just In case there are loading Problems Of the Individual RData files from each set
  #save(lHuntEventFOODfileSrc,file=paste(strDataExportDir,"/filesrcHuntEventsFoodTracks.RData",sep=""))
  #save(lHuntEventTRACKSfileSrc,file=paste(strDataExportDir,"/filesrcHuntEventsFishTracks.RData",sep=""))
  # rle() - Use Run Length Encoding
  
#datAllFrames <- rbindlist(datAllSets);
#datAllFrames = do.call(rbind,datAllSets);
save(datHuntEventMergedFrames,file=paste(strDataExportDir,"datAllHuntEventAnalysisFrames.RData",sep=""))
