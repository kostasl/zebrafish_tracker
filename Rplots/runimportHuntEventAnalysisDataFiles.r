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

source("TrackerDataFilesImport.r")
source("plotTrackScatterAndDensities.r")
#################IMPORT HuntEvent TRACKER FILES # source Tracker Data Files############################### 
lHuntEventFOODfileSrc <- list()
lHuntEventTRACKSfileSrc <- list()

  #n <- n +1

  #lHuntEventTRACKSfileSrc[["LE"]] <- list(getFileSet("LiveFed/Empty/",strTrackeroutPath,"tracks"),"-HuntEvent-LiveFed-Empty")
  #lHuntEventFOODfileSrc[["LE"]] <- list(getFileSet("LiveFed/Empty/",strTrackeroutPath,"food"),"-HuntEventFood-LiveFed-Empty")

  lHuntEventTRACKSfileSrc[["LL"]] <- list(getFileSet("LiveFed/Live",strTrackeroutPath,"tracks"),"-HuntEvent-LiveFed-Empty")
  lHuntEventFOODfileSrc[["LL"]] <- list(getFileSet("LiveFed/Live",strTrackeroutPath,"food"),"-HuntEventFood-LiveFed-Empty")
  
  
  ##OutPutFIleName
  strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis",".RData",sep="") ##To Which To Save After Loading
  message(paste(" Importing to:",strDataFileName))
  ##RUN IMPORT FUNCTION
  datHuntEventFrames <-importTrackerFilesToFrame(lHuntEventTRACKSfileSrc)
  ##Clean Out Non (Fish) Tracks -  Dublicates - Use Template Score TO detect Irrelevant Tracks
  datHuntEventFrames <- datHuntEventFrames[datHuntEventFrames$templateScore > 0.70,] 
  ## IMport Food Tracks And Merge With Fish Tracks In Hunt Events
  datHuntEventMergedFrames<- mergeFoodTrackerFilesToFrame(lHuntEventFOODfileSrc,datHuntEventFrames) ##Load Food Tracking Files And Attach Data On Respective Hunt Event Frames.
  #datHuntEventFrames$dataSet <- idxDataSet ##Identify DataSet
  
  ##CHeck If Exp Ids not found 
  stopifnot(NROW(datHuntEventFrames[which(is.na(datHuntEventFrames$expID)), ]) == 0)
  

  datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))
  
  save(datHuntEventMergedFrames,datTrackedEventsRegister,lHuntEventTRACKSfileSrc,lHuntEventFOODfileSrc,file=strDataFileName) ##Save With Dataset Idx Identifier
  # #### END OF IMPORT HUNT EVENT TRACKER DATA ############
  
  
  ############# LOAD AND PLAYBACK ####
  load(strDataFileName)
  ##Test  PlayBack Plot Hunt Event###  
  ##Make an Updated list of ReTracked Hunt Events that have been imported
  datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))
  names(datTrackedEventsRegister) <- c("expID","eventID","trackID")
  idxH <- 10
  datRenderHuntEvent <- datHuntEventMergedFrames[datHuntEventMergedFrames$expID==datTrackedEventsRegister[idxH,]$expID 
                                                 & datHuntEventMergedFrames$trackID==datTrackedEventsRegister[idxH,]$trackID 
                                                 & datHuntEventMergedFrames$eventID==datTrackedEventsRegister[idxH,]$eventID,]
  renderHuntEventPlayback(datRenderHuntEvent,speed=1)

  #### PROCESS BOUTS ###
  vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
  vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
  vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
  #nNumberOfBouts       <- 
  dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
  dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - Why is dFrame 0?  
  dEventSpeed          <- meanf(vEventPathLength/dframe,5) ##Apply Mean Filter Smooth Out 
  
  #speed_Smoothed <- meanf(dEventSpeed,10)
  ##Replace NA with 0s
#  dEventSpeed[is.na(dEventSpeed)] = 0
  dEventSpeed_smooth <- meanf(dEventSpeed,20)
  dEventSpeed_smooth[is.na(dEventSpeed_smooth)] = 0
  MoveboutsIdx <- find_peaks(dEventSpeed_smooth*100,25)
  ##Reject Peaks Below Half An SD Peak Value - So As to Choose Only Significant Bout Movements # That Are Above the Minimum Speed to Consider As Bout
  MoveboutsIdx_cleaned <- MoveboutsIdx[which(dEventSpeed_smooth[MoveboutsIdx] > sd(dEventSpeed_smooth[MoveboutsIdx])/2 
                                             & dEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
  ##Plot Displacement and Speed(Scaled)
  X11()
  plot(cumsum(vEventPathLength)) ##PLot Total Displacemnt over time
  lines(dEventSpeed_smooth*100,type='l',col="blue")
  points(MoveboutsIdx,dEventSpeed_smooth[MoveboutsIdx]*100,col="black")
  points(MoveboutsIdx_cleaned,dEventSpeed_smooth[MoveboutsIdx_cleaned]*100,col="red")
  message(paste("Number oF Bouts:",length(MoveboutsIdx_cleaned)))

   ##Binarize , Use indicator function 1/0 for frames where Motion Occurs
   vMotionBout <- dEventSpeed_smooth
   vMotionBout[ vMotionBout < G_MIN_BOUTSPEED  ] = 0
   vMotionBout[vMotionBout > G_MIN_BOUTSPEED  ] = 1
   vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
   vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
   vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] ##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
   iPairs <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event
   ##Remove The Motion Regions Where A Peak Was not detected / Only Keep The Bouts with Peaks
   vMotionBout[1:length(vMotionBout)] = 0 ##Reset / Remove All Identified Movement
   for (i in 1:iPairs)
   {
     if (any( MoveboutsIdx_cleaned >= vMotionBout_On[i] & MoveboutsIdx_cleaned < vMotionBout_Off[i] ) == TRUE)
     { ###Motion Interval Does not belong to a detect bout(peak) so remove
       vMotionBout[vMotionBout_On[i]:vMotionBout_Off[i] ] = 1 ##Remove Motion From Vector
     }
     else
     {##Remove the Ones That Do not Have a peak In them
       vMotionBout_On[i] = NA
       vMotionBout_Off[i] = NA
     }
   }

   vMotionBoutDuration <- vMotionBout_Off[1:iPairs]-vMotionBout_On[1:iPairs]
   vMotionBoutDuration <- vMotionBoutDuration[!is.na(vMotionBoutDuration)]


   plot(vMotionBout,type='p')
   points(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red")
   points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green")
   points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="yellow")

   ######### END OF PROCESS BOUT #########


   ###########  Plot Polar Angle to Prey ##############
   X11()
   plot.new()
   polarPlotAngleToPrey(datRenderHuntEvent)
 
   X11()
   plot.new()
   polarPlotAngleToPreyVsDistance(datRenderHuntEvent)
   
   ###################################################
   
##Save the File Sources and all The Frames Combined - Just In case there are loading Problems Of the Individual RData files from each set
  #save(lHuntEventFOODfileSrc,file=paste(strDataExportDir,"/filesrcHuntEventsFoodTracks.RData",sep=""))
  #save(lHuntEventTRACKSfileSrc,file=paste(strDataExportDir,"/filesrcHuntEventsFishTracks.RData",sep=""))
  
#datAllFrames <- rbindlist(datAllSets);
#datAllFrames = do.call(rbind,datAllSets);
save(datHuntEventMergedFrames,file=paste(strDataExportDir,"datAllHuntEventAnalysisFrames.RData",sep=""))
