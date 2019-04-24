###  Kostasl 2018 
##  Validate the capture strike data that was automatically detected from the retracked huntevents :
##  * after running runHuntEpisodeAnalysis.r, a bout detection list for each of the hunt episodes is used to obtain what the undershoot 
##  turn was on firstbout, and link it to the final, capture bout intensity.


source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
## The Original list if the lFirstBout data from runHuntepisode analysis
datMotionBoutsToValidate <-readRDS(file=paste0(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds") ) 

#datHuntEventAllGroupToValidate <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
##Save a backup before changing anything
saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate_backup.rds",sep="")) ##Save With Dataset Idx Identifier
saveRDS(datTrackedEventsRegister, paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate_backup.rds",sep="") ) ## THis is the Processed Register File On 

##Make desired columns if missing 
if (!any(names(datMotionBoutsToValidate) == "MarkValidated"))
  datMotionBoutsToValidate$MarkValidated <- NA
if (!any(names(datMotionBoutsToValidate) == "vd_PreyX"))
  datMotionBoutsToValidate$vd_PreyX <- NA
if (!any(names(datMotionBoutsToValidate) == "vd_PreyY"))
  datMotionBoutsToValidate$vd_PreyY <- NA
if (!any(names(datMotionBoutsToValidate) == "vd_MouthX"))
  datMotionBoutsToValidate$vd_MouthX <- NA
if (!any(names(datMotionBoutsToValidate) == "vd_MouthY"))
  datMotionBoutsToValidate$vd_MouthY <- NA



##Get the Capture strike bout subset
datCaptureBoutsToValidate <- datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ] 

## Exclude the already Validated records
idxToValidate <- datCaptureBoutsToValidate[is.na(datCaptureBoutsToValidate$MarkValidated), ]$RegistarIdx

#idxToValidate <- c(idxReg_DL,idxReg_NL,idxReg_LL)

for (idx in idxToValidate)
{
  rec <- datTrackedEventsRegister[idx,]
  strFormatedEventID <- gsub(" ","0",formatC(format="fg",width=3,as.character(rec$eventID)))
  
  ##Get Respective Video  Filename / Path Uniquely identified by the  expID and event ID
  strVideoFile <- list.files(path =strVideoFilePath, pattern = paste0("_",rec$expID,"_",strFormatedEventID,".mp4") , all.files = FALSE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)  
  if (NROW(strVideoFile) == 0)
    stop(paste("Could not find video file matching,",rec$filenames, " ,in : ",strVideoFilePath ) )
  if (!file.exists(strVideoFile) )
    stop(paste("Video File ",strVideoFile, "does not exist" ) )
  
  recBout <- datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ]
  
  message(paste("\n", row.names(rec) ,". Examining Hunt Event -start:",max(0,rec$startFrame + recBout$vMotionBout_On-1)," -End:",rec$endFrame + recBout$vMotionBout_Off, "ExpID:",rec$expID ) )
  strArgs = paste(" --MeasureMode=1 --HideDataSource=0 --ModelBG=0 --SkipTracked=0 --PolygonROI=1 --invideofile=",strVideoFile," --outputdir=",strTrackeroutPath,"/BoutValidate/",
                  " --startframe=",max(0,rec$startFrame + recBout$vMotionBout_On-1)," --stopframe=",max(0,rec$endFrame + recBout$vMotionBout_Off-1)," --startpaused=1",sep="")
  message(paste(strTrackerPath,"/zebraprey_track",strArgs,sep=""))
  if (!file.exists(paste(strTrackerPath,"/zebraprey_track",sep="")) )
    stop(paste("Tracker software not found in :",strTrackerPath ))
  
  execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout=NULL,stderr =NULL) ## stdout=FALSE stderr = FALSE
  
  if (execres != 0)
    stop(execres) ##Stop If Application Exit Status is not success
  
  flush(con=stdin())
  flush(con=stdout())
  ##Show user Prompt for bout validation input
  message(paste("\n\n ### Event's ", row.names(rec) , "  ####" ) )
  
  strDat <-""
  while ( nchar(strDat) < 3 && substr(strDat,1,1) != 'c')
    strDat <- readline(prompt="# Validation Data [Prey X, Y, Mouth X,Y,Frame Bout On, Bout Off ] :")
  
  
  ## extract validation data on prey location and fish mouth position ##
  strDat <- gsub(strDat,pattern = "[][]",replacement = "")
  arrDat <- strsplit(strDat,"," )[[1]]
  stopifnot(NROW(arrDat) > 5 ) ##Stop if array data not found
  preyX <- as.numeric( arrDat[1]);  preyY  <- as.numeric(arrDat[2]);
  mouthX <- as.numeric(arrDat[3]); mouthY <- as.numeric(arrDat[4]);
  boutOn <- as.numeric(arrDat[5]) - rec$startFrame; boutOff <- as.numeric(arrDat[6]) - rec$startFrame
  fOnSetEyeVergence <- as.numeric(arrDat[7])
  ##update temp record - which push to main data.frame after user mark-validates it
  ## \note Angle to prey Not Recalculated , changes to bout frames  invalidates vMotionBoutDistanceTravelled_mm
  recUpd<- within( datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ],{ ##Multiple vars update
    vd_PreyY <- preyY 
    vd_PreyX <- preyX
    vd_MouthX <- mouthX
    vd_MouthY <- mouthY
    vMotionBout_On <- boutOn
    vMotionBout_Off <- boutOff 
    vMotionBoutDistanceToPrey_mm <- DIM_MMPERPX*(sqrt((preyX-mouthX)^2+(preyY-mouthY)^2 ))
    OnSetEyeVergence <- fOnSetEyeVergence
  })
  
  message(paste("Distance to prey:", prettyNum(digits=3,recUpd$vMotionBoutDistanceToPrey_mm)," eyeV:",prettyNum(digits=3,recUpd$OnSetEyeVergence) ))
  ## Pass data to record   // SAVE
  strKeyC <- readline(prompt="### Mark Validated ? (y/n):")
  if (strKeyC == 'y')
  {
    recUpd$MarkValidated <- 1
    ##update data frame
    datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ] <- recUpd
    
    ##Update the register end frame - is this was a capture (last bout) 
    if (recUpd$boutRank == 1)
      datTrackedEventsRegister[idx,]$endFrame <- rec$startFrame + recUpd$vMotionBout_Off
    
    ## print event label 
    message("Outcome logged on Registry (matched) :",convertToScoreLabel(rec$LabelledScore ))
    print(levels(convertToScoreLabel(rec$LabelledScore )))
    strKeyC <- readline(prompt=paste("### Change label to [",convertToScoreLabel(rec$LabelledScore ),"]:") )
    
    if (is.numeric(strKeyC))
    {
      newLab <- convertToScoreLabel(as.numeric(strKeyC)-1 ) 
      message(paste("new label:",newLab ) )
      datTrackedEventsRegister[idx,]$LabelledScore <- newLab        
    }
      
  }
  
  strKeyC <- readline(prompt="### Move to next or quit ? (n/q):")
  if (strKeyC == 'q')
    break;
  
  
} ## end of loop 

##Update Capture bouts back to original dataframe ##
datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ]  <- datCaptureBoutsToValidate

saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds",sep="")) ##Save With Dataset Idx Identifier
saveRDS(datTrackedEventsRegister, paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 


### Recalc derived bout data / capture speed and angle to prey ###
##datHuntEventMergedFrames
load(file=paste(strDataExportDir,"datAllHuntEventAnalysisFrames_setC.RData",sep=""))
source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r")



idxRegValidated <- datMotionBoutsToValidate[!is.na(datMotionBoutsToValidate$MarkValidated) & datMotionBoutsToValidate$MarkValidated == 1,]$RegistarIdx

for (idx in idxRegValidated )
{
  recReg <- datTrackedEventsRegister[idx,]
  #### PROCESS Validated BOUTS ###
  datPlaybackHuntEvent <- datHuntEventMergedFrames[datHuntEventMergedFrames$expID==recReg$expID 
                                                   & datHuntEventMergedFrames$trackID==recReg$trackID 
                                                   & datHuntEventMergedFrames$eventID==recReg$eventID,]
  
  ## Find the motion bouts of this event that have been validated (so we  can recalc their derived data ie speed)
  datMotionBouts <- datMotionBoutsToValidate[!is.na(datMotionBoutsToValidate$MarkValidated)&
                                               datMotionBoutsToValidate$MarkValidated == 1 &
                                               datMotionBoutsToValidate$RegistarIdx == row.names( recReg), ]
  
  vDeltaDisplacement   <- sqrt(diff(datPlaybackHuntEvent$posX,lag=1,differences=1)^2+diff(datPlaybackHuntEvent$posY,lag=1,differences=1)^2) ## Path Length Calculated As Total Displacement
  
  #nNumberOfBouts       <- 
  dframe               <- diff(datPlaybackHuntEvent$frameN,lag=1,differences=1)
  dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - and Convert To Time sec  
  vFs                   <- datPlaybackHuntEvent$fps[1:NROW(dframe)]
  vEventSpeed_smooth          <- meanf(vDeltaDisplacement[1:NROW(dframe)]/dframe,5) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed_smooth) #meanf(vEventSpeed,100) #
  vEventSpeed_smooth[vEventSpeed_smooth < 0] <- 0 ## Remove -Ve Values As an artefact of Filtering
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventSpeed_smooth_mm <- vFs*vEventSpeed_smooth*DIM_MMPERPX
  
  ## Get peak Capture Speed within capture bout ##
  for (idxBout in row.names(datMotionBouts) )
  {
    oldVal <- datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm
    datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm <- max( vEventSpeed_smooth_mm[datMotionBoutsToValidate[idxBout,]$vMotionBout_On:datMotionBoutsToValidate[idxBout,]$vMotionBout_Off],na.rm=TRUE)
    stopifnot(is.numeric(datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm) | is.infinite((datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm)))
    print(paste(idxBout,"speed",oldVal," new:",datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm))
  }
  ## Get Angle To Prey at onset and offset 
  
  
  ## Do Next Validated rec ##
}


## Recalc First Bout Data based on Validated Info ###
strGroupID <- levels(datTrackedEventsRegister$groupID)
lFirstBoutPoints <- list() ##Add Dataframes Of 1st bout Turns for Each Group
###### PLOT BOUTTURN Vs Prey Angle Coloured with BOUTSEQ ################
for (gp in strGroupID)
{
  groupID <- which(levels(datTrackedEventsRegister$groupID) == gp)
  
  datMotionBoutsToValidate$vMotionBoutDistanceToPrey_mm <- as.numeric(datMotionBoutsToValidate$vMotionBoutDistanceToPrey_mm)
  datMotionBoutCombined <-datMotionBoutsToValidate[datMotionBoutsToValidate$groupID == as.numeric(groupID), ] #Select Group
  
  datMotionBoutCombined$boutRank <- as.numeric(datMotionBoutCombined$boutRank)
  datMotionBoutTurnToPrey <- datMotionBoutCombined[abs(datMotionBoutCombined$OnSetAngleToPrey) >= abs(datMotionBoutCombined$OffSetAngleToPrey) , ]
  datMotionBoutTurnToPrey <- datMotionBoutTurnToPrey[!is.na(datMotionBoutTurnToPrey$RegistarIdx),]
  
  ## Relates First turn to prey to final capture strike parameters  ##
  ##
  lFirstBoutPoints[[gp]] <- cbind(OnSetAngleToPrey = datMotionBoutTurnToPrey[datMotionBoutTurnToPrey$turnSeq == 1 ,]$OnSetAngleToPrey,
                                  Turn= datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$OnSetAngleToPrey - datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1,]$OffSetAngleToPrey
                                  , RegistarIdx=datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx,
                                  CaptureSpeed = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                          datMotionBoutCombined$boutRank == 1 ,]$vMotionPeakSpeed_mm,
                                  DistanceToPrey = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                            datMotionBoutCombined$boutRank == 1 ,]$vMotionBoutDistanceToPrey_mm,
                                  doesCaptureStrike=( datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                               datMotionBoutCombined$boutRank == 1 ,]$vMotionPeakSpeed_mm >= G_THRES_CAPTURE_SPEED ),
                                  CaptureStrikeFrame=( datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                                datMotionBoutCombined$boutRank == 1 ,]$vMotionBout_On  ),
                                  CaptureStrikeEyeVergence = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                                      datMotionBoutCombined$boutRank == 1 ,]$OnSetEyeVergence,
                                  Validated = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                       datMotionBoutCombined$boutRank == 1 ,]$MarkValidated
  )
  
}

##Save List on First Bout Data
saveRDS(lFirstBoutPoints,file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="") ) #Processed Registry on which we add )
