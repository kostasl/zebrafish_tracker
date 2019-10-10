##'  Kostas Lagogiannis 
##' April 2019 
##' Validate the capture strike data that was automatically detected from the retracked huntevents :
##' * after running runHuntEpisodeAnalysis.r, a bout detection list for each of the hunt episodes is used to obtain what the undershoot 
##' turn was on firstbout, and link it to the final, capture bout intensity.
##' these data are then used in statistical models to establish relationships between capture speed, distance to prey, eye vergence and undershoot on 1st turn to prey


source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
## The Original list if the lFirstBout data from runHuntepisode analysis
datMotionBoutsToValidate <-readRDS(file=paste0(strDataExportDir,"huntEpisodeAnalysis_MotionBoutData_peakSpeedExtended_Validated.rds"))  #"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds") ) 

#datHuntEventAllGroupToValidate <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
##Save a backup before changing anything
saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate_backup_fixed.rds",sep="")) ##Save With Dataset Idx Identifier
saveRDS(datTrackedEventsRegister, paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate_backup_fixed.rds",sep="") ) ## THis is the Processed Register File On 

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

#datTrackedEventsRegister$LabelledScore <- convertToScoreLabel( datTrackedEventsRegister$LabelledScore)


##Get the Capture strike bout subset
datCaptureBoutsToValidate <- datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ] 

## Exclude the already Validated records
idxToValidate <- datCaptureBoutsToValidate[is.na(datCaptureBoutsToValidate$MarkValidated), ]$RegistarIdx

#idxToValidate <- {194}
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
  arrDat <- ""
  while ( NROW(arrDat) < 5 && substr(strDat,1,1) != 'c' )
  {
    strDat <- readline(prompt=paste(idx, "# Validation Data [Prey X, Y, Mouth X,Y,Frame Bout On, Bout Off ] :") )
    strDat <- gsub(strDat,pattern = "[][]",replacement = "")
    arrDat <- strsplit(strDat,"," )[[1]]
  }
  
  
  ## extract validation data on prey location and fish mouth position ##
  strDat <- gsub(strDat,pattern = "[][]",replacement = "")
  arrDat <- strsplit(strDat,"," )[[1]]
  
  if (NROW(arrDat) < 5)
    break;
  
  stopifnot(NROW(arrDat) > 5 ) ##Stop if array data not found
  preyX <- as.numeric( arrDat[1]);  preyY  <- as.numeric(arrDat[2]);
  mouthX <- as.numeric(arrDat[3]); mouthY <- as.numeric(arrDat[4]);
  boutOn <- as.numeric(arrDat[5]) - rec$startFrame; boutOff <- as.numeric(arrDat[6]) - rec$startFrame
  fOnSetEyeVergence <- as.numeric(arrDat[7])
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
   message(paste(recUpd$RegistarIdx,"Distance to prey:", prettyNum(digits=3,recUpd$vMotionBoutDistanceToPrey_mm)," eyeV:",prettyNum(digits=3,recUpd$OnSetEyeVergence) ))
  ## Pass data to record   // SAVE
  strKeyC <- readline(prompt="### Mark Validated ? (y/n):")
  if (strKeyC == 'y')
  {
    recUpd$MarkValidated <- 1
    
    ##update data frame
    datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ] <- recUpd
    
    print(paste( recUpd$RegistarIdx, "# Validated") )
    stopifnot(datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ]$MarkValidated == 1)
    ## Save Back ## 
    ##Update Capture bouts back to original dataframe ##
    datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ]  <- datCaptureBoutsToValidate

    ##Update the register end frame - is this was a capture (last bout) 
    if (recUpd$boutRank == 1)
      datTrackedEventsRegister[idx,]$endFrame <- rec$startFrame + recUpd$vMotionBout_Off
    
    
    ## print event label 
    message("Outcome logged on Registry (matched) :",convertToScoreLabel( datTrackedEventsRegister[idx,]$LabelledScore ))
    print(levels(convertToScoreLabel(rec$LabelledScore )))
    strKeyC <- readline(prompt=paste("### Change label to [",convertToScoreLabel(rec$LabelledScore ),"]:") )
    strKeyC <- as.numeric(strKeyC)
    if (is.numeric(strKeyC))
    {
      newLab <- convertToScoreLabel(as.numeric(strKeyC)-1 ) 
      message(paste("new label:",newLab ) )
      datTrackedEventsRegister[idx,]$LabelledScore <- newLab        
    }
      
    ##Save 
    saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds",sep="")) ##Save With Dataset Idx Identifier
    saveRDS(datTrackedEventsRegister, paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
    
    
  } ##Mark validated
  
 
  strKeyC <- readline(prompt="### Move to next or quit ? (n/q):")
  if (strKeyC == 'q')
    break;
  
  
} ## end of loop 
# datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ]  <- datCaptureBoutsToValidate

## Correct Recs: 
#datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx == 12,]$MarkValidated <- NA

#Save 
saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds",sep="")) ##Save With Dataset Idx Identifier
saveRDS(datTrackedEventsRegister, paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 


stop("Done Labelling")

#### Recalculations Of Bout And 1st Turn To Prey Data - Using retracked Hunt Events####
### Recalc derived bout data / capture speed and angle to prey ###
## loaded file contains datHuntEventMergedFrames
load(file=paste(strDataExportDir,"datAllHuntEventAnalysisFrames_setC.RData",sep=""))


source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r")


idxRegValidated <- datMotionBoutsToValidate[!is.na(datMotionBoutsToValidate$MarkValidated) & datMotionBoutsToValidate$MarkValidated == 1,]$RegistarIdx

#idx <- 79
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
  vEventSpeed_smooth          <- meanf(vDeltaDisplacement[1:NROW(dframe)]/dframe,3) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed_smooth) #meanf(vEventSpeed,100) #
  vEventSpeed_smooth[vEventSpeed_smooth < 0] <- 0 ## Remove -Ve Values As an artefact of Filtering
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventSpeed_smooth_mm <- vFs*vEventSpeed_smooth*DIM_MMPERPX
  
  ## Calc Distance TO Prey, by shifting centroid forward to approx where the Mouth Point is.
  bearingRad = pi/180*(datPlaybackHuntEvent$BodyAngle-90)##+90+180 - Body Heading
  posVX = datPlaybackHuntEvent$posX + cos(bearingRad)*DIM_DISTTOMOUTH_PX/2
  posVY = datPlaybackHuntEvent$posY + sin(bearingRad)*DIM_DISTTOMOUTH_PX/2
  
  
  vDistToPrey          <- sqrt( (posVX -datMotionBouts$vd_PreyX  )^2 + (posVY - datMotionBouts$vd_PreyY)^2   )
  
    ## Get peak Capture Speed within capture bout ##
  for (idxBout in row.names(datMotionBouts) )
  {
    oldVal <- datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm
    idxBoutOn <- datMotionBoutsToValidate[idxBout,]$vMotionBout_On
    idxBoutOff <- min(datMotionBoutsToValidate[idxBout,]$vMotionBout_Off,NROW(vEventSpeed_smooth))
    
    distToPreyInBout     <- vDistToPrey[datMotionBoutsToValidate[idxBout,]$vMotionBout_On:datMotionBoutsToValidate[idxBout,]$vMotionBout_Off]
    ##Find 1st frame wITHIN THIS bOUT where larva closest to prey - usuful when looking at capture frames
    idxTouchDown         <- datMotionBoutsToValidate[idxBout,]$vMotionBout_On + head(which(distToPreyInBout == min(distToPreyInBout,na.rm=T) ),1) 
 
    stopifnot(!is.na(idxTouchDown))
    message(idxTouchDown)
    ##Locate end of capture bout:Find the 1st frame when speed goes below threshold after larva centroid has passed prey item
    ## Re-Clip Capture Bout Around pEak Speed
    idxSpeedLowEND    <- idxTouchDown+min(idxBoutOff-idxTouchDown,which(vEventSpeed_smooth_mm[idxTouchDown:idxBoutOff] < G_THRES_MOTION_BOUT_SPEED),na.rm=T)  ##Find frame where larva closest to prey - assume this is the capture frame
    ##Find Start of the Capture Move, if speed doesnt drop below thresh then assume original boutOn
    idxSpeedLowBEGIN  <- max( idxBoutOn,idxBoutOn+which(vEventSpeed_smooth_mm[idxBoutOn:idxTouchDown] < G_THRES_MOTION_BOUT_SPEED),na.rm=T)
                                     ##Find frame where larva closest to prey - assume this is the capture frame
    ## Debug:   
    plot(vEventSpeed_smooth); points(idxTouchDown,vEventSpeed_smooth[idxTouchDown],col="red"); points(idxBoutOn,vEventSpeed_smooth[idxBoutOn],col="red",cex=3,pch=3); 
    points(idxSpeedLowEND,vEventSpeed_smooth[idxSpeedLowEND],col="blue",pch=13,cex=3)
    points(idxSpeedLowBEGIN,vEventSpeed_smooth[idxSpeedLowBEGIN],col="blue",pch=13,cex=3)
    
    if (!is.integer (idxSpeedLow)) ##If low speed thres not hit, then just take the end of the bout
      idxSpeedLow <- idxBoutOff
    
    frameRangeCapt <- (idxBoutOn:   min(idxBoutOff,idxSpeedLowEND+1 ) ) ##Until Speed is low or end of Bout
    ## MAX Speed over then frames until speed goes down ##
    datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm <- max( vEventSpeed_smooth_mm[idxBoutOn:idxSpeedLowEND] ,na.rm=TRUE )
    stopifnot(is.numeric(datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm) | is.infinite((datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm)))
    
    
    ##Save The Frame/time of peak speed ##
    datMotionBoutsToValidate[idxBout,"vMotionFramesToPeakSpeed"] =  which(vEventSpeed_smooth_mm[idxBoutOn:idxSpeedLowEND] == max( vEventSpeed_smooth_mm[idxBoutOn:idxSpeedLowEND] ,na.rm=TRUE ) )

    ##--Replaced: Save The integrated distance travelled At time of Peak Speed sum(vEventSpeed_smooth[idxBoutOn:idxBoutOn+datMotionBoutsToValidate[idxBout,]$vMotionFramesToPeakSpeed])
    ###  Distance as straight line between OnBout point and Peak Speed point 
    datMotionBoutsToValidate[idxBout,"vMotionPeakSpeed_displacement_px"]<- sqrt( (datPlaybackHuntEvent[idxBoutOn,]$posX - datPlaybackHuntEvent[idxBoutOn+datMotionBoutsToValidate[idxBout,]$vMotionFramesToPeakSpeed,]$posX  )^2 + (datPlaybackHuntEvent[idxBoutOn,]$posY - datPlaybackHuntEvent[idxBoutOn+datMotionBoutsToValidate[idxBout,]$vMotionFramesToPeakSpeed,]$posY )^2   ) 
    stopifnot(is.numeric(datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_displacement_px) | is.infinite((datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_displacement_px)))
    
    ##Save the Frame/time Of min Dist to prey ##
    datMotionBoutsToValidate[idxBout,"PreyMinDistance_frame"] <- idxTouchDown
    
    datPlaybackHuntEvent[datMotionBouts$vMotionBout_On:datMotionBouts$vMotionBout_Off, ]$PreyID <- 1
    datPlaybackHuntEvent[datMotionBouts$vMotionBout_On:datMotionBouts$vMotionBout_Off, ]$Prey_X <- datMotionBouts$vd_PreyX  
    datPlaybackHuntEvent[datMotionBouts$vMotionBout_On:datMotionBouts$vMotionBout_Off, ]$Prey_Y <- datMotionBouts$vd_PreyY
    
    ##Calc relative angle of mouth to Prey , and normalize to Vertical axis by removing the  body axis angle
    bodyAngle <- datPlaybackHuntEvent[datMotionBouts$vMotionBout_On,]$BodyAngle 

    #relAngle <- ( ( 180 +  180/pi * atan2(datMotionBouts$vd_PreyX - datMotionBouts$vd_MouthX,datMotionBouts$vd_MouthY - datMotionBouts$vd_PreyY  )) - bodyAngle    ) %% 360 - 180
    relAngle <- 180- 180/pi * atan2(datMotionBouts$vd_PreyX - datMotionBouts$vd_MouthX, datMotionBouts$vd_PreyY-datMotionBouts$vd_MouthY  ) 
    ##Check which direction the angles should be compared
    {
      if ( abs(bodyAngle - relAngle) > 180) ##Body Angle Leads - Substract
        relAngle = ifelse(relAngle < bodyAngle, 360+relAngle-bodyAngle,relAngle - (360+bodyAngle) )
      else 
        relAngle <- relAngle-bodyAngle
    }
    ##Convert atan2 angle to 360 circle
    #relAngle <- ifelse(relAngle < 0,relAngle + 360,relAngle)
    ##Now Get relative Angle to larva body orientation (also measured in 360 to vertical axis)    
    #relAngle <-  ifelse(relAngle > 180 & bodyAngle < 180,(relAngle) - (360+bodyAngle), relAngle) ##If  comparing angles across the vertical 0 line
    #relAngle <-  ifelse(bodyAngle > 180 & relAngle < 180,(360+relAngle) - (bodyAngle), (relAngle) - (bodyAngle)) ##If  comparing angles across the vertical 0 line
    
    if (abs(relAngle) > 160 )
      stop("Angle to prey over 160 degrees")
    
    ##Set Them Back in The Bout INfo
    oldValAngle <- datMotionBoutsToValidate[idxBout,]$OnSetAngleToPrey
    datMotionBoutsToValidate[idxBout,]$OnSetAngleToPrey  <- relAngle
    
    datMotionBoutsToValidate[idxBout,]$OffSetAngleToPrey <- NA # relAngle[datMotionBouts$vMotionBout_Off]
    
    
    oldDist <- datMotionBoutsToValidate[idxBout,]$vMotionBoutDistanceToPrey_mm     
    ##Calc Distance From Mouth
    datMotionBoutsToValidate[idxBout,]$vMotionBoutDistanceToPrey_mm <-  DIM_MMPERPX*(sqrt((datMotionBoutsToValidate[idxBout,]$vd_PreyX-datMotionBoutsToValidate[idxBout,]$vd_MouthX )^2+(datMotionBoutsToValidate[idxBout,]$vd_PreyY-datMotionBoutsToValidate[idxBout,]$vd_MouthY)^2 ))
    if (is.na(datMotionBoutsToValidate[idxBout,]$vMotionBoutDistanceToPrey_mm))
      stop("NA Dist")
    
    print( paste(row.names(recReg),"speed",oldVal,"->",datMotionBoutsToValidate[idxBout,]$vMotionPeakSpeed_mm,
                "Prey Angle: ",oldValAngle,"->",datMotionBoutsToValidate[idxBout,]$OnSetAngleToPrey,
                "Prey Dist: ",oldDist,"->",datMotionBoutsToValidate[idxBout,]$vMotionBoutDistanceToPrey_mm) )
    
  }
  
  ## Do Next Validated rec ##
}

## DUBUG CODE Testing  mouth-Prey Position and Angles 
##Make relative to y axis instead of x - and convert to clockwise counting
datMotionBouts <- datMotionBoutsToValidate[idxBout,]
d <- datMotionBoutsToValidate[idxBout,]$vMotionBoutDistanceToPrey_mm/DIM_MMPERPX
xAngle <- 180-180/pi * atan2(datMotionBouts$vd_PreyX - datMotionBouts$vd_MouthX, datMotionBouts$vd_PreyY-datMotionBouts$vd_MouthY  ) 
##Make atan2 angle go around 360 circle so we can compare to body angle

relAngle <- ifelse(xAngle < 0,xAngle + 360,xAngle)
##Check which direction the angles should be compared
{
if ( abs(bodyAngle - relAngle) > 180) ##Body Angle Leads - Substract
  relAngle = ifelse(relAngle < bodyAngle, 360+relAngle-bodyAngle,relAngle - (360+bodyAngle) )
else 
  relAngle <- relAngle-bodyAngle
}
  
##Now Get relative Angle to larva body orientation (also measured in 360 to vertical axis)    
#relAngle <-  ifelse(relAngle > 180 & bodyAngle < 180,(relAngle) - (360+bodyAngle), relAngle) ##If  comparing angles across the vertical 0 line
#relAngle <-  ifelse(bodyAngle > 180 & relAngle < 180,(360+relAngle) - (bodyAngle), (relAngle) - (bodyAngle)) ##If  comparing angles across the vertical 0 line

x <- (10)*cos(pi/180 * bodyAngle -pi/2) + datMotionBouts$vd_MouthX
y <- (10)*sin(pi/180 * bodyAngle - pi/2) + datMotionBouts$vd_MouthY
xp <- (d)*cos(pi/180 * (bodyAngle+relAngle) - pi/2) + datMotionBouts$vd_MouthX
yp <- (d)*sin(pi/180 * (bodyAngle+relAngle) - pi/2) + datMotionBouts$vd_MouthY

plot(datMotionBouts$vd_PreyX,512-datMotionBouts$vd_PreyY,xlim=c(datMotionBouts$vd_MouthX-50,datMotionBouts$vd_MouthX+50)) ## Prey
points(datMotionBouts$vd_MouthX,512-datMotionBouts$vd_MouthY,col="blue",pch=16) ##Mouth
segments(datMotionBouts$vd_MouthX,512-datMotionBouts$vd_MouthY,x,512-y,cex=1, main="",col="red") ##Body Orientation
segments(datMotionBouts$vd_MouthX,512-datMotionBouts$vd_MouthY,xp,512-yp,cex=1, main="",col="magenta")

##Update Saved Data
saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_peakSpeedExtended_Validated.rds",sep="")) ##Save With Dataset Idx Identifier

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
                                  ColisionFrame = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                          datMotionBoutCombined$boutRank == 1 ,]$PreyMinDistance_frame, ##Min Dist to prey (centroid) Hit time in frames
                                  NFramesToPeakSpeed = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                          datMotionBoutCombined$boutRank == 1 ,]$vMotionFramesToPeakSpeed, ## Time until speed Peaks (Accelleration)
                                  PeakSpeedDistance = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                           datMotionBoutCombined$boutRank == 1 ,]$vMotionPeakSpeed_displacement_px*DIM_MMPERPX, ##Dist travelled until peak speed 
                                  
                                  DistanceToPrey = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                            datMotionBoutCombined$boutRank == 1 ,]$vMotionBoutDistanceToPrey_mm,
                                  doesCaptureStrike=( datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                               datMotionBoutCombined$boutRank == 1 ,]$vMotionPeakSpeed_mm >= G_THRES_CAPTURE_SPEED ),
                                  CaptureBoutStartFrame=( datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                                datMotionBoutCombined$boutRank == 1 ,]$vMotionBout_On  ),
                                  CaptureBoutEndFrame=( datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                                datMotionBoutCombined$boutRank == 1 ,]$vMotionBout_Off  ),
                                  
                                  CaptureStrikeEyeVergence = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                                      datMotionBoutCombined$boutRank == 1 ,]$OnSetEyeVergence,
                                  Validated = datMotionBoutCombined[ datMotionBoutCombined$RegistarIdx %in% datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx  &
                                                                       datMotionBoutCombined$boutRank == 1 ,]$MarkValidated
  )
  
}

##Save List on First Bout Data
saveRDS(lFirstBoutPoints,file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated",".rds",sep="") ) #Processed Registry on which we add )

###Need to Re-Run Merging with Cluster ID


### Load Pre Calc Results
load(file =paste(strDataExportDir,"stat_CaptSpeedCluster_RJags.RData",sep=""))
#### Main Figure 4 - Show Distance Vs Capture speed clusters for all groups - and Prob Of Capture Strike###

##Update The Capture  Bout Data list with the new clustering (huntEpisodeAnalysis_FirstBoutData)
drawClust <- list(NF=draw_NF,LF=draw_LF,DF=draw_DF)
makeCaptureClusteredData(lFirstBoutPoints,drawClust)


###get distance travelled at peak speed - ie Accellaration - 
## Peak speed should be close to prey position - Which We can compare to the distance to prey 




