
### Hunt Episode Data Extraction and Analysis Scripts 
### Kostasl Jul 2018 
## It is aimed at identified hunt event that have been carefully retracked and imported separatelly from these new Tracked CsV files using runImportHuntEpisodeTrackFiles
## It provides the data for identifying Bouts using FishSpeed. 
## Notes: 
## * Gaussian Mixt. Clusters are used to Identify Bout from Speeds, however the BoutOnset-Offset Is then Fixed To Identify When Speed Ramps Up , 
## and final decelpartion is used to obtain  more accurate estimates of BoutDurations-
## * The Distance to Prey Is calculated But Also Interpolated using Fish Displacement where there are missing values - 
#####


library(signal)
library(MASS)
library(mclust,quietly = TRUE)

source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r")
source("TrackerDataFilesImport_lib.r")
source("plotTrackScatterAndDensities.r")


strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis",".RData",sep="") ##To Which To Save After Loading
message(paste(" Importing Retracked HuntEvents from:",strDataFileName))



#    for (i in 1:40) dev.off()

#
############# LOAD AND PLAYBACK OF HUNT EVENTS ####
load(strDataFileName)
##Test  PlayBack Plot Hunt Event###  
##Make an Updated list of ReTracked Hunt Events that have been imported
# datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))

## Setup Filters ## Can Check Bands with freqz(bf_speed)
Fs <- 430; #sampling rate
bf_tail <- butter(1, 0.2,type="low");
bf_eyes <- butter(4, 0.025,type="low",plane="z");
bf_speed <- butter(4, 0.05,type="low");  ##Focus On Low Fq to improve Detection Of Bout Motion and not little Jitter motion
###
nEyeFilterWidth <- nFrWidth*8 ##For Median Filtering

lMotionBoutDat <- list()


#idxH <- 20

for (idxH in 1:NROW(datTrackedEventsRegister))
{
  
  expID <- datTrackedEventsRegister[idxH,]$expID
  trackID<- datTrackedEventsRegister[idxH,]$trackID
  eventID <- datTrackedEventsRegister[idxH,]$eventID
  groupID <- datTrackedEventsRegister[idxH,]$groupID
  selectedPreyID <- datTrackedEventsRegister[idxH,]$PreyIDTarget
  
  message(paste(idxH, ".Process Hunt Event Expid:",expID,"Event:",eventID))
  
  datRenderHuntEvent <- datHuntEventMergedFrames[datHuntEventMergedFrames$expID==expID 
                                                 & datHuntEventMergedFrames$trackID==trackID 
                                                 & datHuntEventMergedFrames$eventID==eventID,]
  
  strFolderName <- paste( strPlotExportPath,"/renderedHuntEvent",expID,"_event",eventID,"_track",trackID,sep="" )
  #dir.create(strFolderName )
  ##Remove NAs

  lMax <- 55
  lMin <- -20
  #spectrum(datRenderHuntEvent$LEyeAngle)
  datRenderHuntEvent$LEyeAngle <- clipEyeRange(datRenderHuntEvent$LEyeAngle,lMin,lMax)
  datRenderHuntEvent$LEyeAngle <-medianf(datRenderHuntEvent$LEyeAngle,nEyeFilterWidth)
  datRenderHuntEvent$LEyeAngle[is.na(datRenderHuntEvent$LEyeAngle)] <- 0
  datRenderHuntEvent$LEyeAngle <-filtfilt(bf_eyes,datRenderHuntEvent$LEyeAngle) # filtfilt(bf_eyes, medianf(datRenderHuntEvent$LEyeAngle,nFrWidth)) #meanf(datHuntEventMergedFrames$LEyeAngle,20)
  
  #X11()
  #lines(medianf(datRenderHuntEvent$LEyeAngle,nFrWidth),col='red')
  #lines(datRenderHuntEvent$LEyeAngle,type='l',col='blue')
  ##Replace Tracking Errors (Values set to 180) with previous last known value
  
  lMax <- 15
  lMin <- -50
  datRenderHuntEvent$REyeAngle <- clipEyeRange(datRenderHuntEvent$REyeAngle,lMin,lMax)
  datRenderHuntEvent$REyeAngle <-medianf(datRenderHuntEvent$REyeAngle,nEyeFilterWidth)
  datRenderHuntEvent$REyeAngle[is.na(datRenderHuntEvent$REyeAngle)] <- 0
  datRenderHuntEvent$REyeAngle <- filtfilt(bf_eyes,datRenderHuntEvent$REyeAngle  ) #meanf(datHuntEventMergedFrames$REyeAngle,20)
  #datRenderHuntEvent$REyeAngle <-medianf(datRenderHuntEvent$REyeAngle,nFrWidth)
  
  
  
  lMax <- +50
  lMin <- -50
  datRenderHuntEvent$DThetaSpine_7 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_7,lMin,lMax) )
  datRenderHuntEvent$DThetaSpine_6 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_6,lMin,lMax) )
  datRenderHuntEvent$DThetaSpine_5 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_5,lMin,lMax) )
  datRenderHuntEvent$DThetaSpine_4 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_4,lMin,lMax))
  datRenderHuntEvent$DThetaSpine_3 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_3,lMin,lMax))
  datRenderHuntEvent$DThetaSpine_2 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_2,lMin,lMax))
  datRenderHuntEvent$DThetaSpine_1 <- filtfilt(bf_tail, clipEyeRange(datRenderHuntEvent$DThetaSpine_1,lMin,lMax))
  
  rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
  r <- c(rfc(8),"#FF0000");
  
  
  
  ## PLAYBACK ####
  #        renderHuntEventPlayback(datRenderHuntEvent,speed=1) #saveToFolder =  strFolderName
  ########################
  
  ############ PREY SELECTION #####
  ## Begin Data EXtraction ###
  ## Remove Multiple Prey Targets ########
  ##Get Number oF Records per Prey
  tblPreyRecord <-table(datRenderHuntEvent$PreyID) 
  if (NROW(tblPreyRecord) > 1)
    warning("Multiple Prey Items Tracked In Hunt Episode-Selecting Longest Track")
  
  
  ##Add PreyTarget ID To Register Save The Prey Target To Register
  if (!any(names(datTrackedEventsRegister) == "PreyIDTarget"))
    datTrackedEventsRegister$PreyIDTarget <- NA

  #selectedPreyID <- max(as.numeric(names(which(tblPreyRecord == max(tblPreyRecord)))))
  ##Check If Assigned OtherWise Automatically Select the longest Track
  if (is.na(datTrackedEventsRegister[idxH,]$PreyIDTarget)) 
  {
    selectedPreyID <-  max(as.numeric(names(which(tblPreyRecord == max(tblPreyRecord)))))
    datTrackedEventsRegister[idxH,]$PreyIDTarget <- selectedPreyID
    save(datHuntEventMergedFrames,datTrackedEventsRegister,lHuntEventTRACKSfileSrc,lHuntEventFOODfileSrc,file=strDataFileName) ##Save With Dataset Idx Identifier
  }
  
  #
  ##Save The Selected Prey Item
  if (selectedPreyID != datTrackedEventsRegister[idxH,]$PreyIDTarget)
  {
    message(paste("Targeted Prey Changed For Hunt Event to ID:",selectedPreyID," - Updating Register...") )
    datTrackedEventsRegister[idxH,]$PreyIDTarget <- selectedPreyID
    save(datHuntEventMergedFrames,datTrackedEventsRegister,lHuntEventTRACKSfileSrc,lHuntEventFOODfileSrc,file=strDataFileName) ##Save With Dataset Idx Identifier
  }
  ################ END OF PREY SELECT / Start Processing ###
  
  
  ##Select Prey Specific Subset
  datRenderHuntEvent <- datRenderHuntEvent[datRenderHuntEvent$PreyID == selectedPreyID,] 
  
  
  #### PROCESS BOUTS ###
  vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
  vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
  vDeltaDisplacement   <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
  #nNumberOfBouts       <- 
  dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
  dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - and Convert To Time sec  
  vEventSpeed          <- meanf(vDeltaDisplacement/dframe,3) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
  #vEventPathLength     <- cumsum(vEventSpeed) ##Noise Adds to Length
  vDistToPrey          <- meanf(sqrt( (datRenderHuntEvent$Prey_X -datRenderHuntEvent$posX )^2 + (datRenderHuntEvent$Prey_Y -datRenderHuntEvent$posY)^2   ),3)
  vSpeedToPrey         <- diff(vDistToPrey,lag=1,differences=1)
  
  #speed_Smoothed <- meanf(vEventSpeed,10)
  ##Replace NA with 0s
  vEventSpeed[is.na(vEventSpeed)] = 0
  vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed) #meanf(vEventSpeed,100) #
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventPathLength <- cumsum(vEventSpeed_smooth)
  MoveboutsIdx <- detectMotionBouts(vEventSpeed)##find_peaks(vEventSpeed_smooth*100,25)
  
  MoveboutsIdx_cleaned <- MoveboutsIdx# which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   ) #MoveboutsIdx# 
  
  ##Distance To PRey
  vDistToPrey_Fixed      <- interpolateDistToPrey(vDistToPrey,vEventSpeed_smooth)
  lMotionBoutDat[[idxH]]  <- calcMotionBoutInfo(MoveboutsIdx_cleaned,vEventSpeed_smooth,vDistToPrey_Fixed)
  rows <- NROW(lMotionBoutDat[[idxH]])
  lMotionBoutDat[[idxH]] <- cbind(lMotionBoutDat[[idxH]] ,RegistarIdx = rep(idxH,rows),expID=rep(expID,rows),eventID=rep(eventID,rows),groupID=rep(as.character(groupID),rows))
}  

datEpisodeMotionBout <- lMotionBoutDat[[1]]
##On Bout Lengths
##Where t=0 is the capture bout, -1 -2 are the steps leading to it


## Plot Durations of Pause/Go
X11()
plot(datEpisodeMotionBout[,"vMotionBoutDuration"],
     xlab="Bout",ylab="msec",xlim=c(0,NROW(datEpisodeMotionBout) ),ylim=c(0,500),
     col="red",main="Bout Duration",pch=16) ##Take Every Bout Length
points(datEpisodeMotionBout[,"vMotionBoutIBI"],col="blue",pch=21) ##Take every period between / Inter Bout Interval
legend(1,400,c("Motion","Pause" ),col=c("red","blue"),pch=c(16,21) )

X11()
plot(datEpisodeMotionBout[,"vMotionBoutDistanceToPrey_mm"],
     xlab="Bout",ylab="mm",xlim=c(0,NROW(datEpisodeMotionBout)),ylim=c(0,3),
     col="red",main="Bout Distance To Prey",pch=16) ##Take Every Bout Length

X11()
plot(datEpisodeMotionBout[,"vMotionBoutDistanceTravelled_mm"],
     xlab="Bout",ylab="mm",xlim=c(0,NROW(datEpisodeMotionBout)),ylim=c(0,2),
     col="red",main="Bout Power",pch=16) ##Take Every Bout Length


X11()
plot(datRenderHuntEvent$frameN,datRenderHuntEvent$LEyeAngle,type='l',col="blue",ylim=c(-60,60),main="Eye Motion ")
lines(datRenderHuntEvent$frameN,datRenderHuntEvent$REyeAngle,type='l',col="magenta")

## Plot The Start Stop Motion Bout Binarized Data
#X11()
#plot(vMotionBout,type='p')
#points(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red")
#points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green",pch=7) ##On
#points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="yellow",pch=21)##Off

######### END OF PROCESS BOUT #########

## ## Tail Curvature 
##Filter The Noise 
#when apply twice with filtfilt, #results in a 0 phase shift  : W * (Fs/2) == half-amplitude cut-off when combined with filtfilt

vTailDir <-  datRenderHuntEvent$DThetaSpine_1 +  datRenderHuntEvent$DThetaSpine_2 + datRenderHuntEvent$DThetaSpine_3 + datRenderHuntEvent$DThetaSpine_4 + datRenderHuntEvent$DThetaSpine_5 + datRenderHuntEvent$DThetaSpine_6 + datRenderHuntEvent$DThetaSpine_7
vTailDisp <-  datRenderHuntEvent$DThetaSpine_6 + datRenderHuntEvent$DThetaSpine_7 #+ datRenderHuntEvent$DThetaSpine_7 #+ datRenderHuntEvent$DThetaSpine_7 #abs(datRenderHuntEvent$DThetaSpine_1) +  abs(datRenderHuntEvent$DThetaSpine_2) + abs(datRenderHuntEvent$DThetaSpine_3) + abs(datRenderHuntEvent$DThetaSpine_4) + abs(datRenderHuntEvent$DThetaSpine_5) + abs(datRenderHuntEvent$DThetaSpine_6) + abs(datRenderHuntEvent$DThetaSpine_7)
vTailDispFilt <- filtfilt(bf_tail, vTailDisp)

#X11()
#layout(matrix(c(1,2), 2, 1, byrow = TRUE))
#plot(vTailDisp,type="l")
#lines(vEventSpeed_smooth*50,type='l',col="blue")
##plot Correlation Of Tail Movement To speed 
#corr_speedVsTail <- ccf(abs(vTailDisp),vEventSpeed_smooth,type="correlation",plot=TRUE)

#llRange <- min(NROW(abs(vEventSpeed_smooth)),NROW(abs(vTailDisp))) 
#cor_TailToSpeed <- cov(abs(vTailDisp[1:llRange]),vEventSpeed_smooth[1:llRange])

#X11()
#plot(abs(vEventSpeed_smooth[1:llRange]) , abs(vTailDisp[1:llRange]) , type="p")

#lines(vTailDir,type='l',col="green")
##END OF CURVATURE ##

##Plot Tail Segments Displacements #
#X11()
#plot(datRenderHuntEvent$DThetaSpine_1,type='l',col=r[1])
#lines(datRenderHuntEvent$DThetaSpine_2,type='l',col=r[2])
#lines(datRenderHuntEvent$DThetaSpine_3,type='l',col=r[3])
#lines(datRenderHuntEvent$DThetaSpine_4,type='l',col=r[4])
#lines(datRenderHuntEvent$DThetaSpine_5,type='l',col=r[5])
#lines(datRenderHuntEvent$DThetaSpine_6,type='l',col=r[6])
#lines(datRenderHuntEvent$DThetaSpine_7,type='l',col=r[7])



###########  Plot Polar Angle to Prey ##############
X11()
plot.new()
polarPlotAngleToPrey(datRenderHuntEvent)
dev.copy(png,filename=paste(strPlotExportPath,"/AngleToPreyVsTime_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
dev.off()

X11()
plot.new()
polarPlotAngleToPreyVsDistance(datRenderHuntEvent)
dev.copy(png,filename=paste(strPlotExportPath,"/AngleToPreyVsDistance_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
dev.off()
###################################################



datMotionBoutCombined <-  data.frame( do.call(rbind,lMotionBoutDat ) )
 
X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDistanceToPrey_mm,main="Distance From Prey",ylab="mm")
boxplot(datMotionBoutCombined$vMotionBoutDistanceToPrey_mm ~ datMotionBoutCombined$boutRank,main="Distance From Prey",ylab="mm",xlab="Bout Sequence (From Capture - Backwards)")

X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDistanceTravelled_mm,main="Distance Of Bout (power)",ylab="mm")
boxplot(datMotionBoutCombined$vMotionBoutDistanceTravelled_mm ~ datMotionBoutCombined$boutRank,main="Distance Of Bout (power)",ylab="mm",xlab="Bout Sequence (From Capture - Backwards)")


X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDuration,main=" Bout Duration",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")
boxplot(datMotionBoutCombined$vMotionBoutDuration ~ datMotionBoutCombined$boutRank,main=" Bout Duration",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")

X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutIBI,main=" Inter Bout Intervals ",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")
boxplot( datMotionBoutCombined$vMotionBoutIBI ~ datMotionBoutCombined$boutRank,main=" Inter Bout Intervals ",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")
# for (i in 1:20) dev.off()

