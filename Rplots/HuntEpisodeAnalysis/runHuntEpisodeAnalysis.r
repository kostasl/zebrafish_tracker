library(signal)
library(MASS)
library(mclust,quietly = TRUE)

source("TrackerDataFilesImport_lib.r")
source("plotTrackScatterAndDensities.r")


strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis",".RData",sep="") ##To Which To Save After Loading
message(paste(" Importing Retracked HuntEvents from:",strDataFileName))

#for (i in 1:20) dev.off()

#
############# LOAD AND PLAYBACK OF HUNT EVENTS ####
load(strDataFileName)
##Test  PlayBack Plot Hunt Event###  
##Make an Updated list of ReTracked Hunt Events that have been imported
# datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))

## Setup Filters ## Can Check Bands with freqz(bf_speed)
Fs <- 430; #sampling rate
bf_tail <- butter(1, c(0.001,0.1),type="pass");
bf_eyes <- butter(4, 0.025,type="low",plane="z");
bf_speed <- butter(4, 0.05,type="low");  
###
nEyeFilterWidth <- nFrWidth*8 ##For Median Filtering


idxH <- 20
expID <- datTrackedEventsRegister[idxH,]$expID
trackID<- datTrackedEventsRegister[idxH,]$trackID
eventID <- datTrackedEventsRegister[idxH,]$eventID
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
#X11()
#spectrum(datRenderHuntEvent$LEyeAngle)
datRenderHuntEvent$LEyeAngle <-filtfilt(bf_eyes,datRenderHuntEvent$LEyeAngle) # filtfilt(bf_eyes, medianf(datRenderHuntEvent$LEyeAngle,nFrWidth)) #meanf(datHuntEventMergedFrames$LEyeAngle,20)
#X11()
#spectrum(datRenderHuntEvent$LEyeAngle)

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
X11()
plot(datRenderHuntEvent$frameN,datRenderHuntEvent$LEyeAngle,type='l',col="blue",ylim=c(-60,60),main="Eye Motion ")
lines(datRenderHuntEvent$frameN,datRenderHuntEvent$REyeAngle,type='l',col="magenta")


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


#renderHuntEventPlayback(datRenderHuntEvent,speed=1) #saveToFolder =  strFolderName

#### PROCESS BOUTS ###
vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2)*DIM_MMPERPX ## Path Length Calculated As Total Displacement
#nNumberOfBouts       <- 
dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
dframe               <- dframe[dframe > 0]/Fs ##Clear Any possible Nan - and Convert To Time sec  
dEventSpeed          <- meanf(vEventPathLength/dframe,3) ##Apply Mean Filter Smooth Out 


#speed_Smoothed <- meanf(dEventSpeed,10)
##Replace NA with 0s
dEventSpeed[is.na(dEventSpeed)] = 0
dEventSpeed_smooth <- filtfilt(bf_speed, dEventSpeed) #meanf(dEventSpeed,100) #
dEventSpeed_smooth[is.na(dEventSpeed_smooth)] = 0
MoveboutsIdx <- detectMotionBouts(dEventSpeed_smooth)##find_peaks(dEventSpeed_smooth*100,25)
##Reject Peaks Below Half An SD Peak Value - So As to Choose Only Significant Bout Movements # That Are Above the Minimum Speed to Consider As Bout
#MoveboutsIdx_cleaned <- MoveboutsIdx[which(dEventSpeed_smooth[MoveboutsIdx] > sd(dEventSpeed_smooth[MoveboutsIdx])/3 
                                           #& dEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
MoveboutsIdx_cleaned <-MoveboutsIdx #[which(dEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]

##Binarize , Use indicator function 1/0 for frames where Motion Occurs
vMotionBout <- dEventSpeed_smooth
vMotionBout[ 1:NROW(vMotionBout) ] = 0
vMotionBout[ MoveboutsIdx_cleaned  ] = 1
vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
vMotionBout_rle <- rle(vMotionBout)
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

##Get Bout Statistics ##
vMotionBoutDuration_msec <- vMotionBout_Off[1:iPairs]-vMotionBout_On[1:iPairs]
vMotionBoutDuration_msec <- 1000*vMotionBoutDuration_msec[!is.na(vMotionBoutDuration_msec)]/Fs
vMotionBoutIntervals_msec <- 1000*(vMotionBout_On[3:(iPairs)] - vMotionBout_Off[2:(iPairs-1)])/Fs


## Take InterBoutIntervals in msec from Last to first
vMotionBout_rle <- rle(vMotionBout)
vMotionBoutIBI <-1000*vMotionBout_rle$lengths[seq(NROW(vMotionBout_rle$lengths)-2,1,-2 )]/Fs
vMotionBoutDuration <-1000*vMotionBout_rle$lengths[seq(NROW(vMotionBout_rle$lengths)-1,2,-2 )]/Fs
stopifnot(vMotionBout_rle$values[NROW(vMotionBout_rle$lengths)] == 0 )
stopifnot(vMotionBout_rle$values[3] == 0 ) ##THe INitial vMotionBoutIBI Is not Actually A pause interval , but belongs to motion!
##On Bout Lengths
##Where t=0 is the capture bout, -1 -2 are the steps leading to it
vMotionPeriod <- 1-(NROW(vMotionBoutIBI)-seq(NROW(vMotionBoutIBI),1 ))
vMotionBoutIBI <- vMotionBoutIBI

## Plot Durations of Pause/Go
X11()
plot(0.0+NROW(vMotionBoutDuration)-seq(NROW(vMotionBoutDuration),1 ),vMotionBoutDuration,
     xlab="Bout",ylab="msec",xlim=c(0,10),ylim=c(0,500),
     col="red",main="Bout Duration",pch=16) ##Take Every Bout Length
points(0.5+NROW(vMotionBoutIBI)-seq(NROW(vMotionBoutIBI),1 ),vMotionBoutIBI,col="blue",pch=21) ##Take every period between / Inter Bout Interval
legend(1,400,c("Motion","Pause" ),col=c("red","blue"),pch=c(16,21) )







## Plot The Start Stop Motion Bout Binarized Data
X11()
plot(vMotionBout,type='p')
points(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red")
points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green") ##On
points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="yellow")##Off

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
#lines(dEventSpeed_smooth*50,type='l',col="blue")
##plot Correlation Of Tail Movement To speed 
#corr_speedVsTail <- ccf(abs(vTailDisp),dEventSpeed_smooth,type="correlation",plot=TRUE)

#llRange <- min(NROW(abs(dEventSpeed_smooth)),NROW(abs(vTailDisp))) 
#cor_TailToSpeed <- cov(abs(vTailDisp[1:llRange]),dEventSpeed_smooth[1:llRange])

#X11()
#plot(abs(dEventSpeed_smooth[1:llRange]) , abs(vTailDisp[1:llRange]) , type="p")

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


##Plot Displacement and Speed(Scaled)
X11()
plot(cumsum(vEventPathLength),ylab="mm/sec",ylim=c(0,max(dEventSpeed_smooth))) ##PLot Total Displacemnt over time
lines(dEventSpeed_smooth,type='l',col="blue")
lines(vTailDispFilt,type='l',col="magenta")
points(MoveboutsIdx,dEventSpeed_smooth[MoveboutsIdx],col="black")
points(MoveboutsIdx_cleaned,dEventSpeed_smooth[MoveboutsIdx_cleaned],col="red")
message(paste("Number oF Bouts:",length(MoveboutsIdx_cleaned)))
dev.copy(png,filename=paste(strPlotExportPath,"/Movement-Bout_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
dev.off()

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

# for (i in 1:20) dev.off()

