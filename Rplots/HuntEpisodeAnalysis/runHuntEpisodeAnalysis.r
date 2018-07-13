library(signal)
#### Analyse Extracted/Labelled and then Retracked Hunt Events ///




strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis",".RData",sep="") ##To Which To Save After Loading
message(paste(" Importing Retracked HuntEvents from:",strDataFileName))


#
############# LOAD AND PLAYBACK OF HUNT EVENTS ####
load(strDataFileName)
##Test  PlayBack Plot Hunt Event###  
##Make an Updated list of ReTracked Hunt Events that have been imported
# datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))

## Setup Filters ## Can Check Bands with freqz(bf_speed)
Fs <- 430; #sampling rate
bf_tail <- butter(4, c(0.02,0.2),type="pass");
bf_eyes <- butter(4, 0.015,type="low",plane="z");
bf_speed <- butter(4, 0.05,type="low");  
###


idxH <- 14
expID <- datTrackedEventsRegister[idxH,]$expID
trackID<- datTrackedEventsRegister[idxH,]$trackID
eventID <- datTrackedEventsRegister[idxH,]$eventID
datRenderHuntEvent <- datHuntEventMergedFrames[datHuntEventMergedFrames$expID==expID 
                                               & datHuntEventMergedFrames$trackID==trackID 
                                               & datHuntEventMergedFrames$eventID==eventID,]

strFolderName <- paste( strPlotExportPath,"/renderedHuntEvent",expID,"_event",eventID,"_track",trackID,sep="" )
#dir.create(strFolderName )
##Remove NAs

X11()
plot(datRenderHuntEvent$LEyeAngle,type='l')
lines(medianf(datRenderHuntEvent$LEyeAngle,nFrWidth),col='red')

#spectrum(datRenderHuntEvent$LEyeAngle)

datRenderHuntEvent$LEyeAngle <-medianf(datRenderHuntEvent$LEyeAngle,nFrWidth)
datRenderHuntEvent$LEyeAngle[is.na(datRenderHuntEvent$LEyeAngle)] <- 0
#X11()
#spectrum(datRenderHuntEvent$LEyeAngle)
datRenderHuntEvent$LEyeAngle <-filtfilt(bf_eyes,datRenderHuntEvent$LEyeAngle) # filtfilt(bf_eyes, medianf(datRenderHuntEvent$LEyeAngle,nFrWidth)) #meanf(datHuntEventMergedFrames$LEyeAngle,20)
#X11()
#spectrum(datRenderHuntEvent$LEyeAngle)

#X11()
lines(datRenderHuntEvent$LEyeAngle,type='l',col='blue')

datRenderHuntEvent$REyeAngle <-medianf(datRenderHuntEvent$REyeAngle,nFrWidth)
datRenderHuntEvent$REyeAngle[is.na(datRenderHuntEvent$REyeAngle)] <- 0
datRenderHuntEvent$REyeAngle <- filtfilt(bf_eyes,datRenderHuntEvent$REyeAngle  ) #meanf(datHuntEventMergedFrames$REyeAngle,20)
datRenderHuntEvent$REyeAngle <-medianf(datRenderHuntEvent$REyeAngle,nFrWidth)

datRenderHuntEvent$DThetaSpine_7 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_7)
datRenderHuntEvent$DThetaSpine_6 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_6)
datRenderHuntEvent$DThetaSpine_5 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_5)
datRenderHuntEvent$DThetaSpine_4 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_4)
datRenderHuntEvent$DThetaSpine_3 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_3)
datRenderHuntEvent$DThetaSpine_2 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_2)
datRenderHuntEvent$DThetaSpine_1 <- filtfilt(bf_tail, datRenderHuntEvent$DThetaSpine_1)

rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
r <- c(rfc(8),"#FF0000");


renderHuntEventPlayback(datRenderHuntEvent,speed=1) #saveToFolder =  strFolderName

#### PROCESS BOUTS ###




vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
#nNumberOfBouts       <- 
dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - Why is dFrame 0?  
dEventSpeed          <- meanf(vEventPathLength/dframe,3) ##Apply Mean Filter Smooth Out 

#speed_Smoothed <- meanf(dEventSpeed,10)
##Replace NA with 0s
dEventSpeed[is.na(dEventSpeed)] = 0
dEventSpeed_smooth <- filtfilt(bf_speed, dEventSpeed) #meanf(dEventSpeed,100) #
dEventSpeed_smooth[is.na(dEventSpeed_smooth)] = 0
MoveboutsIdx <- find_peaks(dEventSpeed_smooth*100,25)
##Reject Peaks Below Half An SD Peak Value - So As to Choose Only Significant Bout Movements # That Are Above the Minimum Speed to Consider As Bout
MoveboutsIdx_cleaned <- MoveboutsIdx[which(dEventSpeed_smooth[MoveboutsIdx] > sd(dEventSpeed_smooth[MoveboutsIdx])/4 
                                           & dEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]


##Binarize , Use indicator function 1/0 for frames where Motion Occurs
vMotionBout <- dEventSpeed_smooth
vMotionBout[ vMotionBout < G_MIN_BOUTSPEED  ] = 0
vMotionBout[vMotionBout > G_MIN_BOUTSPEED  ] = 1
vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
rle(vMotionBout)
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

X11()
plot(vMotionBout,type='p')
points(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red")
points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green") ##On
points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="yellow")##Off

######### END OF PROCESS BOUT #########

## ## Tail Curvature 
##Filter The Noise 
#when apply twice with filtfilt, #results in a 0 phase shift  : W * (Fs/2) == half-amplitude cut-off when combined with filtfilt
X11()
vTailDir <-  datRenderHuntEvent$DThetaSpine_1 +  datRenderHuntEvent$DThetaSpine_2 + datRenderHuntEvent$DThetaSpine_3 + datRenderHuntEvent$DThetaSpine_4 + datRenderHuntEvent$DThetaSpine_5 + datRenderHuntEvent$DThetaSpine_6 + datRenderHuntEvent$DThetaSpine_7
vTailDisp <-  datRenderHuntEvent$DThetaSpine_6 #+ datRenderHuntEvent$DThetaSpine_7 #+ datRenderHuntEvent$DThetaSpine_7 #abs(datRenderHuntEvent$DThetaSpine_1) +  abs(datRenderHuntEvent$DThetaSpine_2) + abs(datRenderHuntEvent$DThetaSpine_3) + abs(datRenderHuntEvent$DThetaSpine_4) + abs(datRenderHuntEvent$DThetaSpine_5) + abs(datRenderHuntEvent$DThetaSpine_6) + abs(datRenderHuntEvent$DThetaSpine_7)
vTailDispFilt <- filtfilt(bf_tail, vTailDisp)

corr_speedVsTail <- ccf(vTailDispFilt,dEventSpeed_smooth,type="correlation",plot=TRUE)

plot(vTailDispFilt,type="l")
lines(dEventSpeed_smooth*100,type='l',col="blue")
#lines(vTailDir,type='l',col="green")
##END OF CURVATURE ##

##Plot Tail Segments Displacements #
X11()
plot(datRenderHuntEvent$DThetaSpine_1,type='l',col=r[1])
lines(datRenderHuntEvent$DThetaSpine_2,type='l',col=r[2])
lines(datRenderHuntEvent$DThetaSpine_3,type='l',col=r[3])
lines(datRenderHuntEvent$DThetaSpine_4,type='l',col=r[4])
lines(datRenderHuntEvent$DThetaSpine_5,type='l',col=r[5])
lines(datRenderHuntEvent$DThetaSpine_6,type='l',col=r[6])
lines(datRenderHuntEvent$DThetaSpine_7,type='l',col=r[7])


##Plot Displacement and Speed(Scaled)
X11()
plot(cumsum(vEventPathLength)) ##PLot Total Displacemnt over time
lines(dEventSpeed_smooth*100,type='l',col="blue")
lines(vTailDispFilt,type='l',col="magenta")
points(MoveboutsIdx,dEventSpeed_smooth[MoveboutsIdx]*100,col="black")
points(MoveboutsIdx_cleaned,dEventSpeed_smooth[MoveboutsIdx_cleaned]*100,col="red")
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

