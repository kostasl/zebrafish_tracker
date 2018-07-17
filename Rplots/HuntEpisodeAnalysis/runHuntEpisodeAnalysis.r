
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


idxH <- 10
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

##PlayBack
#renderHuntEventPlayback(datRenderHuntEvent,speed=1) #saveToFolder =  strFolderName

#### PROCESS BOUTS ###
vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
vDeltaDisplacement   <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2)*DIM_MMPERPX ## Path Length Calculated As Total Displacement
#nNumberOfBouts       <- 
dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
dframe               <- dframe[dframe > 0]/Fs ##Clear Any possible Nan - and Convert To Time sec  
vEventSpeed          <- meanf(vDeltaDisplacement/dframe,3) ##Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
vEventPathLength     <- cumsum(vEventSpeed)
vDistToPrey          <- meanf(sqrt( (datRenderHuntEvent$Prey_X -datRenderHuntEvent$posX )^2 + (datRenderHuntEvent$Prey_Y -datRenderHuntEvent$posY)^2   ),3)
vSpeedToPrey         <- diff(vDistToPrey,lag=1,differences=1)

#speed_Smoothed <- meanf(vEventSpeed,10)
##Replace NA with 0s
vEventSpeed[is.na(vEventSpeed)] = 0
vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed) #meanf(vEventSpeed,100) #
vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
MoveboutsIdx <- detectMotionBouts(vEventSpeed_smooth)##find_peaks(vEventSpeed_smooth*100,25)
#MoveboutsIdx_cleaned <- MoveboutsIdx[which(vEventSpeed_smooth[MoveboutsIdx] > sd(vEventSpeed_smooth[MoveboutsIdx])/3 
                                           #& vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
##Having Identified Where Bouts Are We can Extend to where we believe the motion really begins - Here Unchanged, the GauusianMix Cluster Looks ok on bout Width
MoveboutsIdx_cleaned <-MoveboutsIdx #[which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]

##Binarize , Use indicator function 1/0 for frames where Motion Occurs
#vMotionBout <- vEventSpeed_smooth
vMotionBout[ 1:NROW(vMotionBout) ]   <- 0
vMotionBout[ MoveboutsIdx_cleaned  ] <- 1
vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
vMotionBout_rle <- rle(vMotionBout)

##x10 and Round so as to detect zeroCrossings simply
vEventAccell_smooth <- round(abs(diff(vEventSpeed_smooth,lag=1,difference = 1))*10)
vEventAccell_smooth_Onset <- which(round(vEventAccell_smooth) == 0) ##Where Speed Rises Begin

##Bout On Points Are Found At the OnSet Of the Rise/ inflexion Point - Look for Previous derivative /Accelleration change
vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1


##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] 
iPairs <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event
##Remove The Motion Regions Where A Peak Was not detected / Only Keep The Bouts with Peaks
vMotionBout[1:length(vMotionBout)] = 0 ##Reset / Remove All Identified Movement
for (i in 1:iPairs)
{
  ###Motion Interval belongs to a detect bout(peak)  // Set Frame Indicators vMotionBout To Show Bout Frames
  if (any( MoveboutsIdx_cleaned >= vMotionBout_On[i] & MoveboutsIdx_cleaned < vMotionBout_Off[i] ) == TRUE)
  { 
    ##Fix Bout Onset Using Accelleration To Detect When Bout Actually Began
    ##Find Closest Speed Onset
    ##Calculate TimeDiff Between Detected BoutOnset And Actual Accelleration Onsets - Find the Onset Preceding the Detected Bout 
    OnSetTD <- vMotionBout_On[i] - vEventAccell_smooth_Onset
    ##Shift To Correct Onset Of Speed Increase / Denoting Where Bout Actually Began ##FIX ONSETS 
    if (NROW( min(OnSetTD[OnSetTD > 0  ]) ) > 0) ##If Start Of Accellaration For this Bout Can Be Found / Fix It otherwise Leave it alone
      vMotionBout_On[i] <-  vMotionBout_On[i] - min(OnSetTD[OnSetTD > 0  ]) 
    ##FIX OFFSET to when Decellaration Ends and A new One Begins
    OffSetTD <- vEventAccell_smooth_Onset - vMotionBout_Off[i]  
    if (NROW(OffSetTD[OffSetTD > 0  ]) > 0) ##If An Offset Can Be Found (Last Bout Maybe Runs Beyond Tracking Record)
      vMotionBout_Off[i] <-  vMotionBout_Off[i] + min(OffSetTD[OffSetTD > 0  ]) ##Shift |Forward To The End Of The bout
    
      
    
    vMotionBout[vMotionBout_On[i]:vMotionBout_Off[i] ] = 1 
  }
  else
  {##Remove the Ones That Do not Have a peak In them
    vMotionBout_On[i] = NA 
    vMotionBout_Off[i] = NA
  }
  
  
}

##Get Bout Statistics #### NOt Used / Replaced##
#vMotionBoutDuration_msec <- vMotionBout_Off[1:iPairs]-vMotionBout_On[1:iPairs]
#vMotionBoutDuration_msec <- 1000*vMotionBoutDuration_msec[!is.na(vMotionBoutDuration_msec)]/Fs
#vMotionBoutIntervals_msec <- 1000*(vMotionBout_On[3:(iPairs)] - vMotionBout_Off[2:(iPairs-1)])/Fs
############################

## Distance To Prey Handling  -- Fixing missing Values By Interpolation ##
## Interpolate Missing Values from Fish Speed - Assume Fish Is moving to Prey ##
##Estimate Initial DIstance From Prey Onto Which We Add the integral of Speed, By Looking At Initial PreyDist and adding any fish displacemnt to this in case The initial dist Record Is NA
vDisplacementToPrey <- (cumsum(vSpeedToPrey[!is.na(vDistToPrey)]) ) ##But diff and integration Caused a shift
vDisplacementToPrey[3:NROW(vDisplacementToPrey)] <- vDisplacementToPrey[1:(NROW(vDisplacementToPrey)-3)] ##Fix Time Shift
InitDistance             <- max(vDistToPrey[!is.na(vDistToPrey)]-vDisplacementToPrey,na.rm = TRUE )  ##vDistToPrey[!is.na(vDistToPrey)][1] + sum(vEventSpeed_smooth[(1:which(!is.na(vDistToPrey))[1])])
vSpeedToPrey[is.na(vSpeedToPrey)] <- -vEventSpeed_smooth[is.na(vSpeedToPrey)] ##Complete The Missing Speed Record To Prey By Using ThE fish Speed as estimate
vDistToPrey_Fixed <- abs(InitDistance + (cumsum(vSpeedToPrey))) ## From Initial Distance Integrate the Displacents / need -Ve Convert To Increasing Distance
vMotionBoutDistanceToPrey_mm <- vDistToPrey_Fixed[vMotionBout_On]*DIM_MMPERPX

X11()
plot((cumsum(vSpeedToPrey)+InitDistance),type='l')
lines(vDistToPrey,type='l',col="blue")


## Get Bout Statistics Again Now Using Run Length Encoding Method 
## Take InterBoutIntervals in msec from Last to first - 
vMotionBout_rle <- rle(vMotionBout)
lastBout <- max(which(vMotionBout_rle$values == 1))
firstBout <- min(which(vMotionBout_rle$values[2:lastBout] == 1)+1) ##Skip If Recording Starts With Bout , And Catch The One After the First Pause
vMotionBoutIBI <-1000*vMotionBout_rle$lengths[seq(lastBout-1,1,-2 )]/Fs #' IN msec and in reverse Order From Prey Capture Backwards
vMotionBoutDuration <-1000*vMotionBout_rle$lengths[seq(lastBout,2,-2 )]/Fs
boutSeq <- seq(NROW(vMotionBoutDistanceToPrey_mm),1,-1 )


vMotionBoutDistanceToPrey_mm <- vMotionBoutDistanceToPrey_mm[boutSeq] ##Reverse Order
stopifnot(vMotionBout_rle$values[NROW(vMotionBout_rle$lengths)] == 0 )
stopifnot(vMotionBout_rle$values[firstBout+1] == 0 ) ##THe INitial vMotionBoutIBI Is not Actually A pause interval , but belongs to motion!
datMotionBout <- cbind(boutSeq,vMotionBoutIBI,vMotionBoutDuration,vMotionBoutDistanceToPrey_mm) ##Make Data Frame
##On Bout Lengths
##Where t=0 is the capture bout, -1 -2 are the steps leading to it



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


##Plot Displacement and Speed(Scaled)
X11()
plot(cumsum(vEventPathLength),ylab="mm/sec",ylim=c(0,max(vEventSpeed_smooth))) ##PLot Total Displacemnt over time
lines(vEventSpeed_smooth,type='l',col="blue")
lines(vTailDispFilt,type='l',col="magenta")
points(MoveboutsIdx,vEventSpeed_smooth[MoveboutsIdx],col="black")
points(MoveboutsIdx_cleaned,vEventSpeed_smooth[MoveboutsIdx_cleaned],col="red")
points(vMotionBout_On,vEventSpeed_smooth[vMotionBout_On],col="red",pch=17)
points(vMotionBout_Off,vEventSpeed_smooth[vMotionBout_Off],col="blue",pch=6)
lines(vDistToPrey_Fixed*DIM_MMPERPX,col="purple",lw=2)
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

