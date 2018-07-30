
### Hunt Episode Data Extraction and Analysis Scripts 
### Kostasl Jul 2018 
## It is aimed at identified hunt event that have been carefully retracked and imported separatelly from these new Tracked CsV files using runImportHuntEpisodeTrackFiles
## It provides the data for identifying Bouts using FishSpeed. 
## Notes: 
## * Gaussian Mixt. Clusters are used to Identify Bout from Speeds, however the BoutOnset-Offset Is then Fixed To Identify When Speed Ramps Up , 
## and final decelpartion is used to obtain  more accurate estimates of BoutDurations-
## * The Distance to Prey Is calculated But Also Interpolated using Fish Displacement where there are missing values - 
## \Notes:
## Can use findLabelledEvent( datTrackedEventsRegister[IDXOFHUNT,]) to locate which HuntEvent is associated from the Labelled Set Record, And Retrack it by running 
## main_LabellingBlind.r and providing the row.name as ID 
#####


library(signal)
library(MASS)
library(mclust,quietly = TRUE)

require(Rwave)

source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r")
source("TrackerDataFilesImport_lib.r")
source("plotTrackScatterAndDensities.r")

strDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis",".RData",sep="") ##To Which To Save After Loading
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strDataFileName))

#    for (i in 1:40) dev.off()


rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
r <- c(rfc(8),"#FF0000");


#
############# Analysis AND REPLAY OF HUNT EVENTS ####
load(strDataFileName)
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
##Make an Updated list of ReTracked Hunt Events that have been imported
# datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))

## Setup Filters ## Can Check Bands with freqz(bf_speed) ## These are used in filterEyeTailNoise 
Fs <- 430; #sampling rate
bf_tail <- butter(1, c(0.01,0.3),type="pass"); ##Remove DC
bf_tailClass <- butter(4, c(0.01,0.3),type="pass"); ##Remove DC
bf_tailClass2 <- butter(4, 0.05,type="low"); ##Remove DC
bf_eyes <- butter(4, 0.025,type="low",plane="z");
bf_speed <- butter(4, 0.04,type="low");  ##Focus On Low Fq to improve Detection Of Bout Motion and not little Jitter motion
###
nEyeFilterWidth <- nFrWidth*8 ##For Median Filtering

lMotionBoutDat <- list()


#idxH <- 20
idTo <- 12#NROW(datTrackedEventsRegister)


idxNLSet <- which(datTrackedEventsRegister$groupID == "NL")
idxLLSet <- which(datTrackedEventsRegister$groupID == "LL")
idxTestSet = c(14)


for (idxH in idxTestSet)#NROW(datTrackedEventsRegister)
{
  
  expID <- datTrackedEventsRegister[idxH,]$expID
  trackID<- datTrackedEventsRegister[idxH,]$trackID
  eventID <- datTrackedEventsRegister[idxH,]$eventID
  groupID <- datTrackedEventsRegister[idxH,]$groupID
  selectedPreyID <- datTrackedEventsRegister[idxH,]$PreyIDTarget
  
  message(paste(idxH, ".Process Hunt Event Expid:",expID,"Event:",eventID))
  
  datPlaybackHuntEvent <- datHuntEventMergedFrames[datHuntEventMergedFrames$expID==expID 
                                                 & datHuntEventMergedFrames$trackID==trackID 
                                                 & datHuntEventMergedFrames$eventID==eventID,]
  
  
  
  
  
  strFolderName <- paste( strPlotExportPath,"/renderedHuntEvent",expID,"_event",eventID,"_track",trackID,sep="" )
  
  ############ PREY SELECTION #####
  ## Begin Data EXtraction ###
  ## Remove Multiple Prey Targets ########
  ##Get Number oF Records per Prey
  tblPreyRecord <-table(datPlaybackHuntEvent$PreyID) 
  if (NROW(tblPreyRecord) > 1)
    warning("Multiple Prey Items Tracked In Hunt Episode-Selecting Longest Track")
  
  ##### FILTERS #######
  message("Filtering Fish Motion (on all PreyID replicates)...")
  ldatFish <- list()
  for (p in names(tblPreyRecord))
  {
    message(p)
    ldatFish[[as.character(p)]] <- filterEyeTailNoise(datPlaybackHuntEvent[!is.na(datPlaybackHuntEvent$PreyID) 
                                                                           & datPlaybackHuntEvent$PreyID == p,])
  }
  ldatFish[["NA"]] <- filterEyeTailNoise(datPlaybackHuntEvent[is.na(datPlaybackHuntEvent$PreyID) ,])
  
  ##Recombine the datframes Split By PreyID
  datPlaybackHuntEvent <- do.call(rbind,ldatFish)
  
  
  ##PLAYBACK ####
  #     renderHuntEventPlayback(datPlaybackHuntEvent,selectedPreyID,speed=1) #saveToFolder =  strFolderName
  ########################
  
  
  ##Add PreyTarget ID To Register Save The Prey Target To Register
  if (!any(names(datTrackedEventsRegister) == "PreyIDTarget"))
    datTrackedEventsRegister$PreyIDTarget <- NA
  if (!any(names(datTrackedEventsRegister) == "startFrame"))
    datTrackedEventsRegister$startFrame <- NA
  
  
  #selectedPreyID <- max(as.numeric(names(which(tblPreyRecord == max(tblPreyRecord)))))
  ##Check If Assigned OtherWise Automatically Select the longest Track
  if (is.na(datTrackedEventsRegister[idxH,]$PreyIDTarget)) 
  {
    selectedPreyID <-  max(as.numeric(names(which(tblPreyRecord == max(tblPreyRecord)))))
    datTrackedEventsRegister[idxH,]$PreyIDTarget <- selectedPreyID
    datTrackedEventsRegister[idxH,]$PreyCount    <- NROW(tblPreyRecord)
    datTrackedEventsRegister[idxH,]$startFrame   <- min(datRenderHuntEvent$frameN)
    saveRDS(datTrackedEventsRegister,file=strRegisterDataFileName) ##Save With Dataset Idx Identifier
  }
  
  #
  ##Save The Selected Prey Item
  if (selectedPreyID != datTrackedEventsRegister[idxH,]$PreyIDTarget)
  {
    message(paste("Targeted Prey Changed For Hunt Event to ID:",selectedPreyID," - Updating Register...") )
    datTrackedEventsRegister[idxH,]$PreyIDTarget <- selectedPreyID
    saveRDS(datTrackedEventsRegister,file=strRegisterDataFileName) ##Save With Dataset Idx Identifier
  }
  ################ END OF PREY SELECT / Start Processing ###
  
  
  ## Select Prey Specific Subset
  datFishMotionVsTargetPrey <- ldatFish[[as.character(selectedPreyID)]] #datPlaybackHuntEvent[datPlaybackHuntEvent$PreyID == ,] 
  datRenderHuntEvent <- datFishMotionVsTargetPrey
  
  ### Filter / Process Motion Variables - Noise Removal / 
  #datFishMotionVsTargetPrey <- filterEyeTailNoise(datFishMotionVsTargetPrey)
  ##
  ##Vector Of Vergence Angle
  vEyeV <- datFishMotionVsTargetPrey$LEyeAngle-datFishMotionVsTargetPrey$REyeAngle
  
  
    
  
  
  #### PROCESS BOUTS ###
  vDeltaXFrames        <- diff(datRenderHuntEvent$posX,lag=1,differences=1)
  vDeltaYFrames        <- diff(datRenderHuntEvent$posY,lag=1,differences=1)
  vDeltaDisplacement   <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
  
  
  #nNumberOfBouts       <- 
  dframe               <- diff(datRenderHuntEvent$frameN,lag=1,differences=1)
  dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - and Convert To Time sec  
  vEventSpeed          <- meanf(vDeltaDisplacement/dframe,3) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
  

  vDeltaBodyAngle      <- diffPolar(datRenderHuntEvent$BodyAngle) #(  ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX,datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
  vTurnSpeed           <- meanf(vDeltaBodyAngle[1:NROW(dframe)]/dframe,3)
  vAngleDisplacement   <- cumsum(vDeltaBodyAngle)

    
  #vEventPathLength     <- cumsum(vEventSpeed) ##Noise Adds to Length
  vDistToPrey          <- meanf(sqrt( (datFishMotionVsTargetPrey$Prey_X -datFishMotionVsTargetPrey$posX )^2 + (datFishMotionVsTargetPrey$Prey_Y - datFishMotionVsTargetPrey$posY)^2   ),3)
  vSpeedToPrey         <- diff(vDistToPrey,lag=1,differences=1)

  ## Tail Motion ####
  vTailDir <-  datRenderHuntEvent$DThetaSpine_1 +  datRenderHuntEvent$DThetaSpine_2 + datRenderHuntEvent$DThetaSpine_3 + datRenderHuntEvent$DThetaSpine_4 + datRenderHuntEvent$DThetaSpine_5 + datRenderHuntEvent$DThetaSpine_6 + datRenderHuntEvent$DThetaSpine_7
  vTailDisp <-  datRenderHuntEvent$DThetaSpine_6 + datRenderHuntEvent$DThetaSpine_7 #+ datRenderHuntEvent$DThetaSpine_7 #+ datRenderHuntEvent$DThetaSpine_7 #abs(datRenderHuntEvent$DThetaSpine_1) +  abs(datRenderHuntEvent$DThetaSpine_2) + abs(datRenderHuntEvent$DThetaSpine_3) + abs(datRenderHuntEvent$DThetaSpine_4) + abs(datRenderHuntEvent$DThetaSpine_5) + abs(datRenderHuntEvent$DThetaSpine_6) + abs(datRenderHuntEvent$DThetaSpine_7)
  vTailDisp <- filtfilt(bf_tailClass, clipEyeRange(vTailDisp,-120,120))
  vTailDispFilt <- filtfilt(bf_tailClass2,abs( vTailDisp) )  ##Heavily Filtered and Used For Classifying Bouts

  
  #X11()
  #plot((1000*1:NROW(vTailDisp)/Fs),vTailDisp,type='l')
  
  ##Plot Tail Spectral Density
  #png(filename=paste(strPlotExportPath,"/TailSpectrum",idxH,"_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
  
  #dev.off()
  
  #tmp<-mk.cwt(w,noctave = floor(log2(length(w)))-1,nvoice=10)
  
  #speed_Smoothed <- meanf(vEventSpeed,10)
  ##Replace NA with 0s
  vEventSpeed[is.na(vEventSpeed)] = 0
  vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed) #meanf(vEventSpeed,100) #
  vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
  vEventPathLength <- cumsum(vEventSpeed_smooth)
  
  vTurnSpeed[is.na(vTurnSpeed)] <- 0
  vTurnSpeed <- filtfilt(bf_speed, vTurnSpeed)
  
  vDistToPrey_Fixed_FullRange      <- interpolateDistToPrey(vDistToPrey[1:NROW(vEventSpeed_smooth)],vEventSpeed_smooth)
  ##Find Region Of Interest For Analysis Of Bouts
  ## As the Furthers point Between : Either The Prey Distance Is minimized, or The Eye Vergence Switches Off) 
  regionToAnalyse       <-seq(1,
                              max(which(vDistToPrey_Fixed_FullRange == min(vDistToPrey_Fixed_FullRange)), 
                                    max(which(vEyeV > G_THRESHUNTVERGENCEANGLE) )  )+150
                              ) ##Set To Up To The Minimum Distance From Prey
  vDistToPrey_Fixed      <- interpolateDistToPrey(vDistToPrey_Fixed_FullRange,vEventSpeed_smooth,regionToAnalyse)
  
  #plot(vTailDisp,type='l')
  ## Do Wavelet analysis Of Tail End-Edge Motion Displacements - 
  # Returns List Structure will all Relevant Data including Fq Mode Per Time Unit
  lwlt <- getPowerSpectrumInTime(vTailDisp,Fs)
  
  
  #MoveboutsIdx <- detectMotionBouts(vEventSpeed)##find_peaks(vEventSpeed_smooth*100,25)
  #### Cluster Tail Motion Wtih Fish Speed - Such As to Identify Motion Bouts Idx 
  #vMotionSpeed <- vEventSpeed_smooth + vTurnSpeed
  MoveboutsIdx <- detectMotionBouts2(vEventSpeed_smooth,lwlt$freqMode)
  TurnboutsIdx <- detectTurnBouts(abs(vTurnSpeed),lwlt$freqMode)
  MoveboutsIdx_cleaned <- TurnboutsIdx #c(MoveboutsIdx,TurnboutsIdx)# which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   ) #MoveboutsIdx# 
  MoveboutsIdx_cleaned[MoveboutsIdx_cleaned %in% MoveboutsIdx] <-  NA
  ##Append The MoveBoutsIdx
  MoveboutsIdx_cleaned <- c(MoveboutsIdx_cleaned[!is.na(MoveboutsIdx_cleaned)],MoveboutsIdx)
  
  
  ## Detect Tail Motion Bouts
  vTailActivity <- rep(0,NROW(vTailDispFilt))
  vTailActivity[vTailDispFilt>5] <- 1
  #vTailActivity[vTailDispFilt <=2] <- 0
  
  
  #  plot(vTailDispFilt,type='l')
  #  points(which(vTailActivity==1),vTailDispFilt[which(vTailActivity==1)])
  
  #points(vTailActivity)
  ##Distance To PRey
  ##Length Of Vector Determines Analysis Range For Motion Bout 
  #MoveboutsIdx_cleaned <- which(vTailActivity==1)

  ###PLot Event Detection Summary
  #
  #pdf(paste(strPlotExportPath,"/MotionBoutPage",idxH,"_exp",expID,"_event",eventID,"_track",trackID,".pdf",sep=""),width = 8,height = 12 ,paper = "a4",onefile = TRUE );
  X11()
  par(mar=c(4,4,1.5,1.5))
  
  layout(matrix(c(1,6,2,6,3,7,4,7,5,8), 5, 2, byrow = TRUE))
    t <- seq(1:NROW(vEventSpeed_smooth))/(Fs/1000) ##Time Vector
  
    lMotionBoutDat[[idxH]]  <- calcMotionBoutInfo2(MoveboutsIdx_cleaned,vEventSpeed_smooth,vDistToPrey_Fixed_FullRange,vTailDisp,regionToAnalyse,plotRes = TRUE)
    ##Change If Fish Heading
    plot(t,vAngleDisplacement[1:NROW(t)],type='l',
         xlab="(msec)",
         ylab="Degrees",
         col="blue",main=" Angle Displacement")
    lines(t,cumsum(vTurnSpeed)[1:NROW(t)],type='l',lwd=2,
         xlab=NA,
         ylab=NA,
         col="blue4")
    
    par(new = FALSE)
    plot(t,vDistToPrey_Fixed_FullRange[1:NROW(t)]*DIM_MMPERPX,type='l',
         xlab="(msec)",
         ylab="Distance (mm)",
         col="purple",main="Motion Relative Prey and Eye Angles",lwd=2,ylim=c(0,5))
    axis(side = 2,col="purple",cex=1.2,lwd=2)
    
    ##Add Eye Angles  ##
    par(new = TRUE )
    par(mar=c(4,4,2,2))
    plot(t,datRenderHuntEvent$REyeAngle[1:NROW(t)],axes=F,col="red3",type='l',xlab=NA,ylab=NA,cex=1.2,ylim=c(-55,55))
    axis(side = 4,col="red")
    mtext(side = 4, line = 3, 'Angles (Deg)')
    lines(t,datRenderHuntEvent$LEyeAngle[1:NROW(t)],axes=F,col="blue",type='l',xlab=NA,ylab=NA)
    
    ##Add Angle To Prey OnTop Of Eye Angles##
    Polarrfc <- colorRampPalette(rev(brewer.pal(8,'Dark2')));
    colR <- c(Polarrfc(NROW(tblPreyRecord) ) ,"#FF0000");
    lAngleToPrey <- calcRelativeAngleToPrey(datRenderHuntEvent)
    n<-0
    for (vAngleToPrey in lAngleToPrey)
    {
      l <- min(NROW(t),NROW(vAngleToPrey))
      n<-n+1; lines((vAngleToPrey[1:l,1]-min(datRenderHuntEvent$frameN))/(Fs/1000),vAngleToPrey[1:l,2],type='l',col=colR[n],xlab=NA,ylab=NA)
    }
    legend(max(t)-720,55,c(paste("(mm) Prey",selectedPreyID),"(Deg) R Eye","(Deg) L Eye",paste("(Deg) Prey",names(lAngleToPrey)) ) ,fill=c("purple","red","blue",colR),cex=0.7,box.lwd =0 )
    ###
    plotTailPowerSpectrumInTime(lwlt)
    polarPlotAngleToPreyVsDistance(datPlaybackHuntEvent)
    polarPlotAngleToPrey(datPlaybackHuntEvent)
    plotTailSpectrum(vTailDisp)##Tail Spectrum
    
    
  #dev.off() 
  ##END OF PLOT
  
  ##Tail Fq Mode
  X11()
  plot(1000*1:NROW(lwlt$freqMode)/lwlt$Fs,lwlt$freqMode,type='l',ylim=c(0,50),xlab="msec",ylab="Hz",main="Tail Beat Fq Mode")
  
  ##Calc Angle To Prey Per Bout
  vAngleToPrey <- lAngleToPrey[as.character(selectedPreyID)]
  
  ##Exclude Idx of Bouts for Which We do not have an angle
  BoutOnsetWithinRange <- lMotionBoutDat[[idxH]][,"vMotionBout_On"][ lMotionBoutDat[[idxH]][,"vMotionBout_On"] < NROW(vAngleToPrey[[1]] ) ]
  vAnglesAtOnset <- vAngleToPrey[[1]][BoutOnsetWithinRange,2]
  rows <- NROW(lMotionBoutDat[[idxH]])
  lMotionBoutDat[[idxH]] <- cbind(lMotionBoutDat[[idxH]] ,
                                  AngleToPrey = vAnglesAtOnset,
                                  RegistarIdx = as.numeric(rep(idxH,rows)),
                                  expID=as.numeric(rep(expID,rows)),
                                  eventID=as.numeric(rep(eventID,rows)),
                                  groupID=rep((groupID) ,rows),
                                  PreyCount = rep(NROW(tblPreyRecord),rows))
} ###END OF EACH Hunt Episode Loop 

datEpisodeMotionBout <- lMotionBoutDat[[1]]
##On Bout Lengths
##Where t=0 is the capture bout, -1 -2 are the steps leading to it

# 
# ## Plot Durations of Pause/Go
# X11()
# plot(datEpisodeMotionBout[,"vMotionBoutDuration"],
#      xlab="Bout",ylab="msec",xlim=c(0,NROW(datEpisodeMotionBout) ),ylim=c(0,500),
#      col="red",main="Bout Duration",pch=16) ##Take Every Bout Length
# points(datEpisodeMotionBout[,"vMotionBoutIBI"],col="blue",pch=21) ##Take every period between / Inter Bout Interval
# legend(1,400,c("Motion","Pause" ),col=c("red","blue"),pch=c(16,21) )
# 
# X11()
# plot(datEpisodeMotionBout[,"vMotionBoutDistanceToPrey_mm"],
#      xlab="Bout",ylab="mm",xlim=c(0,NROW(datEpisodeMotionBout)),ylim=c(0,3),
#      col="red",main="Bout Distance To Prey",pch=16) ##Take Every Bout Length
# 
# X11()
# plot(datEpisodeMotionBout[,"vMotionBoutDistanceTravelled_mm"],
#      xlab="Bout",ylab="mm",xlim=c(0,NROW(datEpisodeMotionBout)),ylim=c(0,2),
#      col="red",main="Bout Power",pch=16) ##Take Every Bout Length
# 
# 
# X11()
# plot(datRenderHuntEvent$frameN,datRenderHuntEvent$LEyeAngle,type='l',col="blue",ylim=c(-60,60),main="Eye Motion ")
# lines(datRenderHuntEvent$frameN,datRenderHuntEvent$REyeAngle,type='l',col="magenta")

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


# ###  AnALYSIS #XXX
##Make Vector Of Number oF Bouts Vs Distance
lBoutsVsPreyDistance <- list()
for (rec in lMotionBoutDat)
{
  if (is.null(rec)) next;
  ##Take Distance of the 1st bout Detected (which has the largest #Rank (1 Last, N first))
  lBoutsVsPreyDistance[[rec[1,"RegistarIdx"]]] <- list(nBouts=max(rec[,"boutSeq"]),
                                                       Distance= as.numeric(rec[rec[,"boutRank"] == max(rec[,"boutRank"]),"vMotionBoutDistanceToPrey_mm"]),
                                                       Angle= as.numeric(rec[rec[,"boutRank"] == max(rec[,"boutRank"]),"AngleToPrey"]))
}

datBoutVsPreyDistance <-  data.frame( do.call(rbind,lBoutsVsPreyDistance ) )
X11()
plot(datBoutVsPreyDistance$nBouts,datBoutVsPreyDistance$Distance,
     main = "Initial distance to Prey Vs Bouts Performed",
     ylab="Distance to Prey  (mm)",
     xlab="Number of Tracking Movements",
     ylim=c(0,6),xlim=c(0,max(unlist(datBoutVsPreyDistance$nBouts) )))

X11()
plot(datBoutVsPreyDistance$nBouts,datBoutVsPreyDistance$Angle,
     main = "Initial Bearing to Prey Vs Bouts Performed",
     ylab="Angle to Prey  (mm)",
     xlab="Number of Tracking Movements",
     ylim=c(-180,180),xlim=c(0,max(unlist(datBoutVsPreyDistance$nBouts) )))


### Box Plots Per Bout ##
####Select Subset Of Data To Analyse
datMotionBoutCombinedAll <-  data.frame( do.call(rbind,lMotionBoutDat ) )
#datMotionBoutCombinedAll$groupID <- levels(datTrackedEventsRegister$groupID)[datMotionBoutCombinedAll$groupID]
datMotionBoutCombined <-datMotionBoutCombinedAll#
#datMotionBoutCombinedAll[datMotionBoutCombinedAll$groupID == "DL", ] 

X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDistanceToPrey_mm,main="Distance From Prey",ylab="mm")
boxplot(as.numeric(datMotionBoutCombined$vMotionBoutDistanceToPrey_mm) ~ as.numeric(datMotionBoutCombined$boutRank),main="Distance From Prey",ylab="mm",xlab="Bout Sequence (From Capture - Backwards)")


X11()
boxplot(as.numeric(datMotionBoutCombined$AngleToPrey) ~ as.numeric(datMotionBoutCombined$boutRank),
        main="Bearing To Prey",
        ylab="(Deg)",
        xlab="Bout Sequence (From Capture - Backwards)",
        ylim=c(-180,180))


X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDistanceTravelled_mm,main="Distance Of Bout (power)",ylab="mm")
boxplot(datMotionBoutCombined$vMotionBoutDistanceTravelled_mm ~
datMotionBoutCombined$boutRank,main="Distance Of Bout (power)",ylab="mm",xlab="Bout Sequence (From Capture - Backwards)")


X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDuration,main=" Bout Duration",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")
boxplot(datMotionBoutCombined$vMotionBoutDuration ~ datMotionBoutCombined$boutRank,main=" Bout Duration",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)")

X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutIBI,main=" Inter Bout Intervals ",ylab="msec",xlab="Bout Sequence (From Capture -Backwards)")
boxplot( datMotionBoutCombined$vMotionBoutIBI ~ datMotionBoutCombined$boutRank,main=" Inter Bout Intervals ",ylab="msec",xlab="Bout Sequence (From Capture - Backwards)") 
# for (i in1:20) #dev.off()


##Plot Tail
#X11()
#plot(datRenderHuntEvent$DThetaSpine_1 ,type='l',col=r[1])
#lines(datRenderHuntEvent$DThetaSpine_2 ,type='l',col=r[2])
#lines(datRenderHuntEvent$DThetaSpine_3 ,type='l',col=r[3])
#lines(datRenderHuntEvent$DThetaSpine_4 ,type='l',col=r[4])
#lines(datRenderHuntEvent$DThetaSpine_5 ,type='l',col=r[5])
#lines(datRenderHuntEvent$DThetaSpine_6 ,type='l',col=r[6])
#lines(datRenderHuntEvent$DThetaSpine_7 ,type='l',col=r[7])


