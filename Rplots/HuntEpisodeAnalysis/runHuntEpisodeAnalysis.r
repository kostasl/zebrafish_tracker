
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
r <- c(rfc(11),"#FF0000");


#
############# Analysis AND REPLAY OF HUNT EVENTS ####
load(strDataFileName)
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
##Make an Updated list of ReTracked Hunt Events that have been imported
# datTrackedEventsRegister <- data.frame(unique(cbind(datHuntEventMergedFrames$expID,datHuntEventMergedFrames$eventID,datHuntEventMergedFrames$trackID) ))

## Setup Filters ## Can Check Bands with freqz(bf_speed) ## These are used in filterEyeTailNoise 
Fs <- 430; #sampling rate
bf_tail <- butter(1, c(0.01,0.3),type="pass"); ##Remove DC
bf_tailClass <- butter(4, c(0.01,0.35),type="pass"); ##Remove DC
bf_tailClass2 <- butter(4, 0.05,type="low"); ##Remove DC
bf_eyes <- butter(4, 0.35,type="low",plane="z");
bf_speed <- butter(4, 0.06,type="low");  ##Focus On Low Fq to improve Detection Of Bout Motion and not little Jitter motion
###
nEyeFilterWidth <- nFrWidth*6 ##For Median Filtering

if (!exists("lMotionBoutDat" ,envir = globalenv(),mode="list"))
  lMotionBoutDat <<- list() ##Declared In Global Env


#idxH <- 20
idTo <- 12#NROW(datTrackedEventsRegister)

idxDLSet <- which(datTrackedEventsRegister$groupID == "DL")
idxNLSet <- which(datTrackedEventsRegister$groupID == "NL")
idxLLSet <- which(datTrackedEventsRegister$groupID == "LL")
idxTestSet = 28#(1:NROW(datTrackedEventsRegister))


for (idxH in idxDLSet)#NROW(datTrackedEventsRegister)
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
  
  
  ## PLAYBACK ####
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
  vEventSpeed          <- meanf(vDeltaDisplacement/dframe,5) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
  

  vDeltaBodyAngle      <- diffPolar(datRenderHuntEvent$BodyAngle) #(  ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX,datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
  vTurnSpeed           <- meanf(vDeltaBodyAngle[1:NROW(dframe)]/dframe,5)
  vAngleDisplacement   <- cumsum(vDeltaBodyAngle)

    
  #vEventPathLength     <- cumsum(vEventSpeed) ##Noise Adds to Length
  vDistToPrey          <- meanf(sqrt( (datFishMotionVsTargetPrey$Prey_X -datFishMotionVsTargetPrey$posX )^2 + (datFishMotionVsTargetPrey$Prey_Y - datFishMotionVsTargetPrey$posY)^2   ),3)
  vSpeedToPrey         <- diff(vDistToPrey,lag=1,differences=1)
  lAngleToPrey <- calcRelativeAngleToPrey(datRenderHuntEvent)
  
  ##Calc Angle To Prey Per Bout
  vAngleToPrey <- data.frame(lAngleToPrey[as.character(selectedPreyID)])
  names(vAngleToPrey) = c("frameN","AngleToPrey")
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
  ## which(datFishMotionVsTargetPrey$LEyeAngle > G_THRESHUNTANGLE) , which(datFishMotionVsTargetPrey$REyeAngle < -G_THRESHUNTANGLE) )) -50,
  regionToAnalyse       <-seq(min( c( which(vEyeV > G_THRESHUNTVERGENCEANGLE) ) ) - 50,
                              min(which(vDistToPrey_Fixed_FullRange == min(vDistToPrey_Fixed_FullRange)), 
                                    max(which(vEyeV > G_THRESHUNTVERGENCEANGLE) )  )+100
                              ) ##Set To Up To The Minimum Distance From Prey
  vDistToPrey_Fixed      <- interpolateDistToPrey(vDistToPrey_Fixed_FullRange,vEventSpeed_smooth,regionToAnalyse)
  
  #plot(vTailDisp,type='l')
  ## Do Wavelet analysis Of Tail End-Edge Motion Displacements - 
  # Returns List Structure will all Relevant Data including Fq Mode Per Time Unit
  lwlt <- getPowerSpectrumInTime(vTailDisp,Fs)
  
  #MoveboutsIdx <- detectMotionBouts(vEventSpeed)##find_peaks(vEventSpeed_smooth*100,25)
  #### Cluster Tail Motion Wtih Fish Speed - Such As to Identify Motion Bouts Idx 
  #vMotionSpeed <- vEventSpeed_smooth + vTurnSpeed
  #MoveboutsIdx <- detectMotionBouts2(vEventSpeed_smooth,lwlt$freqMode)

  ########## BOUT DETECTION #################
  TurnboutsIdx <- NA
  MoveboutsIdx <- NA
  TailboutsIdx <- NA
  
  MoveboutsIdx <- detectMotionBouts(vEventSpeed_smooth)
  TailboutsIdx <- detectTailBouts(lwlt$freqMode)
  TurnboutsIdx <- detectTurnBouts(abs(vTurnSpeed),lwlt$freqMode)
  
  MoveboutsIdx  <- c(TailboutsIdx, MoveboutsIdx,TurnboutsIdx )
  ##Score Detected Frames On Overlapping Detectors
  tblMoveboutsScore<- table(MoveboutsIdx[!is.na(MoveboutsIdx)])
  
  MoveboutsIdx_cleaned <- as.numeric(names(tblMoveboutsScore[tblMoveboutsScore>1]))
  
  #######################
  
  
  #  plot(vTailDispFilt,type='l')
  #  points(which(vTailActivity==1),vTailDispFilt[which(vTailActivity==1)])
  
  #points(vTailActivity)
  ##Distance To PRey
  ##Length Of Vector Determines Analysis Range For Motion Bout 
  #MoveboutsIdx_cleaned <- which(vTailActivity==1)

  ###PLot Event Detection Summary
  #
  pdf(paste(strPlotExportPath,"/MotionBoutPage",idxH,"_exp",expID,"_event",eventID,"_track",trackID,".pdf",sep=""),width = 8,height = 12 ,paper = "a4",onefile = TRUE );
  #X11()
  par(mar=c(4,4,1.5,1.5))
  
  layout(matrix(c(1,6,2,6,3,7,4,7,5,8), 5, 2, byrow = TRUE))
    t <- seq(1:NROW(vEventSpeed_smooth))/(Fs/1000) ##Time Vector
  
    lMotionBoutDat[[idxH]]  <- calcMotionBoutInfo2(MoveboutsIdx_cleaned,vEventSpeed_smooth,vDistToPrey_Fixed_FullRange,vAngleToPrey,vTailDisp,regionToAnalyse,plotRes = TRUE)
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
    
    n<-0
    for (vAToPrey in lAngleToPrey)
    {
      l <- min(NROW(t),NROW(vAToPrey))
      n<-n+1; lines((vAToPrey[1:l,1]-min(datRenderHuntEvent$frameN))/(Fs/1000),vAToPrey[1:l,2],type='l',col=colR[n],xlab=NA,ylab=NA)
    }
    legend(max(t)-720,55,c(paste("(mm) Prey",selectedPreyID),"(Deg) R Eye","(Deg) L Eye",paste("(Deg) Prey",names(vAToPrey)) ) ,
           fill=c("purple","red","blue",colR),cex=0.7,box.lwd =0 )
    ###
    plotTailPowerSpectrumInTime(lwlt)
    polarPlotAngleToPreyVsDistance(datPlaybackHuntEvent)
    polarPlotAngleToPrey(datPlaybackHuntEvent)
    plotTailSpectrum(vTailDisp)##Tail Spectrum
    
    
  dev.off() 
  ##END OF PLOT
  
  ##Tail Fq Mode
  #X11()
  #plot(1000*1:NROW(lwlt$freqMode)/lwlt$Fs,lwlt$freqMode,type='l',ylim=c(0,50),xlab="msec",ylab="Hz",main="Tail Beat Fq Mode")
  
  
  ##Exclude Idx of Bouts for Which We do not have an angle -Make Vectors In the Right Sequence 
  BoutOnsetWithinRange <- lMotionBoutDat[[idxH]][,"vMotionBout_On"][ lMotionBoutDat[[idxH]][,"vMotionBout_On"] < NROW(vAngleToPrey ) ][lMotionBoutDat[[idxH]][,"boutSeq"]]
  BoutOffsetWithinRange <- lMotionBoutDat[[idxH]][,"vMotionBout_Off"][ lMotionBoutDat[[idxH]][,"vMotionBout_Off"] < NROW(vAngleToPrey ) ][lMotionBoutDat[[idxH]][,"boutSeq"]]
  vAnglesAtOnset <- vAngleToPrey[BoutOnsetWithinRange ,2]
  vAnglesAtOffset <- vAngleToPrey[BoutOffsetWithinRange,2]
  
  rows <- NROW(lMotionBoutDat[[idxH]])
  lMotionBoutDat[[idxH]] <- cbind(lMotionBoutDat[[idxH]] ,
                                  OnSetAngleToPrey = vAnglesAtOnset[lMotionBoutDat[[idxH]][,"boutSeq"]], ##Reverse THe Order Of Appearance Before Col. Bind
                                  OffSetAngleToPrey = vAnglesAtOffset[lMotionBoutDat[[idxH]][,"boutSeq"]],
                                  RegistarIdx = as.numeric(rep(idxH,rows)),
                                  expID=as.numeric(rep(expID,rows)),
                                  eventID=as.numeric(rep(eventID,rows)),
                                  groupID=rep(as.character(groupID) ,rows),
                                  PreyCount = rep(NROW(tblPreyRecord),rows))
} ###END OF EACH Hunt Episode Loop 

datEpisodeMotionBout <- lMotionBoutDat[[1]]
##On Bout Lengths
##Where Seq is the order Of Occurance, While Rank denotes our custom Ordering From Captcha backwards

# ###  AnALYSIS #XXX
##Make Vector Of Number oF Bouts Vs Distance
lBoutInfoPerEvent <- list()
for (rec in lMotionBoutDat)
{
  if (is.null(rec)) next;
  ##Take Distance of the 1st bout Detected (which has the largest #Rank (1 First Bout, N Last/Capture Bout))
  lBoutInfoPerEvent[[rec[1,"RegistarIdx"]]] <- list(nBouts=as.numeric(max(rec[,"boutSeq"])),
                                                       Distance= as.numeric(rec[rec[,"boutRank"] == max(rec[,"boutRank"]),"vMotionBoutDistanceToPrey_mm"]),
                                                       Angle= as.numeric(rec[rec[,"boutRank"] == max(rec[,"boutRank"]),"OnSetAngleToPrey"]), ##
                                                       groupID = as.character(datTrackedEventsRegister[ as.numeric(rec[1,"RegistarIdx"]),"groupID" ]),
                                                       Duration= as.numeric(rec[rec[,"boutRank"] == min(rec[,"boutRank"]),"vMotionBout_Off"]) - as.numeric(rec[rec[,"boutRank"] == max(rec[,"boutRank"]),"vMotionBout_On"])
                                                       )
}

datBoutVsPreyDistance <-  data.frame( do.call(rbind,lBoutInfoPerEvent ) )
datBoutVsPreyDistance[datBoutVsPreyDistance$groupID == as.character(groupID),] ##Select Group For Analysis / Plotting

###Distance ColourRing 
##Calculate Colour Idx For Each Of the Distances - Based on ncolBands
vUniqDist <- unique(round(unlist(datBoutVsPreyDistance$Distance)*10))
ncolBands <- 15
vcolIdx <-  vector() ; ##Distances Rounded / Indexed 
vdistToPrey <- round(unlist(datBoutVsPreyDistance$Distance)*10)
maxDistanceToPrey <- 50
vdistToPrey[vdistToPrey>maxDistanceToPrey] <- maxDistanceToPrey
vcolBands <- seq(min(vUniqDist),maxDistanceToPrey,length.out = ncolBands) 
for (j in 1:NROW(unlist(datBoutVsPreyDistance$Distance))) 
  vcolIdx[j] <- min(which(vcolBands >=  vdistToPrey[j] ))

# ) 
######Plot Turn Angle Vs Bearing - With DISTANCE Colour Code
Polarrfc <- colorRampPalette(rev(brewer.pal(15,'Spectral')));
colR <- c(Polarrfc( ncolBands ));
####### PLOT Turning Bout Vs Bearing TO Prey - Does the animal estimate turn amount Well?



#X11()
pdf(file= paste(strPlotExportPath,"/DistanceVsBoutCount_",groupID,".pdf",sep=""))
plot(datBoutVsPreyDistance$nBouts,datBoutVsPreyDistance$Distance,
     main = paste("Initial distance to Prey Vs Bouts Performed",unique(datBoutVsPreyDistance$groupID)  ) ,
     ylab="Distance to Prey  (mm)",
     xlab="Number of Tracking Movements",
     ylim=c(0,6),
     xlim=c(0,max(unlist(datBoutVsPreyDistance$nBouts) )),
     col=colR[vcolIdx],
     pch=19
     )
dev.off()


pdf(file= paste(strPlotExportPath,"/BearingVsBoutCount_",groupID,".pdf",sep=""))
plot(datBoutVsPreyDistance$nBouts,datBoutVsPreyDistance$Angle,
     main = paste("Initial Bearing to Prey Vs Bouts Performed",unique(datBoutVsPreyDistance$groupID)  ) ,
     ylab="Angle to Prey  (mm)",
     xlab="Number of Tracking Movements",
     ylim=c(-180,180),xlim=c(0,max(unlist(datBoutVsPreyDistance$nBouts) )),
     col=colR[vcolIdx],
     pch=19)

dev.off()

#X11()
pdf(file= paste(strPlotExportPath,"/DurationVsBoutCount_",groupID,".pdf",sep=""))
plot(unlist(datBoutVsPreyDistance$nBouts),1000*unlist(datBoutVsPreyDistance$Duration)/Fs,
     main = paste("Duration Of Hunt Event Vs Bouts Performed ",unique(datBoutVsPreyDistance$groupID)  )  ,
     ylab="Time  (msec)",
     xlab="Number of Tracking Movements",
     ylim=c(0,5000),xlim=c(0,max(unlist(datBoutVsPreyDistance$nBouts) )),
     col=colR[vcolIdx],
     pch=19)
dev.off()

pdf(file= paste(strPlotExportPath,"/DurationVsDistance_",groupID,".pdf",sep=""))
plot(unlist(datBoutVsPreyDistance$Distance),1000*unlist(datBoutVsPreyDistance$Duration)/Fs,
     main = paste("Duration Of Hunt Event Vs Distance ",unique(datBoutVsPreyDistance$groupID) )  ,
     ylab="Time  (msec)",
     xlab="Distance (mm)",
     ylim=c(0,5000),xlim=c(0,6),
     col=colR[vcolIdx],
     pch=19)

dev.off()


### Box Plots Per Bout ##
####Select Subset Of Data To Analyse
datMotionBoutCombinedAll <-  data.frame( do.call(rbind,lMotionBoutDat ) )
#datMotionBoutCombinedAll$groupID <- levels(datTrackedEventsRegister$groupID)[datMotionBoutCombinedAll$groupID]
datMotionBoutCombined <-datMotionBoutCombinedAll#
#datMotionBoutCombinedAll[datMotionBoutCombinedAll$groupID == "DL", ] 
##Turns Vs Bearing To Prey


####### PLOT Turning Bout Vs Bearing TO Prey - Does the animal estimate turn amount Well?
##Add Angle To Prey OnTop Of Eye Angles##
Polarrfc <- colorRampPalette(rev(brewer.pal(12,'Blues')));
colR <- c("#000000",Polarrfc(max(datMotionBoutCombinedAll$boutRank) ) ,"#FF0000");
ncolBands <- NROW(colR)
#X11()
pdf(file= paste(strPlotExportPath,"/BoutTurnsToPreyWithBoutSeq_",groupID,".pdf",sep=""))
plot(datMotionBoutCombinedAll$OnSetAngleToPrey,datMotionBoutCombinedAll$OnSetAngleToPrey-datMotionBoutCombinedAll$OffSetAngleToPrey,
     main=paste("Turn Size Vs Bearing To Prey ",levels(groupID)[unique(datMotionBoutCombinedAll$groupID)], "+ Bout Number" ),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-90,90),
     ylim=c(-90,90),col=colR[datMotionBoutCombinedAll$boutSeq] ,pch=19) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Punctuate 1st Turn To Prey
datFirstBoutPoints <- cbind(boutSeq = datMotionBoutCombinedAll[datMotionBoutCombinedAll$boutSeq == 1,]$OnSetAngleToPrey,Turn= datMotionBoutCombinedAll[datMotionBoutCombinedAll$boutSeq == 1,]$OnSetAngleToPrey-datMotionBoutCombinedAll[datMotionBoutCombinedAll$boutSeq == 1,]$OffSetAngleToPrey
                            ,RegistarIdx=datMotionBoutCombinedAll[datMotionBoutCombinedAll$boutSeq == 1,]$RegistarIdx)
points(datFirstBoutPoints[,1],
       datFirstBoutPoints[,2],
       pch=9)
points(1:ncolBands,rep(-90, ncolBands), col=colR[1:ncolBands],pch=15) ##Add Legend Head Map
text(-3,-87,labels = paste(min(datMotionBoutCombinedAll$boutRank) ,"#" )  ) ##Heatmap range min
text(ncolBands+4,-87,labels = paste(max(datMotionBoutCombinedAll$boutRank) ,"#" )  )
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=2);
## Plot arrows showing Bout Turns Connecting Bouts From Same Experiment
#for (expID in unique(datMotionBoutCombinedAll$expID) ) 
#{
#  datExp <- datMotionBoutCombinedAll[datMotionBoutCombinedAll$expID ==expID, ]
#  for (ii in max(datExp$boutSeq):2) ##Draw Arrows Showing Sequence Of Bouts
#    arrows(x0=datExp[ii,]$OnSetAngleToPrey, y0= datExp[ii,]$OnSetAngleToPrey-datExp[ii,]$OffSetAngleToPrey ,
#          x1=datExp[ii-1,]$OnSetAngleToPrey, y1=datExp[ii-1,]$OnSetAngleToPrey-datExp[ii-1,]$OffSetAngleToPrey,
#          length = 0.1)
#}
#####
dev.off()

##Calculate Colour Idx For Each Of the Distances - Based on ncolBands
vUniqDist <- unique(round(datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm*10))
ncolBands <- 15
vcolIdx <-  vector() ; ##Distances Rounded / Indexed 
vdistToPrey <- round(datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm*10)
maxDistanceToPrey <- 50
vdistToPrey[vdistToPrey>maxDistanceToPrey] <- maxDistanceToPrey
vcolBands <- seq(min(vUniqDist),maxDistanceToPrey,length.out = ncolBands) 
for (j in 1:NROW(datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm)) 
  vcolIdx[j] <- min(which(vcolBands >  vdistToPrey[j] ))

# ) 
######Plot Turn Angle Vs Bearing - With DISTANCE Colour Code
Polarrfc <- colorRampPalette(rev(brewer.pal(15,'Spectral')));
colR <- c(Polarrfc( ncolBands ));
####### PLOT Turning Bout Vs Bearing TO Prey - Does the animal estimate turn amount Well?

#X11()
pdf(file= paste(strPlotExportPath,"/BoutTurnsToPreyWithPreyDist_",groupID,".pdf",sep=""))
plot(datMotionBoutCombinedAll$OnSetAngleToPrey,datMotionBoutCombinedAll$OnSetAngleToPrey-datMotionBoutCombinedAll$OffSetAngleToPrey,
     main=paste("Turn Size Vs Bearing To Prey ",levels(groupID)[unique(datMotionBoutCombinedAll$groupID)]," + Distance" ),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-90,90),
     ylim=c(-90,90),col=colR[vcolIdx] ,pch=19) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
points(1:ncolBands,rep(-90, ncolBands), col=colR[1:ncolBands],pch=15) ##Add Legend Head Map
text(-3,-57,labels = paste(min(vUniqDist)/10,"mm" )  ) ##Heatmap range min
text(ncolBands+4,-57,labels = paste(maxDistanceToPrey/10,"mm" )  )
points(datMotionBoutCombinedAll$OnSetAngleToPrey,datMotionBoutCombinedAll$OnSetAngleToPrey-datMotionBoutCombinedAll$OffSetAngleToPrey,
     main=paste("Turn Size Vs Bearing To Prey G",levels(groupID)[unique(datMotionBoutCombinedAll$groupID)] ),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-90,90),
     ylim=c(-90,90),col="Black" ,pch=16,cex=0.2) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,  90,0); segments(-90,-90,90,90,lwd=2);

#dev.copy(pdf,file= paste(strPlotExportPath,"/BoutTurnsToPreyWithPreyDist_",groupID,".pdf",sep=""))
dev.off()
################# #### # ##  ## #


X11()
plot(datMotionBoutCombined$boutRank,datMotionBoutCombined$vMotionBoutDistanceToPrey_mm,main="Distance From Prey",ylab="mm")
boxplot(as.numeric(datMotionBoutCombined$vMotionBoutDistanceToPrey_mm) ~ as.numeric(datMotionBoutCombined$boutRank),main="Distance From Prey",ylab="mm",xlab="Bout Sequence (From Capture - Backwards)")


X11()
boxplot(as.numeric(datMotionBoutCombined$OnSetAngleToPrey) ~ as.numeric(datMotionBoutCombined$boutRank),
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


