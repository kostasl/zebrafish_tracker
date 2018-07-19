### Library For Analysis of Hunt Episodes ###

library(signal)
library(MASS)
library(mclust,quietly = TRUE)

#library(sBIC)

citation("mclust")


##Clusters Fish Speed Measurements into Bout And Non Bout
##Use 3 For Better Discrimination When  There Are Exist Bouts Of Different Size
detectMotionBouts <- function(vEventSpeed)
{
  prior_factor2 <- 0.05 ## Adds a prior shift in the threshold Of Classification
  prior_factor1 <- 1.0 ## Adds a prior shift in the threshold Of Classification 
  colClass <- c("#FF0000","#00FF22","#0000FF")
  
  #t <- datRenderHuntEvent$frameN
  
  #BIC <- mclustBIC(dEventSpeed)
  
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  fit <- Mclust(vEventSpeed ,G=3, prior =  priorControl(functionName="defaultPrior", mean=c(1,5,10),shrinkage=0.01 ) )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  summary(fit)
  
  #region <- min(NROW(t),NROW(vEventSpeed))
  X11()
  plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
  rug(vEventSpeed)
  
  #X11()
  #boutClass <- fit$classification
  #plot(dEventSpeed[1:region],type='l',col=colClass[1])
  #points(which(boutClass == 2), dEventSpeed[boutClass == 2],type='p',col=colClass[2])

  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which(fit$classification == 3 ) )
  #return (which( fit$z[,3]> fit$z[,1]*prior_factor1 | fit$z[,3]> fit$z[,2]*prior_factor2    )) #
  
}

## Distance To Prey Handling  -- Fixing missing Values By Interpolation Using Fish Motion##
##
interpolateDistToPrey <- function(vDistToPrey,vEventSpeed_smooth)
{
  ##Calc Speed - And Use it To Merge The Missing Values 
  vSpeedToPrey         <- diff(vDistToPrey,lag=1,differences=1)
  
  ## Interpolate Missing Values from Fish Speed - Assume Fish Is moving to Prey ##
  ##Estimate Initial DIstance From Prey Onto Which We Add the integral of Speed, By Looking At Initial PreyDist and adding any fish displacemnt to this in case The initial dist Record Is NA
  vDisplacementToPrey <- (cumsum(vSpeedToPrey[!is.na(vDistToPrey)]) ) ##But diff and integration Caused a shift
  vDisplacementToPrey[9:(NROW(vDisplacementToPrey))] <- vDisplacementToPrey[1:(NROW(vDisplacementToPrey)-8)] ##Fix Time Shift
  ##Compare Mean Distance Between them - Only Where Orig. PreyDist Is not NA - To obtain Integral Initial Constant (Starting Position)
  InitDistance             <- mean(vDistToPrey[!is.na(vDistToPrey)]-vDisplacementToPrey[!is.na(vDistToPrey)],na.rm = TRUE )  ##vDistToPrey[!is.na(vDistToPrey)][1] + sum(vEventSpeed_smooth[(1:which(!is.na(vDistToPrey))[1])])
  ##Add The Missing Speed Values From The fish Speed, Assumes Prey Remains Fixed And Fish Moves towards Prey
  vSpeedToPrey[is.na(vSpeedToPrey)] <- -vEventSpeed_smooth[is.na(vSpeedToPrey)] ##Complete The Missing Speed Record To Prey By Using ThE fish Speed as estimate
  vDistToPrey_Fixed <- abs(InitDistance +  vDisplacementToPrey)# (cumsum(vSpeedToPrey))) ## From Initial Distance Integrate the Displacents / need -Ve Convert To Increasing Distance
  
  
  #X11() ##Compare Estimated To Recorded Prey Distance
  #plot(vDistToPrey_Fixed,type='l')
  #lines(vDistToPrey,type='l',col="blue")
  #legend()
  
  return(vDistToPrey_Fixed)
}

##Uses The Detected Regions Of Bouts to extract data, on BoutOnset-Offset - Duration, Distance from Prey and Bout Power as a measure of distance moved during bout
## 
calcMotionBoutInfo <- function(MoveboutsIdx,vEventSpeed_smooth,vDistToPrey)
{
  MoveboutsIdx_cleaned <-MoveboutsIdx #[which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
  
  meanBoutSpeed <- median(vEventSpeed_smooth[MoveboutsIdx_cleaned])
  
  ##Binarize , Use indicator function 1/0 for frames where Motion Occurs
  vMotionBout <- vEventSpeed_smooth
  vMotionBout[ 1:NROW(vMotionBout) ]   <- 0
  vMotionBout[ MoveboutsIdx_cleaned  ] <- 1 ##Set Detected BoutFrames As Motion Frames
  
  vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
  #vMotionBout_rle <- rle(vMotionBout)
  
  ##Invert Speed / And Use Peak Finding To detect Bout Edges (Troughs are Peaks in the Inverse image)
  boutEdgesIdx <- find_peaks((max(vEventSpeed_smooth)- vEventSpeed_smooth)*100,15)
  vEventAccell_smooth_Onset <- boutEdgesIdx
  vEventAccell_smooth_Offset <- boutEdgesIdx

  #X11()
 #plot(vEventAccell_smooth, type='l', main="Unprocessed Cut Points")
 #points(vEventAccell_smooth_Onset,vEventAccell_smooth[vEventAccell_smooth_Onset])
 #points(vEventAccell_smooth_Offset,vEventAccell_smooth[vEventAccell_smooth_Offset],pch=6)
  
  
  ##Bout On Points Are Found At the OnSet Of the Rise/ inflexion Point - Look for Previous derivative /Accelleration change
  vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
  
  stopifnot(NROW(vMotionBout_On) > 0) ##No Bouts Detected
  
  ##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
  vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] 
  iPairs <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event
  
  ##Fix Detected Bout Points-
  ##Remove The Motion Regions Where A Peak Was not detected / Only Keep The Bouts with Peaks
  ## Shift the time of Bout to the edges where Start is on the rising foot of speed and stop are on the closest falling foot (This produces dublicates that are removed later)
  
  vMotionBout[1:length(vMotionBout)] = 0 ##Reset / Remove All Identified Movement
  for (i in 1:iPairs)
  {
    ###Motion Interval belongs to a detect bout(peak)  // Set Frame Indicators vMotionBout To Show Bout Frames
    if (any( MoveboutsIdx_cleaned >= vMotionBout_On[i] & MoveboutsIdx_cleaned < vMotionBout_Off[i] ) == TRUE)
    { 
      ##Fix Bout Onset Using Accelleration To Detect When Bout Actually Began
      ##Find Closest Speed Onset
      ##Calculate TimeDiff Between Detected BoutOnset And Actual Accelleration Onsets - Find the Onset Preceding the Detected Bout 
      OnSetTD <- vMotionBout_On[i] - vEventAccell_smooth_Onset[!is.na(vEventAccell_smooth_Onset)]
      ##Shift To Correct Onset Of Speed Increase / Denoting Where Bout Actually Began ##FIX ONSETS 
        ###Leave Out For Now
       if (NROW( (OnSetTD[OnSetTD > 0  ]) )>0) ##If Start Of Accellaration For this Bout Can Be Found / Fix It otherwise Leave it alone
        vMotionBout_On[i] <-  vMotionBout_On[i] - min(OnSetTD[OnSetTD > 0  ]) 
      
      ##FIX OFFSET to The Last MotionBoutIdx Detected Before the next BoutStart (where Decellaration Ends and A new One Begins)
      OffSetTD <- vEventAccell_smooth_Offset[!is.na(vEventAccell_smooth_Offset)] - vMotionBout_Off[i]  
      if (NROW(OffSetTD[OffSetTD > 0  ]) > 0) ##If An Offset Can Be Found (Last Bout Maybe Runs Beyond Tracking Record)
      { ##Find Last Detected Point In bout
        idxMaxEndOfBout <- which(OffSetTD == min(OffSetTD[OffSetTD > 0  ]))
        ##Last Detected Frame Before End is:
        TDNearestBout <- -(vMotionBout_Off - vEventAccell_smooth_Offset[idxMaxEndOfBout]) ##Invert Sign so as to detect TDs preceding the end 
        ##Find Which MotionBout Idx is the The Last One
        idxDetectedLastFrameOfBout <- which(TDNearestBout == min(TDNearestBout[TDNearestBout>0]) )
        vMotionBout_Off[i] <-  vMotionBout_Off[idxDetectedLastFrameOfBout] ##+ min(OffSetTD[OffSetTD > 0  ]) ##Shift |Forward To The End Of The bout
      }
      
        vMotionBout[vMotionBout_On[i]:(vMotionBout_Off[i]) ] = 1 ##Set As Motion Frames
    }
    else
    {##Remove the Ones That Do not Have a peak In them
      vMotionBout_On[i] = NA 
      vMotionBout_Off[i] = NA
    }
    
    
  } ###For Each Pair Of On-Off MotionBout 
  
  ##In Case On/Off Motion Becomes COntigious Then RLE will fail to detect it - So Make Sure Edges are there
  #vMotionBout[vMotionBout_On+1] = 1
  #vMotionBout[vMotionBout_Off] = 0 ##Make Sure Off Remains / For Rle to Work
  
  
  #X11()
  #plot(vEventAccell_smooth,type='l',main="Processed Cut-Points")
  #points(vMotionBout_On,vEventAccell_smooth[vMotionBout_On])
  #points(vMotionBout_Off,vEventAccell_smooth[vMotionBout_Off],pch=6)
  
  
  ##Get Bout Statistics #### NOt Used / Replaced##
  #vMotionBoutDuration_msec <- vMotionBout_Off[1:iPairs]-vMotionBout_On[1:iPairs]
  #vMotionBoutDuration_msec <- 1000*vMotionBoutDuration_msec[!is.na(vMotionBoutDuration_msec)]/Fs
  #vMotionBoutIntervals_msec <- 1000*(vMotionBout_On[3:(iPairs)] - vMotionBout_Off[2:(iPairs-1)])/Fs
  ############################
  
  
  ## Get Bout Statistics Again Now Using Run Length Encoding Method 
  ## Take InterBoutIntervals in msec from Last to first - 
  vMotionBout_rle <- rle(vMotionBout)
  lastBout <- max(which(vMotionBout_rle$values == 1))
  firstBout <- min(which(vMotionBout_rle$values[2:lastBout] == 1)+1) ##Skip If Recording Starts With Bout , And Catch The One After the First Pause
  vMotionBoutIBI <-1000*vMotionBout_rle$lengths[seq(lastBout-1,1,-2 )]/Fs #' IN msec and in reverse Order From Prey Capture Backwards
  ##Now That Indicators Have been integrated On Frames - Redetect On/Off Points
  vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
  vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
  vMotionBout_Off <- which(vMotionBout_OnOffDetect == -1)+1
  vMotionBoutDuration <-1000*vMotionBout_rle$lengths[seq(lastBout,2,-2 )]/Fs
  
  
  ## Denotes the Relative Time of Bout Occurance as a Sequence 1 is first, ... 10th -closer to Prey
  boutSeq <- seq(NROW(vMotionBoutIBI),1,-1 ) 
  
  ## TODO FIx these
  vMotionBoutDistanceToPrey_mm <- vDistToPrey[vMotionBout_On]*DIM_MMPERPX
  vMotionBoutDistanceTravelled_mm <- (vEventPathLength[vMotionBout_Off[1:iPairs]]-vEventPathLength[vMotionBout_On[1:iPairs]])*DIM_MMPERPX ##The Power of A Bout can be measured by distance Travelled
  
  ##Reverse Order 
  vMotionBoutDistanceToPrey_mm <- vMotionBoutDistanceToPrey_mm[boutSeq] 
  vMotionBoutDistanceTravelled_mm <- vMotionBoutDistanceTravelled_mm[boutSeq]
  
  ##Check for Errors
  stopifnot(vMotionBout_rle$values[NROW(vMotionBout_rle$lengths)] == 0 )
  stopifnot(vMotionBout_rle$values[firstBout+1] == 0 ) ##THe INitial vMotionBoutIBI Is not Actually A pause interval , but belongs to motion!
  
  ##Combine and Return
  datMotionBout <- cbind(boutSeq,vMotionBout_On,vMotionBout_Off,vMotionBoutIBI,vMotionBoutDuration,vMotionBoutDistanceToPrey_mm,vMotionBoutDistanceTravelled_mm) ##Make Data Frame
  
  
  #### PLOT DEBUG RESULTS ###
##Make Shaded Polygons
  vEventPathLength_mm<- vEventPathLength*DIM_MMPERPX
  lshadedBout <- list()
  for (i in 1:NROW(vMotionBout_Off))  
  {
  lshadedBout[[i]] <- rbind(
                     cbind(vMotionBout_Off[i],vEventSpeed_smooth[vMotionBout_Off[i]]-1),
                     cbind(vMotionBout_Off[i], max(vEventPathLength_mm) ), #vEventPathLength_mm[vMotionBout_Off[i]]+15),
                     cbind(vMotionBout_On[i], max(vEventPathLength_mm) ),#vEventPathLength_mm[vMotionBout_On[i]]+15),
                     cbind(vMotionBout_On[i], vEventSpeed_smooth[vMotionBout_On[i]]-1)
                     )
  }
  ##Plot Displacement and Speed(Scaled)
  X11()
  plot(vEventPathLength_mm,ylab="mm/sec",ylim=c(-3,max(vEventPathLength_mm[!is.na(vEventPathLength_mm)])  ),type='l',lwd=3) ##PLot Total Displacemnt over time
  lines(vEventSpeed_smooth,type='l',col="blue")
  #lines(vTailDispFilt*DIM_MMPERPX,type='l',col="magenta")
  points(MoveboutsIdx,vEventSpeed_smooth[MoveboutsIdx],col="black")
  points(MoveboutsIdx_cleaned,vEventSpeed_smooth[MoveboutsIdx_cleaned],col="red")
  points(vMotionBout_On,vEventSpeed_smooth[vMotionBout_On],col="blue",pch=17,lwd=3)
  segments(vMotionBout_Off,vEventSpeed_smooth[vMotionBout_Off]-1,vMotionBout_Off,vEventPathLength[vMotionBout_Off]+15,lwd=1.2,col="purple")
  points(vMotionBout_Off,vEventSpeed_smooth[vMotionBout_Off],col="purple",pch=14,lwd=3)
  points(boutEdgesIdx,vEventSpeed_smooth[boutEdgesIdx],col="red",pch=8,lwd=3) 
  segments(vMotionBout_On,vEventSpeed_smooth[vMotionBout_On]-1,vMotionBout_On,vEventPathLength[vMotionBout_On]+15,lwd=0.9,col="green")
  for (poly in lshadedBout)
    polygon(poly,density=3,angle=-45) 
  
  lines(vMotionBoutDistanceToPrey_mm,col="purple",lw=2)
  text(vMotionBout_On+(vMotionBout_Off-vMotionBout_On )/2,max(vEventSpeed_smooth)+3) ##Show Bout Sequence IDs to Debug Identification  
  #legend(1,100,c("PathLength","FishSpeed","TailMotion","BoutDetect","DistanceToPrey" ),fill=c("black","blue","magenta","red","purple") )
  message(paste("Number oF Bouts:",NROW(datMotionBout)))
  #dev.copy(png,filename=paste(strPlotExportPath,"/Movement-Bout_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
  
  #dev.off()
  
  
  ## Plot The Start Stop Motion Bout Binarized Data
  vMotionBout[is.na(vMotionBout)] <- 0
  vMotionBout_On[is.na(vMotionBout_On)] <- 0
  vMotionBout_Off[is.na(vMotionBout_Off)] <- 0
  X11()
  #plot(vMotionBout,type='p',xlim=c(0,max(vMotionBout_Off) )  )
  plot(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red",type='p')
  segments(seq(1:NROW(vMotionBout)),vMotionBout,seq(1:NROW(vMotionBout)),vMotionBout+0.04,lwd=0.2)
  
  points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green",pch=2,cex=2) ##On
  points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="purple",pch=13,cex=2)##Off
  
  
  return(datMotionBout)
}
