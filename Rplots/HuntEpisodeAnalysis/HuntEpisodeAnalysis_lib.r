### Library For Analysis of Hunt Episodes ###

library(signal)
library(MASS)
library(mclust,quietly = TRUE)

#library(sBIC)

citation("mclust")

## Processes The Noise IN the recorded Frames of Fish#'s Eye and Tail Motion
##Filters Fish Records - For Each Prey ID Separatelly
## As Each Row Should Contain a unique fish snapshot and not Repeats for each separate Prey Item - Ie A frame associated with a single preyID
## returns the original data frame, now with L/R Eye angles and tail motion having been filtered
filterEyeTailNoise <- function(datFishMotion)
{
  
  if (NROW(datFishMotion) < 2)
    return(datFishMotion)
  #dir.create(strFolderName )
  ##Remove NAs

    lMax <- 55
    lMin <- -20
    #spectrum(datFishMotion$LEyeAngle)
    datFishMotion$LEyeAngle <- clipEyeRange(datFishMotion$LEyeAngle,lMin,lMax)
    datFishMotion$LEyeAngle <-medianf(datFishMotion$LEyeAngle,nEyeFilterWidth)
    datFishMotion$LEyeAngle[is.na(datFishMotion$LEyeAngle)] <- 0
    datFishMotion$LEyeAngle <-filtfilt(bf_eyes,datFishMotion$LEyeAngle) # filtfilt(bf_eyes, medianf(datFishMotion$LEyeAngle,nFrWidth)) #meanf(datHuntEventMergedFrames$LEyeAngle,20)
    
    #X11()
    #lines(medianf(datFishMotion$LEyeAngle,nFrWidth),col='red')
    #lines(datFishMotion$LEyeAngle,type='l',col='blue')
    ##Replace Tracking Errors (Values set to 180) with previous last known value
    
    lMax <- 15
    lMin <- -50
    datFishMotion$REyeAngle <- clipEyeRange(datFishMotion$REyeAngle,lMin,lMax)
    datFishMotion$REyeAngle <-medianf(datFishMotion$REyeAngle,nEyeFilterWidth)
    datFishMotion$REyeAngle[is.na(datFishMotion$REyeAngle)] <- 0
    datFishMotion$REyeAngle <- filtfilt(bf_eyes,datFishMotion$REyeAngle  ) #meanf(datHuntEventMergedFrames$REyeAngle,20)
    #datFishMotion$REyeAngle <-medianf(datFishMotion$REyeAngle,nFrWidth)
  
    
    ##Fix Angle Circular Distances by DiffPolar Fix on Displacement and then Integrate back to Obtain fixed Angles  
    lMax <- +75; lMin <- -75 ;
    datFishMotion$DThetaSpine_7 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_7))+datFishMotion$DThetaSpine_7[1]  )
    datFishMotion$DThetaSpine_6 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_6))+datFishMotion$DThetaSpine_6[1] )
    datFishMotion$DThetaSpine_5 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_5))+datFishMotion$DThetaSpine_5[1]  )
    datFishMotion$DThetaSpine_4 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_4))+datFishMotion$DThetaSpine_4[1] )
    datFishMotion$DThetaSpine_3 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_3))+datFishMotion$DThetaSpine_3[1] )
    datFishMotion$DThetaSpine_2 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_2))+datFishMotion$DThetaSpine_2[1] )
    datFishMotion$DThetaSpine_1 <- filtfilt(bf_tail, cumsum(diffPolar( datFishMotion$DThetaSpine_1))+datFishMotion$DThetaSpine_1[1] )

  
  return(datFishMotion)
}


## 
#Returns F: corresponding frequencies for the constracuted wavelet scales given #Octaves and Voices
#Fc Assumes Centre Frequency For Wavelet Function (ie morlet)
getfrqscales <- function(nVoices,nOctaves,Fs,w0)
{
  a0 <- 2^(1/nVoices)
  
  #For example, assume you are using the CWT and you set your base to s0=21/12.
  #To attach physical significance to that scale, you must multiply by the sampling interval Δt, 
  #so a scale vector covering approximately four octaves with the sampling interval taken into account is sj_0 Δt j=1,2,..48. 
  #Note that the sampling interval multiplies the scales, it is not in the exponent. For discrete wavelet transforms the base scale is always 2.
  scales <- a0^seq(to=1,by=-1,from=nVoices*nOctaves)*1/Fs
  Fc <- pi/w0 ## Morlet Centre Frequency is 1/2 when w0=2*pi
  Frq <- Fc/(scales )

  return(Frq)
}

plotTailSpectrum <- function(w)
{
  w.spec <- spectrum(w,log="no",span=10,plot=FALSE,method="pgram")
  spx <- w.spec$freq*Fs
  spy <- 2*w.spec$spec #We should also multiply the spectral density by 2 so that the area under the periodogram actually equals the variance   of the time series
  #png(filename=paste(strPlotExportPath,"/TailSpectrum_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
  
  plot(spy~spx,xlab="frequency",ylab="spectral density",type='l',xlim=c(0,60) ) 
}

## Uses Wavelets to obtain the power Fq Spectrum of the tail beat in time
## w input wave 
## returns the original object w, augmented with w.cwt w.coefSq (Power) etc.
## modal Frequencies (w.FqMod) used to detect tail beat frequency
getPowerSpectrumInTime <- function(w,Fs)
{
  ##Can Test Wavelet With Artificial Signal Sample Input Signal
  #t = seq(0,1,len=Fs)
  #w = 2 * sin(2*pi*16*t)*exp(-(t-.25)^2/.001)
  #w= w + sin(2*pi*128*t)*exp(-(t-.55)^2/.001)
  #w= w + sin(2*pi*64*t)*exp(-(t-.75)^2/.001)
  #w = ts(w,deltat=1/Fs)
  N_MODESAMPLES <- 30
  w.Fs <- Fs
  w.nVoices <- 12
  w.nOctaves <- 32
  w.W0 <- 2*pi
  w.cwt <- cwt(w,noctave=w.nOctaves,nvoice=w.nVoices,plot=FALSE,twoD=TRUE,w0=w.W0)
  w.coefSq <- Mod(w.cwt)^2 #Power

  w.Frq <- getfrqscales(w.nVoices,w.nOctaves,w.Fs,w.W0)
  
  ###Make Vector Of Maximum Power-Fq Per Time Unit
  vFqMed <- rep(0,NROW(w.coefSq))
  for (i in 1:NROW(w.coefSq) )
  {
    ##Where is the Max Power at Each TimeStep?
    idxDomFq <- which(w.coefSq[i,NROW( w.Frq):1] == max(w.coefSq[i,NROW(w.Frq):1]))
    FqRank <- which(rank(w.coefSq[i,NROW(w.Frq):1] ) > (NROW(w.Frq)-N_MODESAMPLES)  )
    vFqMed[i] <- median(w.Frq[FqRank]) # sum(w.coefSq[i,NROW(w.Frq):1]*w.Frq)/sum(w.Frq) #/sum(w.coefSq[i,NROW(w.Frq):1]) #w.Frq[idxDomFq] #max(coefSq[i,idxDomFq]*Frq[idxDomFq]) #sum(coefSq[i,NROW(Frq):1]*Frq)/sum(Frq) #lapply(coefSq[,NROW(Frq):1],median)
  }
  X11();
  
  w.FqMod <-vFqMed #
 # X11()
  plot(vFqMed)
  return (list(wavedata=w,nVoices=w.nVoices,nOctaves=w.nOctaves ,MorletFrequency=w.W0, cwt=w.cwt,cwtpower=w.coefSq,Frq=w.Frq,freqMode=w.FqMod,Fs=w.Fs) )
}


## w <- Object containing Filtered Tail Segment motion (Usually the Delta angles of last 2 segments combined )
##     w.cwt <- The Continuous Wavelet Transform 
#     w.nVoices
#     w.nOctaves
## returns
plotTailPowerSpectrumInTime <- function(lwlt)
{
  
  #scales <- a0^seq(to=1,by=-1,from=nVoices*nOctaves)*1/Fs
  #Fa <- 1/2 ## Morlet Centre Frequency is 1/2 when w0=2*pi
  #Frq <- Fa/(scales )
  #Frequencies = cbind(scale=scales*(1/Fs), Frq, Period = 1./Frq)
  
  Frq <- lwlt$Frq
  #
  #plot(raster((  (vTailDisp.cwt)*1/Fs ) ), )
  #print(plot.cwt(tmp,xlab="time (units of sampling interval)"))

  collist<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
  ColorRamp<-colorRampPalette(collist)(10000)
  image(x=(1000*1:NROW(lwlt$cwtpower)/lwlt$Fs),y=Frq,z=lwlt$cwtpower[,NROW(Frq):1]
        ,useRaster=FALSE
        ,main="Frequency Content Of TailBeat"
        ,xlab="Time (msec)"
        ,ylab ="Beat Frequency (Hz)"
        ,ylim=c(0,60)
        ,col=ColorRamp
  )
  #contour(coefSq,add=T)
  #plot(coefSq[,13]   ,type='l') ##Can Plot Single Scale Like So
  
  
}

##Clusters Fish Speed Measurements into Bout And Non Bout
##Use 3 For Better Discrimination When  There Are Exist Bouts Of Different Size
detectMotionBouts <- function(vEventSpeed)
{
  nNumberOfComponents = 17
  nSelectComponents = 6
  colClass <- c("#FF0000","#04A022","#0000FF")
  
  nRec <- NROW(vEventSpeed)
  ##Fix Length Differences
  x  <-  vEventSpeed[1:nRec]
  #t <- datRenderHuntEvent$frameN
  
  #X11();plot(pvEventSpeed,pvTailDispFilt,type='p')
  
  #X11();plot(pvEventSpeed,type='p')
  #BIC <- mclustBIC(dEventSpeed)
  
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  ##prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  #modelNames = "EII"
  
  
  #fitBIC <- mclustBIC(x ,G=1:(2*nNumberOfComponents),prior =  priorControl(functionName="defaultPrior", mean=c(c(0.01),c(0.01),c(0.05),c(0.02),c(0.4),c(1.5)) ,shrinkage=0.1 ) )
  #message(attr(fitBIC,"returnCodes"))
  #plot(fitBIC)
  
  
  fit <- Mclust(x ,G=nNumberOfComponents,modelNames = "V",prior =  priorControl(functionName="defaultPrior", mean=c(c(0.01),c(0.01),c(0.05),c(0.02),c(0.4),c(1.5)),shrinkage=0.1 ) )  
  # "VVV" check out doc mclustModelNames
  #fit <- Mclust(xy ,G=2, ,prior =  priorControl(functionName="defaultPrior", mean=c(c(0.005,0),c(0.5,15)),shrinkage=0.8 ) )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  
  #fit <- Mclust(xy ,G=3 )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  summary(fit)
  
  #  X11()
  #plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
  # rug(xy)
  
  #X11()
  
  #plot(pvEventSpeed[1:nRec],type='l',col=colClass[1])
  #points(which(boutClass == 3), pvEventSpeed[boutClass == 3],type='p',col=colClass[2])
  
  ##Find Which Cluster Contains the Highest Peaks
  boutClass <- fit$classification
  clusterActivity <- vector()
  for (i in unique(boutClass))
    clusterActivity[i] <- max(x[boutClass == i])#,mean(pvEventSpeed[boutClass == 2]),mean(pvEventSpeed[boutClass == 3]))
  #clusterActivity <- c(mean(pvEventSpeed[boutClass == 1]),mean(pvEventSpeed[boutClass == 2]))
  
  clusterActivity[is.na(clusterActivity)] <- 0
  #boutCluster <- which(clusterActivity == max(clusterActivity))
  ##Select the Top nSelectComponents of clusterActivity
  boutCluster <- c(which(rank(clusterActivity) >  (nNumberOfComponents-nSelectComponents) ))   
  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which(fit$classification %in% boutCluster ) )
  #return (which( fit$z[,3]> fit$z[,1]*prior_factor1 | fit$z[,3]> fit$z[,2]*prior_factor2    )) #
  
}


##Clusters Fish Speed Measurements into Bout And Non Bout
##Use 3 For Better Discrimination When  There Are Exist Bouts Of Different Size
detectTailBouts <- function(vTailMotionFq)
{
  nNumberOfComponents = 10
  nSelectComponents = 4
  colClass <- c("#FF0000","#04A022","#0000FF")
  
  nRec <- NROW(vTailMotionFq)
  ##Fix Length Differences
  x  <-  vTailMotionFq[1:nRec]
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  ##prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  #modelNames = "EII"
  ###I can test For Possibility Of Clustering With G=n using mclustBIC returnCodes - When 0 Its succesfull
  #fitBIC <- mclustBIC(x ,G=1:(3*nNumberOfComponents),prior =  priorControl(functionName="defaultPrior", mean=c(c(0.01),c(5),c(10),c(20),c(40),c(8.5)) ,shrinkage=0.1 ) )
  #message(attr(fitBIC,"returnCodes"))
  #plot(fitBIC)
  
  fit <- Mclust(x ,G=nNumberOfComponents,modelNames = "E",prior =  priorControl(functionName="defaultPrior", mean=c(c(0.01),c(0.01),c(0.05),c(0.02),c(0.4),c(1.5)),shrinkage=0.1 ) )  
  # "VVV" check out doc mclustModelNames
  #fit <- Mclust(xy ,G=2, ,prior =  priorControl(functionName="defaultPrior", mean=c(c(0.005,0),c(0.5,15)),shrinkage=0.8 ) )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  
  #fit <- Mclust(xy ,G=3 )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  summary(fit)
  
  #  X11()
  #plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
  # rug(xy)

  ##Find Which Cluster Contains the Highest Peaks
  boutClass <- fit$classification
  clusterActivity <- vector()
  for (i in unique(boutClass))
    clusterActivity[i] <- max(x[boutClass == i])#,mean(pvEventSpeed[boutClass == 2]),mean(pvEventSpeed[boutClass == 3]))
  #clusterActivity <- c(mean(pvEventSpeed[boutClass == 1]),mean(pvEventSpeed[boutClass == 2]))
  
  clusterActivity[is.na(clusterActivity)] <- 0
  #boutCluster <- which(clusterActivity == max(clusterActivity))
  ##Select the Top nSelectComponents of clusterActivity
  boutCluster <- c(which(rank(clusterActivity) >  (nNumberOfComponents-nSelectComponents) ))   
  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which(fit$classification %in% boutCluster ) )
  #return (which( fit$z[,3]> fit$z[,1]*prior_factor1 | fit$z[,3]> fit$z[,2]*prior_factor2    )) #
  
}

detectTurnBouts <- function(vTurnSpeed,vTailDispFilt)
{
  nNumberOfComponents = 8
  nSelectComponents = 3
  
  
  nRec <- min(NROW(vTailDispFilt),NROW(vTurnSpeed))
  ##Fix Length Differences
  pvEventSpeed <-  abs(vTurnSpeed[1:nRec])
  pvTailDispFilt <-  abs(vTailDispFilt[1:nRec])
  #t <- datRenderHuntEvent$frameN
  
  #X11();plot(pvEventSpeed,pvTailDispFilt,type='p')
  
  xy <- cbind(pvEventSpeed,pvTailDispFilt)
  #X11();plot(pvEventSpeed,type='p')
  #BIC <- mclustBIC(dEventSpeed)
  
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  ##prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  fit <- Mclust(xy ,G=nNumberOfComponents,modelNames = "VII", prior =  priorControl(functionName="defaultPrior", mean=c(c(0.05,1),c(0.05,20),c(1.5,15),c(2.5,20)),shrinkage=0.1 ) )  
  summary(fit)
  
  boutClass <- fit$classification
  clusterActivity <- vector()
  for (i in unique(boutClass))
    clusterActivity[i] <- max(pvEventSpeed[boutClass == i])#,mean(pvEventSpeed[boutClass == 2]),mean(pvEventSpeed[boutClass == 3]))
  #clusterActivity <- c(mean(pvEventSpeed[boutClass == 1]),mean(pvEventSpeed[boutClass == 2]))
  
  #boutCluster <- which(clusterActivity == max(clusterActivity))
  boutCluster <- c(which(rank(clusterActivity) >  (nNumberOfComponents-nSelectComponents) ))   
  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which(fit$classification %in% boutCluster ) )
  
  
}
##Clusters Fish Speed Measurements into Bout And Non Bout
##Use 3 For Better Discrimination When  There Are Exist Bouts Of Different Size
detectMotionBouts2 <- function(vEventSpeed,vTailDispFilt)
{
  nNumberOfComponents = 5
  nSelectComponents = 3
  prior_factor2 <- 0.90 ## Adds a prior shift in the threshold Of Classification
  prior_factor1 <- 1.0 ## Adds a prior shift in the threshold Of Classification 
  colClass <- c("#FF0000","#04A022","#0000FF")
  
  nRec <- min(NROW(vTailDispFilt),NROW(vEventSpeed))
  ##Fix Length Differences
  pvEventSpeed <-  vEventSpeed[1:nRec]
  pvTailDispFilt <-  abs(vTailDispFilt[1:nRec])
  #t <- datRenderHuntEvent$frameN
  
  #X11();plot(pvEventSpeed,pvTailDispFilt,type='p')
  
  xy <- cbind(pvEventSpeed,pvTailDispFilt)
  #X11();plot(pvEventSpeed,type='p')
  #BIC <- mclustBIC(dEventSpeed)
  
  ### INcreased to 3 Clusters TO Include Other Non-Bout Activity
  ##prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  #modelNames = "EII"
  fit <- Mclust(xy ,G=nNumberOfComponents,modelNames = "VII",prior =  priorControl(functionName="defaultPrior", mean=c(c(0.01,0.1),c(0.01,5),c(0.05,5),c(0.02,2),c(0.4,20),c(1.5,25)),shrinkage=0.1 ) )  
  # "VVV" check out doc mclustModelNames
  #fit <- Mclust(xy ,G=2, ,prior =  priorControl(functionName="defaultPrior", mean=c(c(0.005,0),c(0.5,15)),shrinkage=0.8 ) )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  
  #fit <- Mclust(xy ,G=3 )  #prior=priorControl(functionName="defaultPrior",shrinkage = 0) modelNames = "V"  prior =  shrinkage = 0,modelName = "VVV"
  summary(fit)
  
#  X11()
#plot(fit, what="density", main="", xlab="Velocity (Mm/s)")
# rug(xy)
  
  #X11()
  
  #plot(pvEventSpeed[1:nRec],type='l',col=colClass[1])
  #points(which(boutClass == 3), pvEventSpeed[boutClass == 3],type='p',col=colClass[2])
  
  ##Find Which Cluster Contains the Highest Peaks
  boutClass <- fit$classification
  clusterActivity <- vector()
  for (i in unique(boutClass))
    clusterActivity[i] <- max(pvEventSpeed[boutClass == i])#,mean(pvEventSpeed[boutClass == 2]),mean(pvEventSpeed[boutClass == 3]))
  #clusterActivity <- c(mean(pvEventSpeed[boutClass == 1]),mean(pvEventSpeed[boutClass == 2]))
  
  #boutCluster <- which(clusterActivity == max(clusterActivity))
  ##Select the Top nSelectComponents of clusterActivity
  boutCluster <- c(which(rank(clusterActivity) >  (nNumberOfComponents-nSelectComponents) ))   
  #points(which( fit$z[,2]> fit$z[,1]*prior_factor ), dEventSpeed[ fit$z[,2]> fit$z[,1]*prior_factor  ],type='p',col=colClass[3])
  ## Add Prior Bias to Selects from Clusters To The 
  return (which(fit$classification %in% boutCluster ) )
  #return (which( fit$z[,3]> fit$z[,1]*prior_factor1 | fit$z[,3]> fit$z[,2]*prior_factor2    )) #
  
}

## Distance To Prey Handling  -- Fixing missing Values By Interpolation Using Fish Motion##
## Can Extend Beyond Last Frame Of Where Prey Was Last Seen , By X Frames
interpolateDistToPrey <- function(vDistToPrey,vEventSpeed_smooth, frameRegion = NA)
{
  if (!is.na(frameRegion))
    recLength <- NROW(frameRegion)
  else
    recLength <- NROW(vDistToPrey) 
  
 # stopifnot(recLength <= NROW(vEventSpeed_smooth)) ##Check For Param Error
  
  vDistToPreyInt <- rep(NA,recLength) ##Expand Dist To Prey To Cover Whole Motion Record
  vDistToPreyInt[1:recLength] <-vDistToPrey[1:recLength] ## ##Place Known Part of the vector
  
    ##Calc Speed - And Use it To Merge The Missing Values 
  vSpeedToPrey         <- c(diff(vDistToPreyInt,lag=1,differences=1),NA)
  
  vSpeedToPrey[is.na(vSpeedToPrey)] <- vEventSpeed_smooth[which(is.na(vSpeedToPrey))] ##Complete The Missing Speed Record To Prey By Using ThE fish Speed as estimate
  
  ## Interpolate Missing Values from Fish Speed - Assume Fish Is moving to Prey ##
  ##Estimate Initial DIstance From Prey Onto Which We Add the integral of Speed, By Looking At Initial PreyDist and adding any fish displacemnt to this in case The initial dist Record Is NA
  vDistToPreyInt[!is.na(vDistToPreyInt)] <- (cumsum(vSpeedToPrey[!is.na(vDistToPreyInt)]) ) ##But diff and integration Caused a shift
  
  vDistToPreyInt[9:(NROW(vDistToPreyInt))] <- vDistToPreyInt[1:(NROW(vDistToPreyInt)-min(8,NROW(vDistToPreyInt)) )] ##Fix Time Shift - If NROW > 8
  ##Compare Mean Distance Between them - Only Where Orig. PreyDist Is not NA - To obtain Integral Initial Constant (Starting Position)
  InitDistance             <- mean(vDistToPrey[!is.na(vDistToPrey)]-vDistToPreyInt[!is.na(vDistToPrey)],na.rm = TRUE )  ##vDistToPrey[!is.na(vDistToPrey)][1] + sum(vEventSpeed_smooth[(1:which(!is.na(vDistToPrey))[1])])
  
  #Add The Missing Speed Values From The fish Speed, Assumes Prey Remains Fixed And Fish Moves towards Prey
  #SpeedToPrey[is.na(vSpeedToPrey)] <- vEventSpeed_smooth[is.na(vSpeedToPrey)] ##Complete The Missing Speed Record To Prey By Using ThE fish Speed as estimate

  #vDistToPreyInt[is.na(vDistToPreyInt)] <- (cumsum(vSpeedToPrey[is.na(vDistToPreyInt)]) ) ##Add The Uknown Bit Using THE Fish's Speed and assuming the Prey Position Remains Fixed
  vDistToPreyInt <- (cumsum(vSpeedToPrey) ) ##Add The Uknown Bit Using THE Fish's Speed and assuming the Prey Position Remains Fixed
  
  vDistToPrey_Fixed <- InitDistance +  vDistToPreyInt# (cumsum(vSpeedToPrey))) ## From Initial Distance Integrate the Displacents / need -Ve Convert To Increasing Distance
  
  
  #X11() ##Compare Estimated To Recorded Prey Distance
  #plot(vDistToPrey_Fixed,type='l')
  #lines(vDistToPrey,type='l',col="blue")
  #legend()
  
  return(vDistToPrey_Fixed)
}


##############################
## Identify Bout Sections and Get Data On Durations etc.
##Uses The Detected Regions Of Bouts to extract data, on BoutOnset-Offset - Duration, Distance from Prey and Bout Power as a measure of distance moved during bout
## Note: Incomplete Bouts At the end of the trajectory will be discarted  
## regionToAnalyse - Sequence of Idx On Which To Obtain Bout Motion Data - Usually Set from 1st to last point of prey capture for a specific Prey Item
calcMotionBoutInfo2 <- function(MoveboutsIdx,vEventSpeed_smooth,vDistToPrey,vTailMotion,regionToAnalyse,plotRes=FALSE)
{
  MoveboutsIdx_cleaned <- MoveboutsIdx[MoveboutsIdx %in% regionToAnalyse]  #[which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
  
  meanBoutSpeed <- median(vEventSpeed_smooth[MoveboutsIdx_cleaned])
  
  ##Binarize , Use indicator function 1/0 for frames where Motion Occurs
  vMotionBout <- vEventSpeed_smooth
  vMotionBout[ 1:NROW(vMotionBout) ]   <- 0
  vMotionBout[ MoveboutsIdx_cleaned  ] <- 1 ##Set Detected BoutFrames As Motion Frames
  
  ##Make Initial Cut So There is always a Bout On/Off 1st & Last Frame Is always a pause
  vMotionBout[1] <- 0
  vMotionBout[2] <- 1
  vMotionBout[NROW(vMotionBout)] <- 0
  
  vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
  ##Detect Speed Minima
  boutEdgesIdx <- find_peaks((max(vEventSpeed_smooth)- vEventSpeed_smooth)*100,Fs/5)
  
  
  ##Bout On Points Are Found At the OnSet Of the Rise/ inflexion Point - Look for Previous derivative /Accelleration change
  vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1

  if (NROW(vMotionBout_On) == 1)
    warning("No Bout Onset Detected")
    

  vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] 
  iPairs <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event

      
  ##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
  ## Get Bout Statistics Again Now Using Run Length Encoding Method 
  ## Take InterBoutIntervals in msec from Last to first - 
  vMotionBout_rle <- rle(vMotionBout)
  ##Filter Out Small Bouts/Pauses -
  idxShort <- which(vMotionBout_rle$lengths < MIN_BOUT_DURATION)
  for (jj in idxShort)
  {
    ##Fill In this Gap
    idxMotionStart <- sum(vMotionBout_rle$length[1:(jj-1)])
    idxMotionEnd <- idxMotionStart + vMotionBout_rle$length[jj]
    
    if( vMotionBout_rle$values[jj] == 1) # If this is a Motion
      vMotionBout[idxMotionStart:idxMotionEnd] <- 0 ##Replace short motion with Pause

    if( vMotionBout_rle$values[jj] == 0) # If this is a Motion
      vMotionBout[idxMotionStart:idxMotionEnd] <- 1 ##Replace short Pause with Motion
    
  }
  ##Make Initial Cut So 1st & Last Frame Is always a pause
  vMotionBout[1] <- 0
  vMotionBout[NROW(vMotionBout)] <- 0
  
  ##Redo Fixed Binary Vector
  vMotionBout_rle <- rle(vMotionBout)

    
  lastBout <- max(which(vMotionBout_rle$values == 1))
  firstBout <- min(which(vMotionBout_rle$values[1:lastBout] == 1)) ##Skip If Recording Starts With Bout , And Catch The One After the First Pause
  if (lastBout > firstBout) ##If More than One Bout Exists
    vMotionBoutIBI <-1000*vMotionBout_rle$lengths[seq(lastBout-1,firstBout,-2 )]/Fs #' IN msec and in reverse Order From Prey Capture Backwards
  else
    vMotionBoutIBI <- 1
  ##Add One Since IBI count is 1 less than the bout count
  vMotionBoutIBI <- c(vMotionBoutIBI,NA)
  
  ##Now That Indicators Have been integrated On Frames - Redetect On/Off Points
  vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
  vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
  vMotionBout_Off <- which(vMotionBout_OnOffDetect == -1)+1
  vMotionBoutDuration <-1000*vMotionBout_rle$lengths[seq(lastBout,firstBout,-2 )]/Fs
  
  vEventPathLength_mm<- vEventPathLength*DIM_MMPERPX
  ## Denotes the Relative Time of Bout Occurance as a Sequence 1 is first, ... 10th -closer to Prey
  boutSeq <- seq(NROW(vMotionBoutDuration),1,-1 ) ##The time Sequence Of Event Occurance (Fwd Time)
  boutRank <- seq(1,NROW(vMotionBoutDuration),1 ) ##Denotes Reverse Order - From Prey Captcha being First going backwards to the n bout
  ## TODO FIx these
  vMotionBoutDistanceToPrey_mm <- vDistToPrey[vMotionBout_On]*DIM_MMPERPX
  vMotionBoutDistanceTravelled_mm <- (vEventPathLength_mm[vMotionBout_Off[1:iPairs]]-vEventPathLength_mm[vMotionBout_On[1:iPairs]]) ##The Power of A Bout can be measured by distance Travelled
  
  ##Reverse Order 
  vMotionBoutDistanceToPrey_mm <- vMotionBoutDistanceToPrey_mm[boutSeq] 
  vMotionBoutDistanceTravelled_mm <- vMotionBoutDistanceTravelled_mm[boutSeq]
  
  ##Check for Errors
  #stopifnot(vMotionBout_rle$values[NROW(vMotionBout_rle$lengths)] == 0 )###Check End With  Pause Not A bout
  stopifnot(vMotionBout_rle$values[firstBout+1] == 0 ) ##THe INitial vMotionBoutIBI Is not Actually A pause interval , but belongs to motion!
  
  ##Combine and Return
  datMotionBout <- cbind(boutSeq,boutRank,vMotionBout_On,vMotionBout_Off,vMotionBoutIBI,vMotionBoutDuration,vMotionBoutDistanceToPrey_mm,vMotionBoutDistanceTravelled_mm) ##Make Data Frame
  
  
  #### PLOT DEBUG RESULTS ###
  ##Make Shaded Polygons
  if (plotRes)
  {
    #vEventSpeed_smooth <- vEventSpeed_smooth*5
    
    lshadedBout <- list()
    t <- seq(1:NROW(vEventPathLength_mm))/(Fs/1000)
    for (i in 1:NROW(vMotionBout_Off))  
    {
      lshadedBout[[i]] <- rbind(
        cbind(t[vMotionBout_Off[i] ],vEventSpeed_smooth[vMotionBout_Off[i]]-1),
        cbind(t[vMotionBout_Off[i] ], max(vEventPathLength_mm) ), #vEventPathLength_mm[vMotionBout_Off[i]]+15),
        cbind(t[vMotionBout_On[i] ], max(vEventPathLength_mm) ),#vEventPathLength_mm[vMotionBout_On[i]]+15),
        cbind(t[vMotionBout_On[i] ], vEventSpeed_smooth[vMotionBout_On[i]]-1)
      )
    }
    
    ##Plot Displacement and Speed(Scaled)
    vTailDispFilt <- filtfilt( bf_tailClass2, abs(filtfilt(bf_tailClass, (vTailMotion) ) ) )
    ymax <- 15 #max(vEventPathLength_mm[!is.na(vEventPathLength_mm)])
    plot(t,vEventPathLength_mm,ylab="mm",
         xlab="msec",
         ylim=c(-0.3, ymax  ),type='l',lwd=3) ##PLot Total Displacemnt over time
    par(new=TRUE) ##Add To Path Length Plot But On Separate Axis So it Scales Nicely
    par(mar=c(4,4,2,2))
    plot(t,vEventSpeed_smooth,type='l',axes=F,xlab=NA,ylab=NA,col="blue",ylim=c(0,1.5))
    axis(side = 4,col="blue")
    mtext(side = 4, line = 3, 'Speed (mm/sec)')
    
    #lines(vTailDispFilt*DIM_MMPERPX,type='l',col="magenta")
    points(t[MoveboutsIdx],vEventSpeed_smooth[MoveboutsIdx],col="black")
    points(t[MoveboutsIdx_cleaned],vEventSpeed_smooth[MoveboutsIdx_cleaned],col="red")
    points(t[vMotionBout_On],vEventSpeed_smooth[vMotionBout_On],col="blue",pch=17,lwd=3)
    segments(t[vMotionBout_Off],vEventSpeed_smooth[vMotionBout_Off]-1,t[vMotionBout_Off],vEventPathLength[vMotionBout_Off]+15,lwd=1.2,col="purple")
    points(t[vMotionBout_Off],vEventSpeed_smooth[vMotionBout_Off],col="purple",pch=14,lwd=3)
    points(t[boutEdgesIdx],vEventSpeed_smooth[boutEdgesIdx],col="red",pch=8,lwd=3) 
    segments(t[vMotionBout_On],vEventSpeed_smooth[vMotionBout_On]-1,t[vMotionBout_On],vEventPathLength[vMotionBout_On]+15,lwd=0.9,col="green")
    for (poly in lshadedBout)
      polygon(poly,density=3,angle=-45) 
    
    #lines(vMotionBoutDistanceToPrey_mm,col="purple",lw=2)
    pkPt <- round(vMotionBout_On+(vMotionBout_Off-vMotionBout_On )/2)
    text(t[pkPt],vEventSpeed_smooth[pkPt]+0.1,labels=boutSeq) ##Show Bout Sequence IDs to Debug Identification  
    #legend(1,100,c("PathLength","FishSpeed","TailMotion","BoutDetect","DistanceToPrey" ),fill=c("black","blue","magenta","red","purple") )
    
    plot(t[1:NROW(vTailMotion)],vTailMotion,type='l',
         xlab="msec",
         col="red",main="Tail Motion")
    lines(t[1:NROW(vTailMotion)],vTailDispFilt,col="black" )
    
  } ##If Plot Flag Is Set 
  
  message(paste("Number oF Bouts:",NROW(datMotionBout)))
  # dev.copy(png,filename=paste(strPlotExportPath,"/Movement-Bout_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
  

  return(datMotionBout)
}




######################OLD CODE ####################
# 
# ## Identify Bout Sections and Get Data On Durations etc.
# ##Uses The Detected Regions Of Bouts to extract data, on BoutOnset-Offset - Duration, Distance from Prey and Bout Power as a measure of distance moved during bout
# ## Note: Incomplete Bouts At the end of the trajectory will be discarted  
# ## regionToAnalyse - Sequence of Idx On Which To Obtain Bout Motion Data - Usually Set from 1st to last point of prey capture for a specific Prey Item
# calcMotionBoutInfo <- function(MoveboutsIdx,vEventSpeed_smooth,vDistToPrey,vTailMotion,regionToAnalyse,plotRes=FALSE)
# {
#   MoveboutsIdx_cleaned <- MoveboutsIdx[MoveboutsIdx %in% regionToAnalyse]  #[which(vEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
#   
#   meanBoutSpeed <- median(vEventSpeed_smooth[MoveboutsIdx_cleaned])
#   
#   ##Binarize , Use indicator function 1/0 for frames where Motion Occurs
#   vMotionBout <- vEventSpeed_smooth
#   vMotionBout[ 1:NROW(vMotionBout) ]   <- 0
#   vMotionBout[ MoveboutsIdx_cleaned  ] <- 1 ##Set Detected BoutFrames As Motion Frames
#   
#   
#   #vMotionBout_rle <- rle(vMotionBout)
#   
#   ##Invert Speed / And Use Peak Finding To detect Bout Edges (Troughs are Peaks in the Inverse image)
#   boutEdgesIdx <- find_peaks((max(vEventSpeed_smooth)- vEventSpeed_smooth)*100,Fs/5)
#   vEventAccell_smooth_Onset  <- boutEdgesIdx
#   vEventAccell_smooth_Offset <- c(boutEdgesIdx,NROW(vEventSpeed_smooth))
#   #vMotionBout[boutEdgesIdx]  <- 0 ##Set Edges As Cut Points XX
#   
#   vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
#   #X11()
#   #plot(vEventSpeed_smooth, type='l', main="Unprocessed Cut Points")
#   #points(vEventAccell_smooth_Onset,vEventSpeed_smooth[vEventAccell_smooth_Onset])
#   #points(vEventAccell_smooth_Offset,vEventSpeed_smooth[vEventAccell_smooth_Offset],pch=6)
#   
#   
#   ##Bout On Points Are Found At the OnSet Of the Rise/ inflexion Point - Look for Previous derivative /Accelleration change
#   vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
#   
#   # stopifnot(NROW(vMotionBout_On) > 0) ##No Bouts Detected
#   
#   ##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
#   vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] 
#   iPairs <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event
#   
#   ##Fix Detected Bout Points-
#   ##Remove The Motion Regions Where A Peak Was not detected / Only Keep The Bouts with Peaks
#   ## Shift the time of Bout to the edges where Start is on the rising foot of speed and stop are on the closest falling foot (This produces dublicates that are removed later)
#   
#   vMotionBout[1:length(vMotionBout)] = 0 ##Reset / Remove All Identified Movement
#   for (i in 1:iPairs)
#   {
#     ###Motion Interval belongs to a detect bout(peak)  // Set Frame Indicators vMotionBout To Show Bout Frames
#     if (any( MoveboutsIdx_cleaned >= vMotionBout_On[i] & MoveboutsIdx_cleaned < vMotionBout_Off[i] ) == TRUE)
#     { 
#       ##Fix Bout Onset Using Accelleration To Detect When Bout Actually Began
#       ##Find Closest Speed Onset
#       ##Calculate TimeDiff Between Detected BoutOnset And Actual Accelleration Onsets - Find the Onset Preceding the Detected Bout 
#       OnSetTD <- vMotionBout_On[i] - vEventAccell_smooth_Onset[!is.na(vEventAccell_smooth_Onset)]
#       ##Shift To Correct Onset Of Speed Increase / Denoting Where Bout Actually Began ##FIX ONSETS 
#       ###Leave Out For Now
#       if (NROW( (OnSetTD[OnSetTD > 0  ]) )>0) ##If Start Of Accellaration For this Bout Can Be Found / Fix It otherwise Leave it alone
#       {
#         idxMinStartOfBout <- which(OnSetTD == min(OnSetTD[OnSetTD > 0  ]))
#         TDNearestBout <- (vMotionBout_On - vEventAccell_smooth_Onset[idxMinStartOfBout]) ##Invert Sign so as to detect TDs preceding the end 
#         idxDetectedFirstFrameOfBout <- max(which(TDNearestBout == min(TDNearestBout[TDNearestBout>0]) )  ) ##max to pick the last one in case duplicate vMotionBout_Off values
#         #vMotionBout_On[i] <-  vMotionBout_On[i] - min(OnSetTD[OnSetTD > 0  ])
#         vMotionBout_On[i] <-  vMotionBout_On[idxDetectedFirstFrameOfBout]
#       }
#       
#       ##FIX OFFSET to The Last MotionBoutIdx Detected Before the next BoutStart (where Decellaration Ends and A new One Begins)
#       OffSetTD <- vEventAccell_smooth_Offset[!is.na(vEventAccell_smooth_Offset)] - vMotionBout_Off[i]  
#       if (NROW(OffSetTD[OffSetTD > 0  ]) > 0) ##If An Offset Can Be Found (Last Bout Maybe Runs Beyond Tracking Record)
#       { ##Find Last Detected Point In bout
#         idxMaxEndOfBout <- which(OffSetTD == min(OffSetTD[OffSetTD > 0  ]))
#         ##Last Detected Frame Before End is:
#         TDNearestBout <- -(vMotionBout_Off - vEventAccell_smooth_Offset[idxMaxEndOfBout]) ##Invert Sign so as to detect TDs preceding the end 
#         ##Find Which MotionBout Idx is the The Last One
#         idxDetectedLastFrameOfBout <- max(which(TDNearestBout == min(TDNearestBout[TDNearestBout>0]) )  ) ##max to pick the last one in case duplicate vMotionBout_Off values
#         # vSpeedHillEdge <- vMotionBout_Off[i] + min(OffSetTD[OffSetTD > 0  ]) ##How Far  is the next Edge
#         ##Choose Closest Either The Last Detected Point Or The End/ Edge Of Speed Hill
#         vMotionBout_Off[i] <-  vMotionBout_Off[idxDetectedLastFrameOfBout] #min(vSpeedHillEdge,vMotionBout_Off[idxDetectedLastFrameOfBout]) ##+ min(OffSetTD[OffSetTD > 0  ]) ##Shift |Forward To The End Of The bout
#       }
#       
#       vMotionBout[vMotionBout_On[i]:(vMotionBout_Off[i]) ] = 1 ##Set As Motion Frames
#     }
#     else
#     {##Remove the Ones That Do not Have a peak In them
#       vMotionBout_On[i] = NA 
#       vMotionBout_Off[i] = NA
#     }
#     
#     
#   } ###For Each Pair Of On-Off MotionBout 
#   
#   ##In Case On/Off Motion Becomes COntigious Then RLE will fail to detect it - So Make Sure Edges are there
#   #vMotionBout[vMotionBout_On+1] = 1
#   vMotionBout[vMotionBout_Off] = 0 ##Make Sure Off Remains / For Rle to Work
#   
#   #X11()
#   #plot(vEventAccell_smooth,type='l',main="Processed Cut-Points")
#   #points(vMotionBout_On,vEventAccell_smooth[vMotionBout_On])
#   #points(vMotionBout_Off,vEventAccell_smooth[vMotionBout_Off],pch=6)
#   
#   
#   ##Get Bout Statistics #### NOt Used / Replaced##
#   #vMotionBoutDuration_msec <- vMotionBout_Off[1:iPairs]-vMotionBout_On[1:iPairs]
#   #vMotionBoutDuration_msec <- 1000*vMotionBoutDuration_msec[!is.na(vMotionBoutDuration_msec)]/Fs
#   #vMotionBoutIntervals_msec <- 1000*(vMotionBout_On[3:(iPairs)] - vMotionBout_Off[2:(iPairs-1)])/Fs
#   ############################
#   
#   
#   ## Get Bout Statistics Again Now Using Run Length Encoding Method 
#   ## Take InterBoutIntervals in msec from Last to first - 
#   vMotionBout_rle <- rle(vMotionBout)
#   lastBout <- max(which(vMotionBout_rle$values == 1))
#   firstBout <- min(which(vMotionBout_rle$values[2:lastBout] == 1)+1) ##Skip If Recording Starts With Bout , And Catch The One After the First Pause
#   vMotionBoutIBI <-1000*vMotionBout_rle$lengths[seq(lastBout-1,1,-2 )]/Fs #' IN msec and in reverse Order From Prey Capture Backwards
#   ##Now That Indicators Have been integrated On Frames - Redetect On/Off Points
#   vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
#   vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
#   vMotionBout_Off <- which(vMotionBout_OnOffDetect == -1)+1
#   vMotionBoutDuration <-1000*vMotionBout_rle$lengths[seq(lastBout,2,-2 )]/Fs
#   
#   vEventPathLength_mm<- vEventPathLength*DIM_MMPERPX
#   ## Denotes the Relative Time of Bout Occurance as a Sequence 1 is first, ... 10th -closer to Prey
#   boutSeq <- seq(NROW(vMotionBoutIBI),1,-1 ) 
#   boutRank <- seq(1,NROW(vMotionBoutIBI),1 ) ##Denotes Reverse Order - From Prey Captcha being First going backwards to the n bout
#   ## TODO FIx these
#   vMotionBoutDistanceToPrey_mm <- vDistToPrey[vMotionBout_On]*DIM_MMPERPX
#   vMotionBoutDistanceTravelled_mm <- (vEventPathLength_mm[vMotionBout_Off[1:iPairs]]-vEventPathLength_mm[vMotionBout_On[1:iPairs]]) ##The Power of A Bout can be measured by distance Travelled
#   
#   ##Reverse Order 
#   vMotionBoutDistanceToPrey_mm <- vMotionBoutDistanceToPrey_mm[boutSeq] 
#   vMotionBoutDistanceTravelled_mm <- vMotionBoutDistanceTravelled_mm[boutSeq]
#   
#   ##Check for Errors
#   stopifnot(vMotionBout_rle$values[NROW(vMotionBout_rle$lengths)] == 0 )
#   stopifnot(vMotionBout_rle$values[firstBout+1] == 0 ) ##THe INitial vMotionBoutIBI Is not Actually A pause interval , but belongs to motion!
#   
#   ##Combine and Return
#   datMotionBout <- cbind(boutSeq,boutRank,vMotionBout_On,vMotionBout_Off,vMotionBoutIBI,vMotionBoutDuration,vMotionBoutDistanceToPrey_mm,vMotionBoutDistanceTravelled_mm) ##Make Data Frame
#   
#   
#   #### PLOT DEBUG RESULTS ###
#   ##Make Shaded Polygons
#   if (plotRes)
#   {
#     #vEventSpeed_smooth <- vEventSpeed_smooth*5
#     
#     lshadedBout <- list()
#     t <- seq(1:NROW(vEventPathLength_mm))/(Fs/1000)
#     for (i in 1:NROW(vMotionBout_Off))  
#     {
#       lshadedBout[[i]] <- rbind(
#         cbind(t[vMotionBout_Off[i] ],vEventSpeed_smooth[vMotionBout_Off[i]]-1),
#         cbind(t[vMotionBout_Off[i] ], max(vEventPathLength_mm) ), #vEventPathLength_mm[vMotionBout_Off[i]]+15),
#         cbind(t[vMotionBout_On[i] ], max(vEventPathLength_mm) ),#vEventPathLength_mm[vMotionBout_On[i]]+15),
#         cbind(t[vMotionBout_On[i] ], vEventSpeed_smooth[vMotionBout_On[i]]-1)
#       )
#     }
#     
#     ##Plot Displacement and Speed(Scaled)
#     vTailDispFilt <- filtfilt( bf_tailClass2, abs(filtfilt(bf_tailClass, (vTailMotion) ) ) )
#     
#     plot(t,vEventPathLength_mm,ylab="mm",
#          xlab="msec",
#          ylim=c(-0.3,max(vEventPathLength_mm[!is.na(vEventPathLength_mm)])  ),type='l',lwd=3) ##PLot Total Displacemnt over time
#     par(new=T) ##Add To Path Length Plot But On Separate Axis So it Scales Nicely
#     par(mar=c(4,4,2,2))
#     plot(t,vEventSpeed_smooth,type='l',axes=F,xlab=NA,ylab=NA,col="blue")
#     axis(side = 4,col="blue")
#     mtext(side = 4, line = 3, 'Speed (mm/sec)')
#     
#     #lines(vTailDispFilt*DIM_MMPERPX,type='l',col="magenta")
#     points(t[MoveboutsIdx],vEventSpeed_smooth[MoveboutsIdx],col="black")
#     points(t[MoveboutsIdx_cleaned],vEventSpeed_smooth[MoveboutsIdx_cleaned],col="red")
#     points(t[vMotionBout_On],vEventSpeed_smooth[vMotionBout_On],col="blue",pch=17,lwd=3)
#     segments(t[vMotionBout_Off],vEventSpeed_smooth[vMotionBout_Off]-1,t[vMotionBout_Off],vEventPathLength[vMotionBout_Off]+15,lwd=1.2,col="purple")
#     points(t[vMotionBout_Off],vEventSpeed_smooth[vMotionBout_Off],col="purple",pch=14,lwd=3)
#     points(t[boutEdgesIdx],vEventSpeed_smooth[boutEdgesIdx],col="red",pch=8,lwd=3) 
#     segments(t[vMotionBout_On],vEventSpeed_smooth[vMotionBout_On]-1,t[vMotionBout_On],vEventPathLength[vMotionBout_On]+15,lwd=0.9,col="green")
#     for (poly in lshadedBout)
#       polygon(poly,density=3,angle=-45) 
#     
#     #lines(vMotionBoutDistanceToPrey_mm,col="purple",lw=2)
#     text(t[round(vMotionBout_On+(vMotionBout_Off-vMotionBout_On )/2)],max(vEventSpeed_smooth)+3,labels=boutSeq) ##Show Bout Sequence IDs to Debug Identification  
#     #legend(1,100,c("PathLength","FishSpeed","TailMotion","BoutDetect","DistanceToPrey" ),fill=c("black","blue","magenta","red","purple") )
#     
#     plot(t[1:NROW(vTailMotion)],vTailMotion,type='l',
#          xlab="msec",
#          col="red",main="Tail Motion")
#     lines(t[1:NROW(vTailMotion)],vTailDispFilt,col="black" )
#     
#   } ##If Plot Flag Is Set 
#   
#   message(paste("Number oF Bouts:",NROW(datMotionBout)))
#   # dev.copy(png,filename=paste(strPlotExportPath,"/Movement-Bout_exp",expID,"_event",eventID,"_track",trackID,".png",sep="") );
#   
#   #  dev.off()
#   
#   
#   ## Plot The Start Stop Motion Bout Binarized Data
#   #vMotionBout[is.na(vMotionBout)] <- 0
#   #vMotionBout_On[is.na(vMotionBout_On)] <- 0
#   #vMotionBout_Off[is.na(vMotionBout_Off)] <- 0
#   
#   #X11()
#   #plot(vMotionBout,type='p',xlim=c(0,max(vMotionBout_Off) )  )
#   #plot(MoveboutsIdx_cleaned,vMotionBout[MoveboutsIdx_cleaned],col="red",type='p')
#   #segments(seq(1:NROW(vMotionBout)),vMotionBout,seq(1:NROW(vMotionBout)),vMotionBout+0.04,lwd=0.2)
#   
#   #points(vMotionBout_On,vMotionBout[vMotionBout_On],col="green",pch=2,cex=2) ##On
#   #points(vMotionBout_Off,vMotionBout[vMotionBout_Off],col="purple",pch=13,cex=2)##Off
#   
#   
#   return(datMotionBout)
# }
######################################## END OF V1 #############

