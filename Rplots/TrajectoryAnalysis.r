# ################
# ## \todo Code to Characterise Fish Trajectories Events - Locate Videos And Start Frames - ###
# ## Rank Files With Most Hunting

# M. Meyer sent doc with ideas on 30-11-17, suggests:
# * Some additional analyses
# 1. Number of tracks containing vergence events/total number of tracks= probability that a group will switch into hunting mode
# 2. For all tracks with eye vergence events calculate the average duration of vergence. 
#KL:
# Since nLarva can vary, I suggest looking at statistics of HuntRates and Durations as per fish in Group

#################

source("TrackerDataFilesImport_lib.r")
source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r") ##For Filter Initialization

meanf <- function(t,k) {n=length(t);tproc=t;k=min(k,n); for(i in (k/2):n) tproc[i]=mean(t[max(1,i-k/2): min(n, i+k/2) ]);return(tproc) }

### We are looking to detect Exploration/Exploitation (as in Marquez et al. 2020)  using a measure
### of spatial dispersion - calculated for each tracked frame, calculated as the spatial dispersion of trajectory of the preceding X secods
## Trajectory Dispersion - as min radius that can encompass the whole trajectory of last twindowSec sec.  ##"
calcTrajectoryDispersionAndLength <- function(datEventFrames,twindowSec=5)
{
  start.time <- Sys.time()
  ##Sort by frameRow - so as to Calc recent Distances in correct time sequence  
  datEventFrames <- datEventFrames[order( as.numeric(datEventFrames$frameN) ),]
  ## Note, could use combn gen All Combh of positions in trajectory - Here outer is app. faster
  #message(paste("## Trajectory Dispersion - as min radius that can encompass the whole trajectory of last twindowSec sec.  ##" ) )
  stopifnot(NROW(datEventFrames) > 0)
  ##datEventFrames <- datAllFrames[datAllFrames$expID == 218 & datAllFrames$eventID == 2 & datAllFrames$posX != 0,]
  vDispersionPerFrame <- vector()
  vDistanceTravelledToFrame <- vector()
  vMSD <- vector() #Mean Square Displacement of the twindowSec trajectory
  vSD <- vector() # Square Displacement to the start of the twindowSec trajectory
  ## Estimate number of frames for Tsec of video
  nfrm <- head(datEventFrames$fps,1)*twindowSec
  nSpace <- nfrm/4
  ## Loop Through Each Frame and calc Dispersion
  datEventFrames$posX <- meanf(datEventFrames$posX,10)
  datEventFrames$posY <- meanf(datEventFrames$posY,10)
  
  #plot(datEventFrames$posX,datEventFrames$posY)
  if (nfrm > NROW(datEventFrames))
  {
    warning("Event:",head(datEventFrames$eventID,1)," does not have enough frames to estimate dispersion \n");
    return(list(Dispersion=NA,Length=NA,MSD=NA,SD=NA))
  }
  ## \TODO: Do not have to sample all points along trajectory - Sub Sample Spaced out points to Get Dispersion Estimate
  for (i in nfrm:NROW(datEventFrames) ) ##
  {
    ##Skip some intermediates frames to optimize speed - sample equally distant points in  trajectory / instead of every point
    vIdx <- seq(from=(i-nfrm),by=nSpace,to=i)
    vX <- head(datEventFrames[vIdx,"posX"],nfrm)
    vY <- head(datEventFrames[vIdx,"posY"],nfrm)
    ##Find min Radius of circle that could encompass whole trajectory as half of the maximum distance between any two points of the trajectory
    
    ## Use outer to compute All to All point X  traj. differences / 
    mat_posDX <- outer(vX,vX,'-')
    mat_posDY <- outer(vY,vY,'-') #and Y 
    # Combine to Find Distances DX DY between all point of trajectory
    mat_ptDist <- sqrt(mat_posDX^2 + mat_posDY^2)
    # Radius of Circle Encompassing the two most distant parts of trajectory section
    vDispersionPerFrame[i] <- max(mat_ptDist)/2 
    #Calc Path Distance by Summing Successive Point Difference / on Upper Off-Diagonal of distance matrix
    vDistanceTravelledToFrame[i] <- sum(mat_ptDist[row(mat_ptDist) == (col(mat_ptDist) - 1)])
    #Calc Mean Squared Displacement from initial Point x0,y0 being the mean squared sum of the 1st distance Matrix Row
    vMSD[i] <- mean((DIM_MMPERPX*mat_ptDist[1,1:ncol(mat_ptDist)]) ^2)
    vSD[i]  <- (DIM_MMPERPX*mat_ptDist[1,ncol(mat_ptDist)])^2 #Squared Distance To point twindowSec  ago
  }
  
#  vRes100 <- vDispersionPerFrame*DIM_MMPERPX
#  plot(vDispersionPerFrame*DIM_MMPERPX)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  lRet <- list(Dispersion=vDispersionPerFrame*DIM_MMPERPX,Length=vDistanceTravelledToFrame*DIM_MMPERPX,MSD=vMSD,SD=vSD)
  return( lRet)
}
  
##Provide a data structure with organized processed data extracted from each Recording Event
calcRecordingEventSpeed <- function(datAllFrames,vexpID,vdatasetID)
{
  message(paste("## Calculate Speed Statistics for each event in dataFrame  ##" ) )
  idx <- 1
  lEventSpeed <- list()
  ## For Each Exp
  for (e in vexpID)
  {
    stopifnot(is.numeric(e) & e > 0)
    vEventID = unique((datAllFrames[datAllFrames$expID == e,]$eventID))
    ##For Each Event
   for (v in  vEventID)
   {
    datRecordingEvent <- datAllFrames[datAllFrames$expID == e & datAllFrames$eventID == v,]
    
    meanfps <-  unique(datRecordingEvent$fps)
    groupID <- as.character(unique(datRecordingEvent$groupID) )
    
    message(paste("ExpID:",e,"EventID:",v,"fps:",meanfps ) )
    
    ##Need to Identify TrackLet Units, Avoid speed calc errors due to fish going in and out of view
    #### PROCESS TrackLets ###
    vTracklets <- unique(datRecordingEvent$trackletID)
    for (t in vTracklets)
    {
      if (t == 0) ##Periodic FoodDensity Recording /. Not A tracklet
        next()
      datRecordingEvent    <- datRecordingEvent[datRecordingEvent$trackletID == t,]
      ###No Record
      if (NROW(datRecordingEvent )< 3 )
      {
        warning("Tracklet too short")
        lEventSpeed[[idx]] <- list(groupID=groupID,
                                   expID=e,eventID=v,
                                   trackletID=t,
                                   Duration_sec=0,
                                   Length_mm=0) #density(vEventSpeed_smooth,from=0,to=1,kernel="gaussian")
        idx <- idx + 1
        
        next()
      }
      vDeltaXFrames        <- diff(datRecordingEvent$posX,lag=1,differences=1)
      vDeltaYFrames        <- diff(datRecordingEvent$posY,lag=1,differences=1)
      vDeltaDisplacement   <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
    
    
      #nNumberOfBouts       <- 
      dframe               <- diff(datRecordingEvent$frameN,lag=1,differences=1)
      dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - and Convert To Time sec  
      vEventSpeed          <- meanf(vDeltaDisplacement/dframe,5) ##IN (mm) Divide Displacement By TimeFrame to get Instantentous Speed, Apply Mean Filter Smooth Out 
      
      vEventSpeed[is.na(vEventSpeed)] = 0
      vEventSpeed_smooth <- filtfilt(bf_speed, vEventSpeed) #meanf(vEventSpeed,100) #
      vEventSpeed_smooth[vEventSpeed_smooth < 0] <- 0 ## Remove -Ve Values As an artefact of Filtering
      vEventSpeed_smooth[is.na(vEventSpeed_smooth)] = 0
      vEventPathDisplacement_mm <- cumsum(vEventSpeed_smooth)*DIM_MMPERPX
      ##Plot Displacement Vs Time in Sec
      #plot(cumsum(dframe)/G_APPROXFPS,vEventPathDisplacement_mm,col="red")
      ##TODO - obtain Actual FPS of each DAtaset
      lEventSpeed[[idx]] <- list(groupID=groupID,
                                 expID=e,
                                 eventID=v,
                                 trackletID=t,
                                 Duration_sec=sum(dframe)/meanfps,
                                 Length_mm=max(vEventPathDisplacement_mm)) #density(vEventSpeed_smooth,from=0,to=1,kernel="gaussian")
      idx <- idx + 1
      
    }##For Each Tracklet
   }
    ##END For Each Event
  ##END For Each Exp
  }
  datEventSpeed <- data.frame(do.call(rbind,lEventSpeed))
  
  return(datEventSpeed )
  
}


### Returns summary measurements for overall trajectories for each larva 
## Sinuoisity, Overall distance travelled, Speed (Avg,peak) etc
#########################################################
calcMotionStat <- function(datAllFrames,vexpID,vdatasetID)
{
  message("## Start Trajectory Analysis For Group ##")
  lGroupMotion         <- list()
  lTEventDisplacement  <- list() ##Event Trajectory  Path Length
  lTEventSpeed         <- list() ##Event Trajectory  Speed
  lTEventSinuosity     <- list()

  idx             <- 0
  nLarva  	      <- length(vexpID)
  nEventsAnalysed <- 0
  nTotalFrames    <- 0
  
  for (i in vexpID)
  {
    idx                 <- idx+1; ## Only increment if Record Is to be added
    
    datLarvaFrames      <- datAllFrames[datAllFrames$expID == i,]
    vEventID            <- unique((datLarvaFrames$eventID))
    ##Reset Variables ##
    nEventsAnalysed        <- 0
    nLarvalAnalysedFrames  <- 0 ## Number of Motion Frames  Analysed For this Larva
    movementRatio          <- 0
    
    
    DataSetID             <- ifelse(any(names(datLarvaFrames) == "dataSet"),unique(datLarvaFrames$dataSet),0 )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID )
    larvaID               <- unique((datLarvaFrames$larvaID))
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFrames[datLarvaFrames$REyeAngle > -35 & datLarvaFrames$REyeAngle <  35 & datLarvaFrames$LEyeAngle > -35 & datLarvaFrames$LEyeAngle <  35,]
    
    
    lTEventDisplacement  <-list() ##Event Trajectory  Path Length
    lTEventSpeed         <-list() ##Event Trajectory  Speed
    lTEventSinuosity     <- list()

    vTDisplacement    <- list()
    vTSpeed           <- list()
    vTSinuosity       <- list()
    
    
    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
  for (k in vEventID)
  {
      ##Error Captcha ##
      if (is.na(k) | k < 1 )
      {
        stop(paste("Invalid Event Idx:",k) )
      }

	    ##Select Motion Frames of this Event ## 
      ## Filter Out Tracking Loses Set to poxX,Y=0
	    datEventFrames       <- datLarvaFrames[datLarvaFrames$eventID == k & datLarvaFrames$posX > 0 & datLarvaFrames$posY > 0,]

	    ##Consider Only if there are at least 100 frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
	    if  (NROW(datEventFrames$frameN) < 100)
	    {
	      ## Episode Trajectory is too Short Insufficient ##	
	      warning(paste("Trajectory Too short - Ignoring Trajectory  of Experiment:",i,"event",k," larva",larvaID) )
	      
	      lTEventDisplacement[[k]] <- 0
	      lTEventSpeed[[k]]        <- 0
	      lTEventSinuosity[[k]]    <- 0
	      next #Skip To Next Event
	    }
	      
	      nEventsAnalysed               <- nEventsAnalysed + 1
	      nLarvalAnalysedFrames         <- nLarvalAnalysedFrames + length(datEventFrames$frameN)
	      ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
	      ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
	      ## Detects End Events as +ve differences
	      vDeltaXFrames        <- diff(datEventFrames$posX,lag=1,differences=1)
	      vDeltaYFrames        <- diff(datEventFrames$posY,lag=1,differences=1)
	      vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
	      
	        
	      dframe               <- diff(datEventFrames$frameN,lag=1,differences=1)
	      dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - Why is dFrame 0?  
	      #dEventSpeed          <- vEventPathLength/dframe
	      dEventSpeed          <- meanf(vEventPathLength/dframe,5) ##Apply Mean Filter Smooth Out 

	      if (any(is.nan(dEventSpeed)))
	      {
	        stop(paste("Speed is NaN, expID:",i,",eventID",k) )
	      }
	      ##Check If Some Video Frame Stuck Error Has Occured 
	      if (sum(dEventSpeed,na.rm=TRUE) == 0 & NROW(datEventFrames$frameN) > 2000)
	      {
	        warning(paste("Speed is 0 on ",NROW(datEventFrames$frameN),"frames video , expID:",i,",eventID",k) )
	        message((paste("Speed is 0 on ",NROW(datEventFrames$frameN),"frames video , expID:",i,",eventID",k) ))
	      }
	      ###                                              ######
        ### Process Speed/ Extract Bouts Via Peak Speed    ####
	      ###                                              ######
	      
	      ## TODO Add The Butterworth Filters Here filters ## 
	      dEventSpeed_smooth   <- meanf(dEventSpeed,20)
	      dEventSpeed_smooth[is.na(dEventSpeed_smooth)] = 0 ##Remove NA
	      
	      ### FInd Peaks In Speed to Identify Bouts
	      MoveboutsIdx <- find_peaks(dEventSpeed_smooth,25)
	      ##Reject Peaks Below Half An SD Peak Value - So As to Choose Only Significant Bout Movements # That Are Above the Minimum Speed to Consider As Bout
	      MoveboutsIdx_cleaned <- MoveboutsIdx[which(dEventSpeed_smooth[MoveboutsIdx] > sd(dEventSpeed_smooth[MoveboutsIdx])/2 
	                                                 & dEventSpeed_smooth[MoveboutsIdx] > G_MIN_BOUTSPEED   )  ]
	      
	      ## MOTION  BOUT Detection ##########
	      #nNumberOfBouts       <- length(MoveboutsIdx_cleaned)
	      #### Identify Bouts, Count, Duration 
	      ###Binarize , Use indicator function 1/0 for frames where Motion Occurs
	      vMotionBout <- dEventSpeed_smooth
	      vMotionBout[ vMotionBout < G_MIN_BOUTSPEED  ] = 0
	      vMotionBout[vMotionBout > G_MIN_BOUTSPEED  ] = 1
	      vMotionBout_OnOffDetect <- diff(vMotionBout) ##Set 1n;s on Onset, -1 On Offset of Bout
	      vMotionBout_On <- which(vMotionBout_OnOffDetect == 1)+1
	      if (any(is.nan(vMotionBout_On)) | any(is.na(vMotionBout_On)))
	      {
	        stop("vMotionBout_On Had Non Valid Entries")
	      }
	      vMotionBout_Off <- vMotionBout_On
	      if (NROW(vMotionBout_On) > 0) 
	        vMotionBout_Off <- which(vMotionBout_OnOffDetect[vMotionBout_On[1]:length(vMotionBout_OnOffDetect)] == -1)+vMotionBout_On[1] ##Ignore An Odd, Off Event Before An On Event, (ie start from after the 1st on event)
	      
	      if (NROW(vMotionBout_Off) == 0 ) ##No End FOund So Motion Bout Is Invalid? Or Does it just go to the end of this Frame
	      {
	        warning(paste("Motion Bout End Point not identified, removing Bout") )
	      }
	   
	      
	      nNumberOfBouts <- min(length(vMotionBout_On),length(vMotionBout_Off)) ##We can Only compare paired events, so remove an odd On Or Off Trailing Event
	   
	      
	         ##Remove The Motion Regions Where A Peak Was not detected / Only Keep The Bouts with Peaks 
	      vMotionBout[1:length(vMotionBout)] =0 ##Reset / Remove All Identified Movement
	      
	      for (i in 1:nNumberOfBouts)
	      {
	        if (nNumberOfBouts == 0)
	          break
	        
	        if (any(MoveboutsIdx_cleaned >= vMotionBout_On[i] & MoveboutsIdx_cleaned < vMotionBout_Off[i] ) == TRUE)
	        { ###Motion Interval Does not belong to a detect bout(peak) so remove
	          vMotionBout[vMotionBout_On[i]:vMotionBout_Off[i] ] = 1 ##Remove Motion From Vector
	        }else
	        {##Remove the Ones That Do not Have a peak In them
	          vMotionBout_On[i] = NA
	          vMotionBout_Off[i] = NA
	        }
	      }
	      vMotionBoutDuration <- vMotionBout_Off[1:nNumberOfBouts]-vMotionBout_On[1:nNumberOfBouts]
	      vMotionBoutDuration <- vMotionBoutDuration[!is.na(vMotionBoutDuration)]
	      ##### END OF MOTION BOUTS ###
	     
	      ## TODO TAIL Info ##
	      
	      ####
	      
	      ### MOTION Path Statistic ###
	      dEventTotalDistance       <- sum(vEventPathLength,na.rm=TRUE)
	      ##Straight Line From start to end 
	      dShortestPathDisplacement <- sqrt(((datEventFrames[1,]$posX-datEventFrames[NROW(datEventFrames),]$posX)^2+(datEventFrames[1,]$posY-datEventFrames[NROW(datEventFrames),]$posY)^2 ))
	      
	      ##
	      
	      ###(Actual Path Length Over Shortest Path) 
	      #SI < 1.05: almost straight
	      #1.05 ≤ SI <1.25: winding
	      #1.25 ≤ SI <1.50: twisty
	      #1.50 ≤ SI: meandering
	      if (dEventTotalDistance > 0)
	        dEventSinuosity          <- dEventTotalDistance/dShortestPathDisplacement 
	      else
	        dEventSinuosity <- 0
	      
	      if (any(is.nan(dEventSinuosity)))
	      {
	        message(paste("Sinuosity is NaN, expID:",i,",eventID",k, " Path Length:",dEventTotalDistance, " Shortest Path:",dShortestPathDisplacement) )
	        
	      }
	      
	      #stopifnot(any(is.nan(dEventSinuosity))==FALSE)
	     
	      
	      lTEventDisplacement[[k]] <- dEventTotalDistance
	      lTEventSpeed[[k]]        <- dEventSpeed
	      lTEventSinuosity[[k]]    <- dEventSinuosity

	      
  }## For each Event Of this Larva    
    
      nTotalFrames <- nTotalFrames + nLarvalAnalysedFrames
    
    
      ##Collect Results from All Events  of this Larva ##
      vTDisplacement    <- unlist(lTEventDisplacement)
      vTSpeed           <- unlist(lTEventSpeed)
      vTSinuosity       <- unlist(lTEventSinuosity)
    
    
    if (nLarvalAnalysedFrames > 0 )#(length(vTDisplacement) > 0)
    {
	    
      totalDisplacement <- sum(vTDisplacement)
      meanDisplacement  <- mean(vTDisplacement)
	    sdDisplacement   <-   sd(vTDisplacement)
	    medDisplacement  <- median(vTDisplacement)
	    maxDisplacement  <- max(vTDisplacement)
	    minDisplacement  <- min(vTDisplacement)
    	
	    ##Speed ##
	    meanSpeed <- mean(vTSpeed)
	    sdSpeed <-   sd(vTSpeed)
	    medSpeed  <- median(vTSpeed)
	    maxSpeed  <- max(vTSpeed)
	    minSpeed  <- min(vTSpeed)
	    ##Movement Ratio ##
	    ##Number of frames This Larva's Speed was above 0 versus total number of frames
	    
	    movementRatio <- ifelse(nLarvalAnalysedFrames==0,0, NROW(vTSpeed[vTSpeed > 0])/nLarvalAnalysedFrames)
	    
      #stopifnot(any(is.infinite(unlist(movementRatio))) == FALSE )
	    
	    stopifnot(is.nan(movementRatio) == FALSE )
	    stopifnot(is.infinite(movementRatio) == FALSE )
	    
	    ##Sinuosity / Curvature - Meander
	    meanSinuosity<- mean(vTSinuosity)
	    sdSinuosity <-   sd(vTSinuosity)
	    medSinuosity  <- median(vTSinuosity)
	    maxSinuosity  <- max(vTSinuosity)
	    minSinuosity  <- min(vTSinuosity)
	    
	
    }
    else
    {
      totalDisplacement = sdDisplacement = medDisplacement = maxDisplacement = minDisplacement = 0
      meanSpeed = sdSpeed = medSpeed = maxSpeed = minSpeed = 0
      meanSinuosity = sdSinuosity =  medSinuosity = maxSinuosity = minSinuosity=  0
	  }
    
  
    ## Colllect Larval Totals into Group Data
    ##+ Number of Hunting Events Counted as points where eyevergence events are more than 300 frames apart  ##
    lGroupMotion[[idx]] <- list(expID           =factor(i,levels=vexpID), 
                                   larvaID            =factor(larvaID,levels=seq(1:4)),
                                   dataSetID          = factor(DataSetID,levels=vdatasetID),
                                   totalframes        =length(datLarvaFrames$frameN),
                                   totalAnalysedframes=nLarvalAnalysedFrames,
                                   numberOfEvents     =nEventsAnalysed,
                                   vDisplacements     =unlist(vTDisplacement),
                                   totalDisplacement  =totalDisplacement,
                                   vSpeed             = unlist(vTSpeed),
                                   meanSpeed          =meanSpeed,
                                   sdSpeed            =meanSpeed,
                                   maxSpeed           =maxSpeed,
                                   medSpeed           =medSpeed,
                                   movementRatio      =movementRatio, 
                                   meanSinuosity      = meanSinuosity,
                                   sdSinuosity       = sdSinuosity,
                                   medSinuosity     = medSinuosity,
                                   maxSinuosity     =  maxSinuosity,
                                   minSinuosity    = minSinuosity
                               )

    
    
    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
  } ### For Each Larva In Group ####################
  ###############################################
  ### Collect LArva Stats to Make GROUP Stats ##
  ##############################################
  datGroupMotion = (do.call(rbind,lGroupMotion)) ##data.frame
#  datGroupHunting = as.data.frame(lGroupHunting)
  #vDisplacements = as.vector(do.call(rbind,lapply(datGroupMotion,"[[,","vDisplacements")))
  vDataSetID            <- unlist(datGroupMotion[,"dataSetID"])
  vexpID          <- unlist(datGroupMotion[,"expID"])
  vlarvaID          <- unlist(datGroupMotion[,"larvaID"])
  vDisplacements <- as.vector(unlist(datGroupMotion[,"vDisplacements"])) ##Distance Per Event
  vSpeed <- as.vector(unlist(datGroupMotion[,"vSpeed"]))
  vmeanSpeed <- as.vector(unlist(datGroupMotion[,"meanSpeed"])) ##Mean Speed Of Each Larva In Group
  #vSpeed = as.vector(do.call(rbind,lapply(datGroupMotion,"[[","vSpeed")))
  
  
  vEventCounts    <-  unlist(datGroupMotion[,"numberOfEvents"])
  nMotionEvents   <- sum(vEventCounts)
  vMovementRatio  <-  unlist(datGroupMotion[,"movementRatio"])
  vPathLengths    <- unlist(datGroupMotion[,"totalDisplacement"])
  vmuSinuosity    <- unlist(datGroupMotion[,"meanSinuosity"])
  vsdSinuosity    <- unlist(datGroupMotion[,"sdSinuosity"])
  
  nAnalysedFrames    <- sum(unlist(datGroupMotion[,"totalAnalysedframes"]))

  nMotionFrames    <- NROW(vSpeed[vSpeed > 0])
  groupMotionRatio <- ifelse(nTotalFrames==0,0, nMotionFrames/nTotalFrames)
  
  stopifnot(is.nan(groupMotionRatio) == FALSE)
  stopifnot(is.infinite(groupMotionRatio) == FALSE)
  
  meanMotionRatio   <- mean(vMovementRatio)
  sdMotionRatio     <- sd(vMovementRatio) 
  medmotionRatio    <- median(vMovementRatio)
  maxmotionRatio    <- max(vMovementRatio)
  minmotionRatio    <- min(vMovementRatio)
  
  ##Debug
  
  message("Motion Events For This Group:",nMotionEvents)

  totalEventCount <- sum(vEventCounts)
  meanEventCount  <- mean(vEventCounts)
  sdEventCount  <- sd(vEventCounts)
  
  ### Length Paths from All LArva Of This Group ###    
  meanGroupPathLength <- mean(vPathLengths)
  sdGroupPathLength   <- sd(vPathLengths)
  medGroupPathLength  <- median(vPathLengths)
  maxGroupPathLength <- max(vPathLengths)
  minGroupPathLength  <- min(vPathLengths)

  ### Speed From All Paths ###
  meanGroupSpeed <- mean(vSpeed)
  sdGroupSpeed   <- sd(vSpeed)
  medGroupSpeed  <- median(vSpeed)
  minGroupSpeed  <- min(vSpeed)
  maxGroupSpeed  <- max(vSpeed)

  ### Sinuosity ### Meander
  meanGroupSinuosity <- mean(vmuSinuosity)
  sdGroupSinuosity <- mean(sdSinuosity)

  if (nMotionEvents < 1)
   {
     warning("No Motion Events Detected in Group!")
      message("-* No Motion Events Detected in Group! *- ")
    #
    #meanHuntRatio       = 0
    #sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
    
    meanMotionRatio = 0;

    }
  ##Exp ID LookUp Table So Ican locate the same larva Across Empty->Live Condition
  ## 
  datLT <- data.frame(cbind(larvaID=levels(unlist(datGroupMotion$larvaID) )[unlist(datGroupMotion$larvaID)],
                            expID=levels(unlist(datGroupMotion$expID) )[unlist(datGroupMotion$expID)],
                            dataSetID=levels(unlist(datGroupMotion$dataSetID) )[unlist(datGroupMotion$dataSetID )]))
  udatLT <- unique(datLT)
  
  
    ###Return Results As List of Variable
    lGroupMotionStats <- list(nLarva=nLarva,
                           vDataSetID          = vDataSetID,
                           vIDLookupTable                = udatLT,
                            totalFrames   = nAnalysedFrames,
                            totalMotionFrames = nMotionFrames,
                            vEventCounts   = vEventCounts, ##Vectors Containing Means Per Larva In Group
                            vPathLengths   = vPathLengths,
                            vMovementRatio = vMovementRatio,
                            vSpeed         = vmeanSpeed, ##Mean Speed Of Each LArva
                            vSinuosity   = vmuSinuosity,#####
                            totalEventCount=totalEventCount,
                            meanEventCount= meanEventCount,
                            sdEventCount= sdEventCount,
                            seEventCount=sdEventCount/sqrt(nLarva),
                            motionRatio=groupMotionRatio,
                            meanMotionRatio=meanMotionRatio,##Mean among Larvae - Not the Mean Event
                            sdMotionRatio     = sdMotionRatio, 
                            medmotionRatio    = medmotionRatio,
                            maxmotionRatio    = maxmotionRatio,
                            minmotionRatio    = minmotionRatio,
                            meanPathLength    = meanGroupPathLength, ## Mean Of All Paths from The Group
                            sdPathLength = sdGroupPathLength,
                            sePathLength = sdGroupPathLength/sqrt(nLarva),
                            medPathLength = medGroupPathLength,
                            maxPathLength = maxGroupPathLength,
                            minPathLength = minGroupPathLength,
                            meanSpeed     = meanGroupSpeed,
                            sdSpeed     = sdGroupSpeed,
                            seSpeed     = sdGroupSpeed/sqrt(nLarva),
                            medSpeed  = medGroupSpeed,
                            minSpeed  = minGroupSpeed,
                            maxSpeed  = maxGroupSpeed,
                            meanSinuosity = meanGroupSinuosity,
                            sdSinuosity = sdGroupSinuosity,
                            seSinuosity = sdSinuosity/sqrt(nLarva)
                            )
  
    message(paste("Number of Motion Events detected:",totalEventCount, " Mean:", meanEventCount, " SD:",sdEventCount));
    message(paste(" Motion Ratio :",groupMotionRatio, " Event Mean:", meanMotionRatio, " SD:",sdMotionRatio));

    return (lGroupMotionStats)
}

















### Returns summary measurements for overall trajectories for each larva 
## Sinuoisity, Overall distance travelled, Speed (Avg,peak) etc
#########################################################

calcMotionStat2 <- function(datAllFrames,vexpID,vdatasetID)
{
  message("##Start Trajectory Analysis For Group##")
  lGroupMotion    <- list()
  lTEventDisplacement  <-list() ##Event Trajectory  Path Length
  lTEventSpeed         <-list() ##Event Trajectory  Speed
  lTEventSinuosity     <- list()
  
  idx             <- 0
  nLarva  	      <- length(vexpID)
  nEventsAnalysed <- 0
  nTotalFrames    <- 0
  
  for (i in vexpID)
  {
    idx                 <- idx+1; ## Only increment if Record Is to be added
    
    datLarvaFrames      <- datAllFrames[datAllFrames$expID == i,]
    vEventID            <- unique((datLarvaFrames$eventID))
    ##Reset Variables ##
    nEventsAnalysed        <- 0
    nLarvalAnalysedFrames  <- 0 ## Number of Motion Frames  Analysed For this Larva
    movementRatio          <- 0
    
    
    DataSetID             <- ifelse(any(names(datLarvaFrames) == "dataSet"),unique(datLarvaFrames$dataSet),0 )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID )
    larvaID               <- unique((datLarvaFrames$larvaID))
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFrames[datLarvaFrames$REyeAngle > -35 & datLarvaFrames$REyeAngle <  35 & datLarvaFrames$LEyeAngle > -35 & datLarvaFrames$LEyeAngle <  35,]
    
    
    lTEventDisplacement  <-list() ##Event Trajectory  Path Length
    lTEventSpeed         <-list() ##Event Trajectory  Speed
    lTEventSinuosity     <- list()
    
    vTDisplacement    <- list()
    vTSpeed           <- list()
    vTSinuosity       <- list()
    
    
    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    for (k in vEventID)
    {
      ##Error Captcha ##
      if (is.na(k) | k < 1 )
      {
        stop(paste("Invalid Event Idx:",k) )
      }
      
      ##Select Motion Frames of this Event ## 
      ## Filter Out Tracking Loses Set to poxX,Y=0
      datEventFrames       <- datLarvaFrames[datLarvaFrames$eventID == k & datLarvaFrames$posX > 0 & datLarvaFrames$posY > 0,]
      
      
      ##Consider Only if there are at least 100 frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
      if  (NROW(datEventFrames$frameN) > 100)
      {
        
        nEventsAnalysed               <- nEventsAnalysed + 1
        nLarvalAnalysedFrames         <- nLarvalAnalysedFrames + length(datEventFrames$frameN)
        ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
        ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
        ## Detects End Events as +ve differences
        vDeltaXFrames        <- diff(datEventFrames$posX,lag=1,differences=1)
        vDeltaYFrames        <- diff(datEventFrames$posY,lag=1,differences=1)
        vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2)
        dframe               <- diff(datEventFrames$frameN,lag=1,differences=1)
        dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - Why is dFrame 0?  
        dEventSpeed          <- vEventPathLength/dframe
        
        if (any(is.nan(dEventSpeed)))
        {
          stop(paste("Speed is NaN, expID:",i,",eventID",k) )
        }
        ##Check If Some Video Frame Stuck Error Has Occured 
        if (sum(dEventSpeed,na.rm = TRUE) == 0 & NROW(datEventFrames$frameN) > 2000)
        {
          warning(paste("Speed is 0 on ",NROW(datEventFrames$frameN),"frames video , expID:",i,",eventID",k) )
          message((paste("Speed is 0 on ",NROW(datEventFrames$frameN),"frames video , expID:",i,",eventID",k) ))
        }
        
        
        dEventTotalDistance       <- sum(vEventPathLength)
        ##Straight Line From start to end 
        dShortestPathDisplacement <- sqrt(((datEventFrames[1,]$posX-datEventFrames[NROW(datEventFrames),]$posX)^2+(datEventFrames[1,]$posY-datEventFrames[NROW(datEventFrames),]$posY)^2 ))
        
        ###(Actual Path Length Over Shortest Path) 
        #SI < 1.05: almost straight
        #1.05 ≤ SI <1.25: winding
        #1.25 ≤ SI <1.50: twisty
        #1.50 ≤ SI: meandering
        dEventSinuosity          <- dEventTotalDistance/dShortestPathDisplacement 
        
        if (any(is.nan(dEventSinuosity)))
        {
          message(paste("Sinuosity is NaN, expID:",i,",eventID",k) )
          
        }
        stopifnot(any(is.nan(dEventSinuosity))==FALSE)
        
        
        lTEventDisplacement[[k]] <- dEventTotalDistance
        lTEventSpeed[[k]]        <- dEventSpeed
        lTEventSinuosity[[k]]    <- dEventSinuosity
        
      }else{
        ## Episode Trajectory is too Short Insufficient ##	
        warning(paste("Trajectory Too short - Ignoring Trajectory  of Experiment:",i,"event",k," larva",larvaID) )
        
        lTEventDisplacement[[k]] <- 0
        lTEventSpeed[[k]]        <- 0
        lTEventSinuosity[[k]]    <- 0
      }
      
    }##For each Event Of this Larva    
    
    nTotalFrames <- nTotalFrames + nLarvalAnalysedFrames
    
    
    ##Collect Results from All Events  of this Larva ##
    vTDisplacement    <- unlist(lTEventDisplacement)
    vTSpeed           <- unlist(lTEventSpeed)
    vTSinuosity       <- unlist(lTEventSinuosity)
    
    
    if (nLarvalAnalysedFrames > 0 )#(length(vTDisplacement) > 0)
    {
      
      totalDisplacement <- sum(vTDisplacement)
      meanDisplacement  <- mean(vTDisplacement)
      sdDisplacement   <-   sd(vTDisplacement)
      medDisplacement  <- median(vTDisplacement)
      maxDisplacement  <- max(vTDisplacement)
      minDisplacement  <- min(vTDisplacement)
      
      ##Speed ##
      meanSpeed <- mean(vTSpeed)
      sdSpeed <-   sd(vTSpeed)
      medSpeed  <- median(vTSpeed)
      maxSpeed  <- max(vTSpeed)
      minSpeed  <- min(vTSpeed)
      ##Movement Ratio ##
      ##Number of frames This Larva's Speed was above 0 versus total number of frames
      
      movementRatio <- ifelse(nLarvalAnalysedFrames==0,0, NROW(vTSpeed[vTSpeed > 0])/nLarvalAnalysedFrames)
      
      #stopifnot(any(is.infinite(unlist(movementRatio))) == FALSE )
      
      stopifnot(is.nan(movementRatio) == FALSE )
      stopifnot(is.infinite(movementRatio) == FALSE )
      
      ##Sinuosity / Curvature - Meander
      meanSinuosity<- mean(vTSinuosity)
      sdSinuosity <-   sd(vTSinuosity)
      medSinuosity  <- median(vTSinuosity)
      maxSinuosity  <- max(vTSinuosity)
      minSinuosity  <- min(vTSinuosity)
      
      
    }
    else
    {
      totalDisplacement = sdDisplacement = medDisplacement = maxDisplacement = minDisplacement = 0
      meanSpeed = sdSpeed = medSpeed = maxSpeed = minSpeed = 0
      meanSinuosity = sdSinuosity =  medSinuosity = maxSinuosity = minSinuosity=  0
    }
    
    
    ## Colllect Larval Totals into Group Data
    ##+ Number of Hunting Events Counted as points where eyevergence events are more than 300 frames apart  ##
    lGroupMotion[[idx]] <- list(expID           =factor(i,levels=vexpID), 
                                larvaID            =factor(larvaID,levels=seq(1:4)),
                                dataSetID          = factor(DataSetID,levels=vdatasetID),
                                totalframes        =length(datLarvaFrames$frameN),
                                totalAnalysedframes=nLarvalAnalysedFrames,
                                numberOfEvents     =nEventsAnalysed,
                                vDisplacements     =unlist(vTDisplacement),
                                totalDisplacement  =totalDisplacement,
                                vSpeed             = unlist(vTSpeed), ##Displacement / Nframes - this is also used for mean
                                meanSpeed          =meanSpeed,
                                sdSpeed            =meanSpeed,
                                maxSpeed           =maxSpeed,
                                medSpeed           =medSpeed,
                                movementRatio      =movementRatio, 
                                meanSinuosity      = meanSinuosity,
                                sdSinuosity       = sdSinuosity,
                                medSinuosity     = medSinuosity,
                                maxSinuosity     =  maxSinuosity,
                                minSinuosity    = minSinuosity
    )
    
    
    
    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
  } ### For Each Larva In Group ####################
  ###############################################
  ### Collect LArva Stats to Make GROUP Stats ##
  ##############################################
  datGroupMotion = do.call(rbind,lGroupMotion)
  #  datGroupHunting = as.data.frame(lGroupHunting)
  #vDisplacements = as.vector(do.call(rbind,lapply(datGroupMotion,"[[,","vDisplacements")))
  vDataSetID            <- unlist(datGroupMotion[,"dataSetID"])
  vexpID          <- unlist(datGroupMotion[,"expID"])
  vlarvaID          <- unlist(datGroupMotion[,"larvaID"])
  vDisplacements <- as.vector(unlist(datGroupMotion[,"vDisplacements"])) ##Distance Per Event
  vSpeed <- as.vector(unlist(datGroupMotion[,"vSpeed"]))
  vmeanSpeed <- as.vector(unlist(datGroupMotion[,"meanSpeed"])) ##Mean Speed Of Each Larva In Group
  #vSpeed = as.vector(do.call(rbind,lapply(datGroupMotion,"[[","vSpeed")))
  
  
  vEventCounts    <-  unlist(datGroupMotion[,"numberOfEvents"])
  nMotionEvents   <- sum(vEventCounts)
  vMovementRatio  <-  unlist(datGroupMotion[,"movementRatio"])
  vPathLengths    <- unlist(datGroupMotion[,"totalDisplacement"])
  vmuSinuosity    <- unlist(datGroupMotion[,"meanSinuosity"])
  vsdSinuosity    <- unlist(datGroupMotion[,"sdSinuosity"])
  
  nAnalysedFrames    <- sum(unlist(datGroupMotion[,"totalAnalysedframes"]))
  
  nMotionFrames    <- NROW(vSpeed[vSpeed > 0])
  groupMotionRatio <- ifelse(nTotalFrames==0,0, nMotionFrames/nTotalFrames)
  
  stopifnot(is.nan(groupMotionRatio) == FALSE)
  stopifnot(is.infinite(groupMotionRatio) == FALSE)
  
  meanMotionRatio   <- mean(vMovementRatio)
  sdMotionRatio     <- sd(vMovementRatio) 
  medmotionRatio    <- median(vMovementRatio)
  maxmotionRatio    <- max(vMovementRatio)
  minmotionRatio    <- min(vMovementRatio)
  
  ##Debug
  
  message("Motion Events For This Group:",nMotionEvents)
  
  totalEventCount <- sum(vEventCounts)
  meanEventCount  <- mean(vEventCounts)
  sdEventCount  <- sd(vEventCounts)
  
  ### Length Paths from All LArva Of This Group ###    
  meanGroupPathLength <- mean(vPathLengths)
  sdGroupPathLength   <- sd(vPathLengths)
  medGroupPathLength  <- median(vPathLengths)
  maxGroupPathLength <- max(vPathLengths)
  minGroupPathLength  <- min(vPathLengths)
  
  ### Speed From All Paths ###
  meanGroupSpeed <- mean(vSpeed)
  sdGroupSpeed   <- sd(vSpeed)
  medGroupSpeed  <- median(vSpeed)
  minGroupSpeed  <- min(vSpeed)
  maxGroupSpeed  <- max(vSpeed)
  
  ### Sinuosity ### Meander
  meanGroupSinuosity <- mean(vmuSinuosity)
  sdGroupSinuosity <- mean(sdSinuosity)
  
  if (nMotionEvents < 1)
  {
    warning("No Motion Events Detected in Group!")
    message("-* No Motion Events Detected in Group! *- ")
    #
    #meanHuntRatio       = 0
    #sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
    
    meanMotionRatio = 0;
    
  }
  ##Exp ID LookUp Table So Ican locate the same larva Across Empty->Live Condition
  ## 
  datLT <- data.frame(cbind(larvaID=levels(datGroupMotion$larvaID)[datGroupMotion$larvaID],expID=levels(datGroupMotion$expID)[datGroupMotion$expID],dataSetID=levels(datGroupMotion$dataSetID)[datGroupMotion$dataSetID]))
  udatLT <- unique(datLT)
  
  
  ###Return Results As List of Variable
  lGroupMotionStats <- list(nLarva=nLarva,
                            vDataSetID          = vDataSetID,
                            vIDLookupTable                = udatLT,
                            totalFrames   = nAnalysedFrames,
                            totalMotionFrames = nMotionFrames,
                            vEventCounts   = vEventCounts, ##Vectors Containing Means Per Larva In Group
                            vPathLengths   = vPathLengths,
                            vMovementRatio = vMovementRatio,
                            vSpeed         = vmeanSpeed, ##Mean Speed Of Each LArva
                            vSinuosity   = vmuSinuosity,#####
                            totalEventCount=totalEventCount,
                            meanEventCount= meanEventCount,
                            sdEventCount= sdEventCount,
                            seEventCount=sdEventCount/sqrt(nLarva),
                            motionRatio=groupMotionRatio,
                            meanMotionRatio=meanMotionRatio,##Mean among Larvae - Not the Mean Event
                            sdMotionRatio     = sdMotionRatio, 
                            medmotionRatio    = medmotionRatio,
                            maxmotionRatio    = maxmotionRatio,
                            minmotionRatio    = minmotionRatio,
                            meanPathLength    = meanGroupPathLength, ## Mean Of All Paths from The Group
                            sdPathLength = sdGroupPathLength,
                            sePathLength = sdGroupPathLength/sqrt(nLarva),
                            medPathLength = medGroupPathLength,
                            maxPathLength = maxGroupPathLength,
                            minPathLength = minGroupPathLength,
                            meanSpeed     = meanGroupSpeed,
                            sdSpeed     = sdGroupSpeed,
                            seSpeed     = sdGroupSpeed/sqrt(nLarva),
                            medSpeed  = medGroupSpeed,
                            minSpeed  = minGroupSpeed,
                            maxSpeed  = maxGroupSpeed,
                            meanSinuosity = meanGroupSinuosity,
                            sdSinuosity = sdGroupSinuosity,
                            seSinuosity = sdSinuosity/sqrt(nLarva)
  )
  
  message(paste("Number of Motion Events detected:",totalEventCount, " Mean:", meanEventCount, " SD:",sdEventCount));
  message(paste(" Motion Ratio :",groupMotionRatio, " Event Mean:", meanMotionRatio, " SD:",sdMotionRatio));
  
  return (lGroupMotionStats)
} ###Version 2 Of Function ##





##Filter Files With Durable EyeVergences above N number of frames
#vframePerFile = table(datHuntFrames$fileIdx)
#vfileIdx  = as.numeric(names(vframePerFile[vframePerFile > 300 ]))
#vtargetFiles = vfileIdx;

# 
# dThetaREye = diff(filtereddatAllFrames$REyeAngle[filtereddatAllFrames$fileIdx %in% vtargetFiles ],differences=1)
# dThetaLEye = diff(filtereddatAllFrames$LEyeAngle[filtereddatAllFrames$fileIdx %in% vtargetFiles],differences=1)
# 
# filteredREye = convolve(filtereddatAllFrames$REyeAngle[filtereddatAllFrames$fileIdx %in% vtargetFiles],fstepRDetect,type="filter")
# filteredLEye = convolve(filtereddatAllFrames$LEyeAngle[filtereddatAllFrames$fileIdx %in% vtargetFiles],fstepLDetect,type="filter")
# filteredFrames = filtereddatAllFrames$frameN[filtereddatAllFrames$fileIdx %in% vtargetFiles]

## Hunting Events of file 161 Based On R Eye Only ##
# cfilterThreshold = 2000;
# temp[vtargetFiles]
# huntFrames <- filteredFrames[ filteredREye > cfilterThreshold & filteredLEye > cfilterThreshold];
# 
# for (i in vtargetFiles)
# {
#   message(paste("Plot File: ",temp[i]) )
#   
#   expName = strsplit (temp[i],'/')
#   expName = file_path_sans_ext(expName[[1]][length(expName[[1]])])
#   
#   
#   xl <- filtereddatAllFrames$frameN[filtereddatAllFrames$fileIdx == i & filtereddatAllFrames$frameN %in% huntFrames]
#   yl <- filtereddatAllFrames$LEyeAngle[filtereddatAllFrames$fileIdx == i & filtereddatAllFrames$frameN %in% huntFrames]
#   xr <- filtereddatAllFrames$frameN[filtereddatAllFrames$fileIdx == i & filtereddatAllFrames$frameN %in% huntFrames]
#   yr <- filtereddatAllFrames$REyeAngle[filtereddatAllFrames$fileIdx == i & filtereddatAllFrames$frameN %in% huntFrames]
#   
#   if (length(xl) < 10 | length(xr) < 10 | length(yl) < 10 | length(yr) < 10 )
#   {
#     #
#     message(c(length(xl), length(yl),length(xr),length(yr)) )
#     next
#   }
#   
#   setEPS();
#   postscript(paste("./plots/timeseries/",expName,".eps"),width=8,height=8)
#   plot(xl, yl, xlab="Frame #",ylab="Eye Angle", main=expName, pch=24,ylim = c(-40,40), col='green' )
#   lines(filtereddatAllFrames$frameN[filtereddatAllFrames$fileIdx == i],filtereddatAllFrames$LEyeAngle[filtereddatAllFrames$fileIdx == i])
#   
#   points(xr, yr, xlab="Frame #",ylab="R Eye Angle",col='red',pch=25)
#   lines(filtereddatAllFrames$frameN[filtereddatAllFrames$fileIdx == i],filtereddatAllFrames$REyeAngle[filtereddatAllFrames$fileIdx == i])
#   title(expName)
#   dev.off()
#   
#   #break
# }
