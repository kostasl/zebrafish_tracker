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


calcMotionStat <- function(datAllFrames,vexpID,vdatasetID)
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
	      vEventPathLength     <- sqrt(vDeltaXFrames^2+vDeltaYFrames^2) ## Path Length Calculated As Total Displacement
	      #nNumberOfBouts       <- 
	      dframe               <- diff(datEventFrames$frameN,lag=1,differences=1)
	      dframe               <- dframe[dframe > 0] ##Clear Any possible Nan - Why is dFrame 0?  
	      dEventSpeed          <- vEventPathLength/dframe
	      
	      if (any(is.nan(dEventSpeed)))
	      {
	        stop(paste("Speed is NaN, expID:",i,",eventID",k) )
	      }
	      ##Check If Some Video Frame Stuck Error Has Occured 
	      if (sum(dEventSpeed) == 0 && NROW(datEventFrames$frameN) > 2000)
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
}


















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
        if (sum(dEventSpeed) == 0 && NROW(datEventFrames$frameN) > 2000)
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
