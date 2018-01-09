# ################
# ## \todo Code to Detect Prey Hunting Events - Locate Videos And Start Frames - ###
# ## Rank Files With Most Hunting

# M. Meyer sent doc with ideas on 30-11-17, suggests:
# * Some additional analyses
# 1. Number of tracks containing vergence events/total number of tracks= probability that a group will switch into hunting mode
# 2. For all tracks with eye vergence events calculate the average duration of vergence. 
#KL:
# Since nLarva can vary, I suggest looking at statistics of HuntRates and Durations as per fish in Group

#################

##Find Corresponding Video File to Track File###
getVideofilePath <- function(x,strVideoFilePath,init="Auto",end="_tracks",postfix=".mp4"){
  z=regexpr('Auto.*tracks',x,T)
  
  z<-paste(substr(x,z,z+attr(z,"match.length")-1-nchar(end)),postfix,sep="")
  ##In lack of a better method to locate fullpaths of a filename list 
  for (f in z)
  {
    ##Replace with Full Path
    message(f)
    strfilepath <- list.files(path =strVideoFilePath, pattern = f, all.files = FALSE,
                              full.names = TRUE, recursive = TRUE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    message(strfilepath)
    if (length(strfilepath) > 0)
    {
      #z[which(z==f)] <- strfilepath  ##Replace With Full Path If Found
      z[which(z==f)] <- basename(path = strfilepath)  ##Replace With File Name Only If Found (We search Again Upon Labelling)
    }
    else
    {
      warning(paste("FILE NOT FOUND :",strfilepath) )
    }
  }
  
  
  
  return(z)
}

calcHuntStat <- function(datAllFrames,vlarvaID)
{
  message("##Start Hunting Analysis For Group##")
  lGroupHunting    <- list()
  lHuntingEvents    <- list()
  #lHuntingDuration <-list()
  #vHuntStartFrames <- list()  
  #vHuntEndFrames <- list()  

  try(rm("vHuntStartFrames","vHuntEndFrames","vHuntDeltaFrames","lHuntingDuration"),silent=TRUE)
  
  
  nHuntingEventsForLarva   <- 0
  meanHuntingEvents <- 0;
  stdDevHuntingEvents <- 0;
  idx = 0;
  idxHuntRec = 0;
  nLarva  	<- length(vlarvaID)
  

  for (i in vlarvaID)
  {
    idx = idx+1;
    datLarvaFrames <- datAllFrames[datAllFrames$larvaID == i,]

    #Used to Identify Experimet Set From Which Data COmes from - Plotted AS different colour
    DataSetID             <- ifelse(any(names(datLarvaFrames) == "dataSet"),unique(datLarvaFrames$dataSet),0 )
    stopifnot(DataSetID >= 0 )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID ) 
    
    vEventID = unique((datLarvaFrames$eventID))
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFrames[datLarvaFrames$REyeAngle > -35 & datLarvaFrames$REyeAngle <  35 & datLarvaFrames$LEyeAngle > -35 & datLarvaFrames$LEyeAngle <  35,]
    
    nTotalHuntFrames = 0;
    nHuntingEventsForLarva <- 0;
    meanHuntDuration <- sdHuntDuration  <- medHuntDuration <- maxHuntDuration<-minHuntDuration <- 0
    
    lHuntingDuration <-list() ##Empty List Of Event Durations
    
    
    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    for (k in vEventID)
    {

	    ##Select Hunt Frames of this Event ## 
      ## Note Criterion - Left Eye Angle >  Min & Right Eye Angle <  Min (Turned Inwards) & Vergence Angle > 40
      ## Warning - This criterion Is repeated in the plot Function of plotTrackScatterAndDensities - Warning 
	    datHuntFrames    <- datLarvaFrames[datLarvaFrames$eventID == k & 
	                                         datLarvaFrames$REyeAngle < -G_THRESHUNTANGLE &
	                                         datLarvaFrames$LEyeAngle > G_THRESHUNTANGLE & 
	                                         abs(datLarvaFrames$LEyeAngle-datLarvaFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE,]
      nTotalHuntFrames <- nTotalHuntFrames + length(datHuntFrames$frameN)
	    ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
	    ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
	    ## Detects End Events as +ve differences
	    vHuntDeltaFrames <- diff(datHuntFrames$frameN,lag=1,differences=1)
		
	    #Find Start Hunt Event the point leading a big gap in frameN, ie reverse diff s[(1+1):10] - s[(1):9]
	    ##ShifHuntDelteFrames
	    vsHuntDeltaFrames     <- (1:(NROW(datHuntFrames$frameN)+2))
	    vHuntStartFrames      <- list()
	    vHuntEndFrames        <- list()
	    lHuntingDuration[[k]] <- 0
	    minHuntDuration       <- 0
	    
	    
	    # Debug #
	    #message(paste(i,".",k,". Detected Hunt #Frames=", NROW(datHuntFrames$frameN) ))
	            
	    
	    ##Consider Only if there are at least 100 frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
	    if  (NROW(datHuntFrames$frameN) > G_MINEPISODEDURATION)
	    {
	      ##Add (Imaginary) Edge Frame  Numbers  (preceding hunting and after last hunting Frame)
	      vsHuntDeltaFrames[1] <- datHuntFrames$frameN[1]-G_MINGAPBETWEENEPISODES
	      vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- datHuntFrames$frameN[NROW(datHuntFrames$frameN)]+G_MINGAPBETWEENEPISODES
		    ##Copy into Shifted right (lagged) Position
	 	    vsHuntDeltaFrames[2:(NROW(datHuntFrames$frameN)+1)] <- datHuntFrames$frameN
		    ##Do Rev Diff - ie Nabla taking s_n = x_n-x_{n-1} # Detect Start Events as +ve diff
		    vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] <- datHuntFrames$frameN - vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN))]
		    ##\TODO:  Perhaps Throw away Hunts That Began Before Recording Started or End After It Ends - remove Edge Points ##

		    ##Find Hunt Event STARTs ## Needs A gap of At least #MINGAP Between Events To Be a new Episode
		    ##attempting to Remove delay introduced by filtering 
		    vHuntStartFrames <- datHuntFrames$frameN[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
		    ##Find Hunt Event ENDs - Needs To Be atleast G_MINGAPBETWEENEPISODES Away from next start
		    ##Shift (expanded) DFrame vector to pickup on adjacent Hunt-Frame number signifying the end of the Hunt Episode
		    vsHuntDeltaFrames[2:NROW(vsHuntDeltaFrames)] <- vsHuntDeltaFrames[1:NROW(vsHuntDeltaFrames)-1]
        ##Pickout End frames by selected frameN locations shifted by 1+
		    vHuntEndFrames  <- datHuntFrames$frameN[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
		    ###

		    #message(paste("Found i=",length(vHuntEndFrames)," Hunt End points and k=",length(vHuntStartFrames), "Start Points"))
		    stopifnot(length(vHuntEndFrames) == length(vHuntStartFrames))
		    lHuntingDuration[[k]]   <- vHuntEndFrames - vHuntStartFrames
		    message(paste(" Event Hunt Duration is  ",sum(lHuntingDuration[[k]])," in ",length(vHuntStartFrames)," episodes for larvaID:",i," eventID:",k))
		    
		    # idxHuntRec = idxHuntRec + 1;
		    # lHuntingEvents[[idxHuntRec]] <- list(larvaID =rep(i,length(vHuntEndFrames)),
		    #                                      eventID = rep(k,length(vHuntEndFrames)),
		    #                                      dataSetID = rep(DataSetID,length(vHuntEndFrames)),
		    #                                      fileIdx = rep(unique(datHuntFrames$fileIdx),length(vHuntEndFrames)),
		    #                                      startFrame=unlist(vHuntStartFrames),
		    #                                      endFrame=unlist(vHuntEndFrames),
		    #                                      huntScore=rep(0,length(vHuntEndFrames))) ##0 To Be used for Manual Labelling Of Hunting Success
		    # 
		    idxHuntRec = idxHuntRec + 1;
		    lHuntingEvents[[idxHuntRec]] <- data.frame(larvaID =i,
                                   eventID = k,
                                   dataSetID = DataSetID,
                                   fileIdx = unique(datHuntFrames$fileIdx),
                                   startFrame=unlist(vHuntStartFrames),
                                   endFrame=unlist(vHuntEndFrames),
                                   huntScore=0) ##0 To Be used for Manual Labelling Of Hunting Success

		    if( lHuntingDuration[[k]] < 0) 
		      stop(paste("Negative Hunt Duration detected LarvaID:",i," eventID:",k," Possible Duplicate tracker file" ) ) ##Catch Error / Can be caused by Duplicate Tracked Video - Check )
		    
		    minHuntDuration  <- min(lHuntingDuration[[k]])
		    
		    nHuntingEventsForLarva  <- nHuntingEventsForLarva + length(vHuntStartFrames) ##increment Event Count
	    }else{
       ##Hunting Episode Insufficient ##	
       # vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN)+1)] <- 0  ##Copy into Shifted Position
            warning(paste("Ignoring low Hunt Frame Number",NROW(datHuntFrames$frameN)," of Larva:",i,"event",k) )
            if (is.na(k) | k < 1 )
            {
              stop(paste("Invalid Event Idx:",k) )
            }
	      minHuntDuration       <-0
	      lHuntingDuration[[k]] <- 0 ##Otherwise Unused Event IDs end up NULL
	    }
	    
	    
	    #### ERROR Hanlders /  DEbug ###
	    
	    if (is.na(minHuntDuration)) 
	    {  
	      message(paste("Warning - Min Episode Duration is  ",minHuntDuration," for larvaID:",i," eventID:",k))
	      warning(paste("Warning - Min Episode Duration is  ",minHuntDuration," for larvaID:",i," eventID:",k))
      # message("** Probably dublicate track files present **")
	    }else if (minHuntDuration < 0) 
	    {  
	      
	      message(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
	      warning(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
	      message("** Probably dublicate track files present **")
	    }
	   
	    if( any( is.null( lHuntingDuration ) ) )  
	      stop(paste("Null Hunting duration for larvaID:",i," eventID:",k,". Possible Missing Event IDs",  as.character( unlist(which(sapply(lHuntingDuration, is.null) == TRUE)) )) )

     }##For each Event Of this Larva  ##########  
  #####
    vHuntingDuration = unlist(lHuntingDuration) ##Vector With All Event Durations
    if (length(vHuntingDuration) > 0)
    {
	    meanHuntDuration <- mean(vHuntingDuration)
	    sdHuntDuration <-   sd(vHuntingDuration)
	    medHuntDuration  <- median(vHuntingDuration)
	    maxHuntDuration  <- max(vHuntingDuration)
	    minHuntDuration  <- min(vHuntingDuration)
    	##Debug
    }
    else
    {
      meanHuntDuration = medHuntDuration = sdHuntDuration = maxHuntDuration = minHuntDuration = 0
    }
  
    if (is.nan(meanHuntDuration) | meanHuntDuration < 0)
    {
      stop(paste("Invalid Hunt Duration for LarvaId:",i," eventID:",k) )
    }
    
    if (meanHuntDuration == 0)
    {
      warning(paste("ZERO Hunt Duration for LarvaId:",i," All events" ) )
      message(paste("ZERO Hunt Duration for LarvaId:",i," All events") )
    }
    
    ##Total Duration of Hunting For This Larva ##
    ##+ Number of Hunting Events Counted as points where eyevergence events are more than 300 frames apart  ##
    lGroupHunting[[idx]] = list(larvaID=i,
                                dataSetID=DataSetID,
                                count=nHuntingEventsForLarva,
                                huntframes=nTotalHuntFrames,
                        				vHuntDurations=vHuntingDuration, ## Hunt Duration of each Episode  (All data)    	                        
                        				meanHDuration=meanHuntDuration, ##Mean Duration for this Larva (per Episode)
                        				sdHDuration=sdHuntDuration,
                        				seHDuration=sdHuntDuration/sqrt(nLarva),
                                medHDuration=medHuntDuration,
                  			        maxHDuration = maxHuntDuration,
                  			        minHDuration = minHuntDuration,
                                totalframes=length(datLarvaFrames$frameN)
				)
    
    
    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
  } ### For Each Larva In Group ###
##### Note That vHuntStartFrames is invalid beyond this point #

  datGroupHunting  = do.call(rbind,lGroupHunting)
  datHuntingEvents = do.call(rbind,lHuntingEvents)
#  datGroupHunting = as.data.frame(lGroupHunting)
  ntotalRecFrames               <- unlist(datGroupHunting[,"totalframes"])
  vHuntEpisodeDuration          <- as.vector(do.call(rbind,lapply(lGroupHunting,"[[","vHuntDurations"))) ## Of All Episodes
  
  vHuntRatio                    <- unlist(datGroupHunting[,"huntframes"])/ntotalRecFrames 
  ##Where LArva Did not Hunt Then Set Ratio to 0 SO this shows up in data
  vHuntRatio[is.nan(vHuntRatio)== TRUE] <- 0 
  #vHuntRatio <- vHuntRatio[is.nan(vHuntRatio)== FALSE]  #Remove Nan Entries So The mean Is not Affected / Produced When Total Frames Is 0 for some LArva ID
  
  
  vHuntingDuration          <- unlist(datGroupHunting[,"huntframes"])##Vector with Total for each Larva
  vmeanLarvaEpisodeHuntingDuration <- unlist(datGroupHunting[,"meanHDuration"]) ## Vector with Mean Episode Duration of eacg Larva
  vnHuntingEvents           <- unlist(datGroupHunting[,"count"])
  vDataSetID                <- unlist(datGroupHunting[,"dataSetID"])
  nHuntingDuration          <- sum(vHuntingDuration)
  nHuntingEvents            <- sum(vnHuntingEvents);
  ntotalFrames              <- sum(unlist(datGroupHunting[,"totalframes"]))
  meanHuntRatio      <- mean(vHuntRatio) ##Group
  sdHuntRatio        <- sd(vHuntRatio) 
  medHuntRatio       <- median(vHuntRatio)
  maxHuntRatio       <- max(vHuntRatio)
  minHuntRatio       <- min(vHuntRatio)
  
  stopifnot(any(meanHuntRatio >= 0 ) )
  ##Debug
  message("Hunting Event For This Group:",nHuntingEvents)

    ### Hunt Duration Per Larva ###    
    meanLarvaHuntDuration <- mean(vHuntingDuration)
    sdLarvaHuntDuration   <- sd(vHuntingDuration)
    medLarvaHuntDuration  <- median(vHuntingDuration)
    maxLarvaHuntDuration  <- max(vHuntingDuration)
    minLarvaHuntDuration  <- min(vHuntingDuration)
	
    ### Hunt Duration Per Hunting Episode ###
    meanHuntEpisodeDuration <- mean(vHuntEpisodeDuration)
    sdHuntEpisodeDuration   <- sd(vHuntEpisodeDuration)
    medHuntEpisodeDuration  <- median(vHuntEpisodeDuration)
    minHuntEpisodeDuration  <- min(vHuntEpisodeDuration)
    maxHuntEpisodeDuration  <- max(vHuntEpisodeDuration)



    ### Number of Hunting Events Per Larva ##
    meanHuntingEvents  <- mean(vnHuntingEvents);
    sdDevHuntingEvents <- sd(vnHuntingEvents);
    medHuntingEvents   <- median(vnHuntingEvents);
    maxHuntingEvents   <- max(vnHuntingEvents);
    minHuntingEvents   <- min(vnHuntingEvents);
   if (nHuntingEvents > 1)
  {
      
  }
  else
  {
    warning("No Hunting Events Detected in Group!")
    message("-* No Hunting Events Detected in Group! *- ")
    #
    #meanHuntRatio       = 0
    #sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
    
    meanHuntingEvents = 0;
    sdDevHuntingEvents = 0;
    
  }
  
   
    
    ###Return Results As List of Variable
    lGroupHuntStats <- list(nLarva=nLarva,
                            totalFrames         =  ntotalFrames,
                            vDataSetID          = vDataSetID,
                            vLarvaID            = vlarvaID,
                            vHDurationPerLarva = vHuntingDuration,
                            vmeanHEpisodeDurationPerLarva	= vmeanLarvaEpisodeHuntingDuration,
                            vHLarvaEventCount   = vnHuntingEvents,
                            vLarvaHRatio        = vHuntRatio,
                            huntFrames=nHuntingDuration,
                            meanDuration=meanLarvaHuntDuration,
                  			    sdDuration=sdLarvaHuntDuration,
                  			    seDuration=sdHuntDuration/sqrt(nLarva),
                  			    medDuration             =medLarvaHuntDuration,
                  			    maxDuration             = maxLarvaHuntDuration,
                  			    minDuration             = minLarvaHuntDuration,
                  			    meanEpisodeDuration=meanHuntEpisodeDuration,
                  			    sdEpisodeDuration=sdHuntEpisodeDuration,
                  			    seEpisodeDuration=sdHuntDuration/sqrt(nLarva),
                  			    medEpisodeDuration=medHuntEpisodeDuration,
                  			    maxEpisodeDuration = maxHuntEpisodeDuration,
                  			    minEpisodeDuration = minHuntEpisodeDuration,
                            groupHuntRatio=nHuntingDuration/ntotalFrames,
                            meanHuntRatioPerLarva=meanHuntRatio,
                            stdHuntRatioPerLarva=sdHuntRatio,
                            seHuntRatioPerLarva=sdHuntRatio/sqrt(nLarva),
                            medHuntRatioPerLarva=medHuntRatio,
                            maxHuntRatioPerLarva=maxHuntRatio,
                            minHuntRatioPerLarva=minHuntRatio,
                            groupHuntEvents=nHuntingEvents,
                            meanHuntingEventsPerLarva=meanHuntingEvents,
                            stdHuntingEventsPerLarva=sdDevHuntingEvents,
                            seHuntingEventsPerLarva=sdDevHuntingEvents/sqrt(nLarva),
                            seHuntingEventsPerLarva=sdDevHuntingEvents/sqrt(nLarva),
                            medHuntingEventsPerLarva=medHuntingEvents,
                            maxHuntingEventsPerLarva=maxHuntingEvents,
                            minHuntingEventsPerLarva=minHuntingEvents,
                  			    vHuntingEventsList= lHuntingEvents,
                  			    vdHuntingEventsList=datHuntingEvents
                            )
  
    message(paste("Number of Hunting Events detected:",nHuntingEvents, " Mean:", meanHuntingEvents, " SD:",sdDevHuntingEvents));
    message(paste("Hunting Episode Duration :",meanLarvaHuntDuration, " Episode Mean:", meanHuntEpisodeDuration, " SD:",sdHuntEpisodeDuration));
  
    message(paste("Logged ",NROW(datHuntingEvents)," hunt events") )
    
    return (lGroupHuntStats)
}


##Focus on extracting and Identifying Hunting events - Eye Vergence
##Return list of HuntEvents / With start and End Frame / and Video FileName
getHuntEvents <- function(datAllFrames,vlarvaID)
{
  message("## Extract Hunting Events For Group ##")
  lGroupHunting    <- list()
  lHuntingEvents    <- list()
  #lHuntingDuration <-list()
  #vHuntStartFrames <- list()  
  #vHuntEndFrames <- list()  
  
  try(rm("vHuntStartFrames","vHuntEndFrames","vHuntDeltaFrames","lHuntingDuration"),silent=TRUE)
  
  
  nHuntingEventsForLarva   <- 0
  meanHuntingEvents <- 0;
  stdDevHuntingEvents <- 0;
  idx = 0;
  idxHuntRec = 0;
  nLarva  	<- length(vlarvaID)
  
  
  for (i in vlarvaID)
  {
    idx = idx+1;
    datLarvaFrames <- datAllFrames[datAllFrames$larvaID == i,]
    
    #Used to Identify Experimet Set From Which Data COmes from - Plotted AS different colour
    DataSetID             <- ifelse(any(names(datLarvaFrames) == "dataSet"),unique(datLarvaFrames$dataSet),0 )
    stopifnot(DataSetID >= 0 )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID ) 
    
    vEventID = unique((datLarvaFrames$eventID))
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFrames[datLarvaFrames$REyeAngle > -35 & datLarvaFrames$REyeAngle <  35 & datLarvaFrames$LEyeAngle > -35 & datLarvaFrames$LEyeAngle <  35,]
    
    nTotalHuntFrames = 0;
    nHuntingEventsForLarva <- 0;
    meanHuntDuration <- sdHuntDuration  <- medHuntDuration <- maxHuntDuration<-minHuntDuration <- 0
    
    lHuntingDuration <-list() ##Empty List Of Event Durations

    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    for (k in vEventID)
    {
      
      ##Select Hunt Frames of this Event ## 
      ## Note Criterion - Left Eye Angle >  Min & Right Eye Angle <  Min (Turned Inwards) & Vergence Angle > 40
      ## Warning - This criterion Is repeated in the plot Function of plotTrackScatterAndDensities - Warning 
      datHuntFrames    <- datLarvaFrames[datLarvaFrames$eventID == k & 
                                           datLarvaFrames$REyeAngle < -G_THRESHUNTANGLE &
                                           datLarvaFrames$LEyeAngle > G_THRESHUNTANGLE & 
                                           abs(datLarvaFrames$LEyeAngle-datLarvaFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE,]
      nTotalHuntFrames <- nTotalHuntFrames + length(datHuntFrames$frameN)
      ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
      ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
      ## Detects End Events as +ve differences
      vHuntDeltaFrames <- diff(datHuntFrames$frameN,lag=1,differences=1)
      
      #Find Start Hunt Event the point leading a big gap in frameN, ie reverse diff s[(1+1):10] - s[(1):9]
      ##ShifHuntDelteFrames
      vsHuntDeltaFrames     <- (1:(NROW(datHuntFrames$frameN)+2))
      vHuntStartFrames      <- list()
      vHuntEndFrames        <- list()
      lHuntingDuration[[k]] <- 0
      minHuntDuration       <- 0
      
      
      # Debug #
      #message(paste(i,".",k,". Detected Hunt #Frames=", NROW(datHuntFrames$frameN) ))

      ##Consider Only if there are at least 100 frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
      if  (NROW(datHuntFrames$frameN) > G_MINEPISODEDURATION)
      {
        ##Add (Imaginary) Edge Frame  Numbers  (preceding hunting and after last hunting Frame)
        vsHuntDeltaFrames[1] <- datHuntFrames$frameN[1]-G_MINGAPBETWEENEPISODES
        vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- datHuntFrames$frameN[NROW(datHuntFrames$frameN)]+G_MINGAPBETWEENEPISODES
        ##Copy into Shifted right (lagged) Position
        vsHuntDeltaFrames[2:(NROW(datHuntFrames$frameN)+1)] <- datHuntFrames$frameN
        ##Do Rev Diff - ie Nabla taking s_n = x_n-x_{n-1} # Detect Start Events as +ve diff
        vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] <- datHuntFrames$frameN - vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN))]
        ##\TODO:  Perhaps Throw away Hunts That Began Before Recording Started or End After It Ends - remove Edge Points ##
        
        ##Find Hunt Event STARTs ## Needs A gap of At least #MINGAP Between Events To Be a new Episode
        ##attempting to Remove delay introduced by filtering 
        vHuntStartFrames <- datHuntFrames$frameN[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
        ##Find Hunt Event ENDs - Needs To Be atleast G_MINGAPBETWEENEPISODES Away from next start
        ##Shift (expanded) DFrame vector to pickup on adjacent Hunt-Frame number signifying the end of the Hunt Episode
        vsHuntDeltaFrames[2:NROW(vsHuntDeltaFrames)] <- vsHuntDeltaFrames[1:NROW(vsHuntDeltaFrames)-1]
        ##Pickout End frames by selected frameN locations shifted by 1+
        vHuntEndFrames  <- datHuntFrames$frameN[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
        ###
        
        #message(paste("Found i=",length(vHuntEndFrames)," Hunt End points and k=",length(vHuntStartFrames), "Start Points"))
        stopifnot(length(vHuntEndFrames) == length(vHuntStartFrames))
        lHuntingDuration[[k]]   <- vHuntEndFrames - vHuntStartFrames
        message(paste(" Event Hunt Duration is  ",sum(lHuntingDuration[[k]])," in ",length(vHuntStartFrames)," episodes for larvaID:",i," eventID:",k))
        
        # idxHuntRec = idxHuntRec + 1;
        # lHuntingEvents[[idxHuntRec]] <- list(larvaID =rep(i,length(vHuntEndFrames)),
        #                                      eventID = rep(k,length(vHuntEndFrames)),
        #                                      dataSetID = rep(DataSetID,length(vHuntEndFrames)),
        #                                      fileIdx = rep(unique(datHuntFrames$fileIdx),length(vHuntEndFrames)),
        #                                      startFrame=unlist(vHuntStartFrames),
        #                                      endFrame=unlist(vHuntEndFrames),
        #                                      huntScore=rep(0,length(vHuntEndFrames))) ##0 To Be used for Manual Labelling Of Hunting Success
        # 
        idxHuntRec = idxHuntRec + 1;
        lHuntingEvents[[idxHuntRec]] <- data.frame(larvaID =i,
                                                   eventID = k,
                                                   dataSetID = DataSetID,
                                                   fileIdx = unique(datHuntFrames$fileIdx),
                                                   startFrame=unlist(vHuntStartFrames),
                                                   endFrame=unlist(vHuntEndFrames),
                                                   huntScore=0) ##0 To Be used for Manual Labelling Of Hunting Success
        
        if( lHuntingDuration[[k]] < 0) 
          stop(paste("Negative Hunt Duration detected LarvaID:",i," eventID:",k," Possible Duplicate tracker file" ) ) ##Catch Error / Can be caused by Duplicate Tracked Video - Check )
        
        minHuntDuration  <- min(lHuntingDuration[[k]])
        
        nHuntingEventsForLarva  <- nHuntingEventsForLarva + length(vHuntStartFrames) ##increment Event Count
      }else{
        ##Hunting Episode Insufficient ##	
        # vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN)+1)] <- 0  ##Copy into Shifted Position
        warning(paste("Ignoring low Hunt Frame Number",NROW(datHuntFrames$frameN)," of Larva:",i,"event",k) )
        if (is.na(k) | k < 1 )
        {
          stop(paste("Invalid Event Idx:",k) )
        }
        minHuntDuration       <-0
        lHuntingDuration[[k]] <- 0 ##Otherwise Unused Event IDs end up NULL
      }
      
      
      #### ERROR Handlers /  DEbug ###
      
      if (is.na(minHuntDuration)) 
      {  
        message(paste("Warning - Min Episode Duration is  ",minHuntDuration," for larvaID:",i," eventID:",k))
        warning(paste("Warning - Min Episode Duration is  ",minHuntDuration," for larvaID:",i," eventID:",k))
        # message("** Probably dublicate track files present **")
      }else if (minHuntDuration < 0) 
      {  
        
        message(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
        warning(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
        message("** Probably dublicate track files present **")
      }
      
      if( any( is.null( lHuntingDuration ) ) )  
        stop(paste("Null Hunting duration for larvaID:",i," eventID:",k,". Possible Missing Event IDs",  as.character( unlist(which(sapply(lHuntingDuration, is.null) == TRUE)) )) )
      
    }##For each Event Of this Larva  ##########  
    #####

    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
  } ### For Each Larva In Group ###
  ##### Note That vHuntStartFrames is invalid beyond this point #
  
  datHuntingEvents = do.call(rbind,lHuntingEvents)
  #  datGroupHunting = as.data.frame(lGroupHunting)

  nHuntingEvents            <- sum(vnHuntingEvents);
  message("Hunting Event For This Group:",nHuntingEvents)
  

  if (nHuntingEvents < 1)
  {
    warning("No Hunting Events Detected in Group!")
    message("-* No Hunting Events Detected in Group! *- ")
    #
    #meanHuntRatio       = 0
    #sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
  }
  
  message(paste("Logged ",NROW(datHuntingEvents)," hunt events") )
  
  return (datHuntingEvents)
}

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
