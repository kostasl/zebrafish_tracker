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

calcHuntStat <- function(datAllFrames,vlarvaID)
{
  lGroupHunting= list()
  lHuntingDuration = list()
  nHuntingEvents=0;
  meanHuntingEvents = 0;
  stdDevHuntingEvents = 0;
  idx = 0;
  nLarva  	<- length(vlarvaID)
  
  for (i in vlarvaID)
  {
    idx = idx+1;
    datLarvaFrames <- datAllFrames[datAllFrames$larvaID == i,]
    vEventID = unique(datLarvaFrames$eventID)
    nTotalHuntFrames = 0;
    
    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    for (k in vEventID)
    {

	    ##Select Hunt Frames of this Event ##
	    datHuntFrames <- datLarvaFrames[datLarvaFrames$eventID == k & datLarvaFrames$REyeAngle < -20 & datLarvaFrames$LEyeAngle > 20,]
            nTotalHuntFrames = nTotalHuntFrames + length(datHuntFrames$frameN)
	    ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
	    ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
	    ## Detects End Events as +ve differences
	    vHuntDeltaFrames <- diff(datHuntFrames$frameN,lag=1,differences=1)
		
	    #Find Start Hunt Event the point leading a big gap in frameN, ie reverse diff s[(1+1):10] - s[(1):9]
	    ##ShifHuntDelteFrames
	    vsHuntDeltaFrames <- (1:(NROW(datHuntFrames$frameN)+1))
	    if  (NROW(datHuntFrames$frameN)>100)
	    {
		   #Make Sure We count Single Hunt Events occuring for a larva -
		   vHuntDeltaFrames[NROW(vHuntDeltaFrames)] = 400 ##Mark End Of Last Hunting Episode
		   vsHuntDeltaFrames[1] = 400 #Mark Hunt Start as First Diff Point (400 ie above Threshold)
		    ##Copy into Shifted right (lagged) Position
	 	   vsHuntDeltaFrames[2:(NROW(datHuntFrames$frameN)+1)] <- datHuntFrames$frameN
		    ##Do Rev Diff - ie Nabla taking s_n = x_n-x_{n-1} # Detect Start Events as +ve diff
		    vsHuntDeltaFrames <- datHuntFrames$frameN - vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN))]

		    
		    ##Find Hunt Event STARTs
		    vHuntStartFrames <-  datHuntFrames$frameN[vsHuntDeltaFrames> 300 ] ##Note Ignores hunting event at the very beginning of Track (<300 frames away from start)
		    ##Find Hunt Event ENDs
		    vHuntEndFrames <- datHuntFrames$frameN[vHuntDeltaFrames > 300 ] ##Note Ignores hunting event at the very beginning of Track (<300 frames away from start)

		    #message(paste("Found i=",length(vHuntEndFrames)," Hunt End points and k=",length(vHuntStartFrames), "Start Points"))
		    
		    lHuntingDuration[[k]] <- vHuntEndFrames - vHuntStartFrames
	    }else{
                 ##Hunting Episode Insufficient ##	
#	 	   vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN)+1)] <- 0  ##Copy into Shifted Position
	            lHuntingDuration[[k]] <- 0
	    }
     }##For each Event Of this Larva    
   
    vHuntingDuration = unlist(lHuntingDuration)
    if (length(vHuntingDuration) > 0)
    {
	    meanHuntDuration <- mean(vHuntingDuration)
	    sdHuntDuration <-   sd(vHuntingDuration)
	    medHuntDuration  <- median(vHuntingDuration)
	    maxHuntDuration  <- max(vHuntingDuration)
	    minHuntDuration  <- min(vHuntingDuration)
	##Debug
	if (minHuntDuration < 0) 
	{
		message(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
		warning(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for larvaID:",i," eventID:",k))
	}
    }else {
	medHuntDuration = sdHuntDuration = maxHuntDuration = minHuntDuration = 0
	}
    
    ##Total Duration of Hunting 
    ##+ Number of Hunting Events Counted as points where eyevergence events are more than 300 frames apart  ##
    lGroupHunting[[idx]] = list(larvaID=i,
                                count=length(vHuntStartFrames),
                                huntframes=nTotalHuntFrames,
				vHuntDurations=vHuntingDuration,      	                        
				meanHDuration=meanHuntDuration,
				sdHDuration=sdHuntDuration,
				seHDuration=sdHuntDuration/sqrt(nLarva),
                                medHDuration=medHuntDuration,
			        maxHDuration = maxHuntDuration,
			        minHDuration = minHuntDuration,
                                totalframes=length(datLarvaFrames$frameN))
  } #For Each Larva In Group

  datGroupHunting = do.call(rbind,lGroupHunting)
#  datGroupHunting = as.data.frame(lGroupHunting)
  vHuntEpisodes = as.vector(do.call(rbind,lapply(lGroupHunting,"[[","vHuntDurations")))

  vHuntRatio         <- unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"]) 
  vHuntingDuration   <- unlist(datGroupHunting[,"huntframes"])
  vHuntingEvents   <- unlist(datGroupHunting[,"count"])
  nHuntingDuration = sum(vHuntingDuration)
  nHuntingEvents = sum(vHuntingEvents);
  ntotalFrames  <- sum(unlist(datGroupHunting[,"totalframes"]))



  if (nHuntingEvents > 1)
  {
    meanHuntRatio      <- mean(vHuntRatio)
    sdHuntRatio        <- sd(vHuntRatio) 
    
    medHuntRatio       <- median(vHuntRatio)
    maxHuntRatio    <- max(vHuntRatio)
    minHuntRatio    <- min(vHuntRatio)

    ### Hunt Duration Per Larva ###    
    meanLarvaHuntDuration <- mean(vHuntingDuration)
    sdLarvaHuntDuration   <- sd(vHuntingDuration)
    medLarvaHuntDuration  <- median(vHuntingDuration)
    maxLarvaHuntDuration  <- max(vHuntingDuration)
    minLarvaHuntDuration  <- min(vHuntingDuration)
	
    ### Hunt Duration Per Hunting Episode ###
    meanHuntEpisodeDuration <- mean(vHuntEpisodes)
    sdHuntEpisodeDuration   <- sd(vHuntEpisodes)
    medHuntEpisodeDuration  <- median(vHuntEpisodes)
    minHuntEpisodeDuration  <- min(vHuntEpisodes)
    maxHuntEpisodeDuration  <- max(vHuntEpisodes)



    ### Number of Hunting Events Per Larva ##
    meanHuntingEvents  <- mean(vHuntingEvents);
    sdDevHuntingEvents <- sd(vHuntingEvents);
    medHuntingEvents   <- median(vHuntingEvents);
    maxHuntingEvents   <- max(vHuntingEvents);
    minHuntingEvents   <- min(vHuntingEvents);
  }
  else
  {
    warning("No Hunting Events Detected!")
    meanHuntRatio       = NA
    sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
    
    meanHuntingEvents = 0;
    sdDevHuntingEvents = 0;
    
  }
    ###Return Results As List of Variable
    lGroupHuntStats <- list(nLarva=nLarva,
                            totalFrames=ntotalFrames,
                            huntFrames=nHuntingDuration,
			    meanDuration=meanLarvaHuntDuration,
			    sdDuration=sdLarvaHuntDuration,
			    seDuration=sdHuntDuration/sqrt(nLarva),
			    medDuration=medLarvaHuntDuration,
			    maxDuration = maxLarvaHuntDuration,
			    minDuration = minLarvaHuntDuration,
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
                            minHuntingEventsPerLarva=minHuntingEvents
                            )
  
    message(paste("Number of Hunting Events detected:",nHuntingEvents, " Mean:", meanHuntingEvents, " SD:",sdDevHuntingEvents));
    message(paste("Hunting Episode Duration :",meanLarvaHuntDuration, " Episode Mean:", meanHuntEpisodeDuration, " SD:",sdHuntEpisodeDuration));

    return (lGroupHuntStats)
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
