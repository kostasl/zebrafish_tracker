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
    #message(f)
    strfilepath <- list.files(path =strVideoFilePath, pattern = f, all.files = FALSE,
                              full.names = TRUE, recursive = TRUE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    #message(strfilepath)
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


loadHuntEvents <- function(strCondTags,dataSetsToProcess)
{
  
  lHuntStat <- list();
  for (i in strCondTags)
  {
    message(paste("#### ProcessGroup ",i," ###############"))
    strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
    
    load(file=paste(strDataFileName,".RData",sep="" )) ##Loads datHuntEvent
    lHuntStat[[i]] <- calcHuntStat2(datHuntEvent)
    
  }
  
  datHuntStat = do.call(rbind,lHuntStat)
  
  return(datHuntStat)
  
}


######################################################### EXTRACT HUNTING EVENTS ##################################
##Focus on extracting and Identifying Hunting events - Eye Vergence
##Return list of HuntEvents / With start and End Frame / and Video FileName
detectHuntEvents <- function(datAllFrames,vexpID,vdatasetID)
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
  nLarva  	<- length(vexpID)
  
  
  for (i in vexpID)
  {
    idx = idx+1;
    datLarvaFramesRaw <- datAllFrames[which(datAllFrames$expID == i),]
    
    #Used to Identify Experimet Set From Which Data COmes from - Plotted AS different colour
    DataSetID             <- ifelse(any(names(datLarvaFramesRaw) == "dataSet"),unique(datLarvaFramesRaw$dataSet),0 )
    stopifnot(DataSetID >= 0 )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID ) 
    
    vEventID = unique((datLarvaFramesRaw$eventID))
    groupID <- unique((datLarvaFramesRaw$group))
    larvaID = unique((datLarvaFramesRaw$larvaID))
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFramesRaw[which(datLarvaFramesRaw$REyeAngle > -35 & datLarvaFramesRaw$REyeAngle <  35 &
                                                datLarvaFramesRaw$LEyeAngle > -35 & datLarvaFramesRaw$LEyeAngle <  35),]
    
    nTotalHuntFrames      <- 0
    nTotalRecordedFrames  <- NROW(datLarvaFrames) ##Will be calculated as diff in frame N for each Event (Each Event FrameN starts from 0 in tracker)
    
    nHuntingEventsForLarva <- 0;
    meanHuntDuration <- sdHuntDuration  <- medHuntDuration <- maxHuntDuration<-minHuntDuration <- 0
    
    lHuntingDuration <-list() ##Empty List Of Event Durations

    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    nInitialPrey = NA ##Count THe Total Prey Count On 1st Event
    for (k in vEventID)
    {
      
        ##Select Hunt Frames of this Event ## 
        ## Note Criterion - Left Eye Angle >  Min & Right Eye Angle <  Min (Turned Inwards) & Vergence Angle > 40
        ## Warning - This criterion Is repeated in the plot Function of plotTrackScatterAndDensities - Warning
        ## Warning Without the use of which - NA entries Appear!
       ## Take Raw Frames Of Event So we can obtain the Number of Rotifers
        datEventFrames   <- datLarvaFramesRaw[which(datLarvaFramesRaw$eventID == k),]
        
        
        datHuntFrames    <- datLarvaFrames[which(datLarvaFrames$eventID == k &
                                          (datLarvaFrames$REyeAngle < -G_THRESHUNTANGLE |
                                           datLarvaFrames$LEyeAngle > G_THRESHUNTANGLE) &
                                        abs(datLarvaFrames$LEyeAngle-datLarvaFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE),]
        
        nTotalHuntFrames <- nTotalHuntFrames + NROW(datHuntFrames$frameN)
        
        ##On 1st Event Save THe Initial mean Prey Count
        if (k==1 && NROW(datEventFrames) > 0)
          nInitialPrey <- mean(datEventFrames$PreyCount)
        
        ##Add to Total Number of frames - Using Actual Frame Number of steps Instead of simpler NRow(datLarvaFrames) 
        ##(ideally these should match tho) but tracker can loose object for a few frames and skip
        #nTotalRecordedFrames <- nTotalRecordedFrames + sum(diff(datLarvaFrames[datLarvaFrames$eventID == k,] ))
        
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
  
        ##Consider Only if there are at a min of Hunting frames in this Event  - Surragate to proper Hunt Length Filter frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
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
          #Add ENd Edge of DS Vector 
          vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- G_MINGAPBETWEENEPISODES
          ##Pickout End frames by selected frameN locations shifted by 1+ 
          ##End Frames Selected By FrameDIff > Than MinGapBetween Episodes
          ## But it can happen that a single episode is shorter than G_MINGAPBETWEENEPISODES if G_MINEPISODEDURATION < G_MINGAPBETWEENEPISODES
          vHuntEndFrames  <- datHuntFrames$frameN[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
          ###
          
          #message(paste("Found i=",length(vHuntEndFrames)," Hunt End points and k=",length(vHuntStartFrames), "Start Points"))
          if (length(vHuntEndFrames) != length(vHuntStartFrames))
          {
            stop(paste("length(vHuntEndFrames)",length(vHuntEndFrames), "!= length(vHuntStartFrames)",length(vHuntStartFrames)," - check CSVs for expID:",i," event",k," DatasetId:",DataSetID, " fileID:",unlist(unique(datHuntFrames$fileIdx) ) )  )
          }
          if (length(unique(datHuntFrames$fileIdx) ) > 1)
          {
          
            stop(paste("Duplicate Track File for same Event - check CSVs for expID:",i," event",k," DatasetId:",DataSetID,  " fileIds:",unlist(unique(datHuntFrames$fileIdx) ) )  )
    
          }
          
          lHuntingDuration[[k]]   <- vHuntEndFrames - vHuntStartFrames
          message(paste(" Event Hunt Duration is  ",sum(lHuntingDuration[[k]])," in ",length(vHuntStartFrames)," episodes for expID:",i," eventID:",k))
          
  
          idxHuntRec = idxHuntRec + 1;
          
          lHuntingEvents[[idxHuntRec]] <- data.frame(expID      = factor(i,levels=vexpID),
                                                     eventID    = k,
                                                     dataSetID  = factor(DataSetID,levels=vdatasetID),
                                                     larvaID    = factor(larvaID,levels=seq(1:4)) , ##Identifies Larva Between Empty And Live Test Conditions
                                                     groupID    = groupID,
                                                     fileIdx    = unique(datHuntFrames$fileIdx),
                                                     startFrame = unlist(vHuntStartFrames),
                                                     endFrame   = unlist(vHuntEndFrames),
                                                     nExpFrames = nTotalRecordedFrames,
                                                     InitPreyCount =  nInitialPrey, ##Mean Prey Count Across Hunt Frames
                                                     PreyCount =  mean(datHuntFrames$PreyCount), ##Mean Prey Count Across Hunt Frames
                                                     huntScore  = 0,
                                                     stringsAsFactors = FALSE) ##0 To Be used for Manual Labelling Of Hunting Success
          
          if( lHuntingDuration[[k]] < 0) 
            stop(paste("Negative Hunt Duration detected expID:",i," eventID:",k," Possible Duplicate tracker file" ) ) ##Catch Error / Can be caused by Duplicate Tracked Video - Check )
          
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
          message(paste("Warning - Min Episode Duration is  ",minHuntDuration," for expID:",i," eventID:",k))
          warning(paste("Warning - Min Episode Duration is  ",minHuntDuration," for expID:",i," eventID:",k))
          # message("** Probably dublicate track files present **")
        }else if (minHuntDuration < 0) 
        {  
          
          message(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for expID:",i," eventID:",k))
          warning(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for expID:",i," eventID:",k))
          message("** Probably dublicate track files present **")
        }
        
        if( any( is.null( lHuntingDuration ) ) )  
          stop(paste("Null Hunting duration for expID:",i," eventID:",k,". Possible Missing Event IDs",  as.character( unlist(which(sapply(lHuntingDuration, is.null) == TRUE)) )) )
        
        
      
    }##For each Event Of this Larva  ##########  
    #####

    ##For Larva That Did not register any sufficient Hunting Events - Add An Empty Record To Acknowledge 
    if (nHuntingEventsForLarva == 0)
    {
      nmeanPreyCount = NA
      idxHuntRec = idxHuntRec + 1;
      
      if (NROW(datEventFrames) > 2)
        nmeanPreyCount <- mean(datEventFrames$PreyCount)
      
      if (NROW(datEventFrames) == 1)
        nmeanPreyCount = datEventFrames$PreyCount

      lHuntingEvents[[idxHuntRec]] <- data.frame(expID      = factor(i,levels=vexpID),
                                                 eventID    = 0,
                                                 dataSetID  = factor(DataSetID,levels=vdatasetID),
                                                 larvaID    = factor(larvaID,levels=seq(1:4)),
                                                 groupID    = groupID,
                                                 fileIdx    = 0,
                                                 startFrame = 0,
                                                 endFrame   = 0,
                                                 nExpFrames = nTotalRecordedFrames,
                                                 InitPreyCount = nInitialPrey,
                                                 PreyCount  = nmeanPreyCount,
                                                 huntScore  = 0,
                                                 stringsAsFactors = FALSE) ##0 To Be used for Manual Labelling Of Hunting Success
      
    }
    
    
    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
    
  } ### For Each Larva In Group ###
  ##### Note That vHuntStartFrames is invalid beyond this point #
  
  datHuntingEvents = do.call(rbind,lHuntingEvents)
  
  
  #stopifnot(any(is.nan(datHuntingEvents) ) ==FALSE )
  #  datGroupHunting = as.data.frame(lGroupHunting)

  nHuntingEvents            <- length(datHuntingEvents);
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

####### Calc Hunt Statistics Giving list of Hunting Events of a group - 
## Note: 
calcHuntStat2 <- function(datHuntEvent)
{

  message(paste("##Calculate Hunting Statitistics for Group ",unique(datHuntEvent$groupID), " ##" ) )
  
  stopifnot(length(unique(datHuntEvent$groupID) ) ==1 ) ##Only One Condition Should Be analysed at a time - OtherWise LarvaID may mix results between conditions 
  
  ##This Method Produces The Vector WIth The zero Values for Non Hunting Larvae - but mean And SD are correct with aggregate Method
  tblHuntDurationPerLarva  <- tapply(datHuntEvent$endFrame-datHuntEvent$startFrame, datHuntEvent$expID,sum)
  tblHuntDurationPerLarva  <- replace(tblHuntDurationPerLarva,is.na(tblHuntDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  
  ###But  Mean Hunting Duration Per Larva is also correct if using aggregate  ##
  #tblHuntDurationsPerLarva <- aggregate(x=datHuntEvent$endFrame-datHuntEvent$startFrame, by=list(expID=datHuntEvent$expID), FUN=sum)
  
  ## Hunting Duration Per Episode ##
  
  ##Replace 0s with NA - Ignore 0 Duration Episodes
  datHuntEventNonZeroEpi <- datHuntEvent[datHuntEvent$endFrame-datHuntEvent$startFrame > 0,] 
  ##Calc number of hunting events stat per Experiment
  tblHuntsCounts<-table(datHuntEventNonZeroEpi$expID)
  ##Given An Episode Occured
  tblMeanEpisodeDurationPerLarva  <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,mean)
  #tblMeanEpisodeDurationPerLarva  <- replace(tblMeanEpisodeDurationPerLarva,is.na(tblMeanEpisodeDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  tblMedianEpisodeDurationPerLarva  <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,median)
  ##Min Should Return Smallest Episode Length - ie collects from where episodes did happen to get the shortest one
  tblMinEpisodeDurationPerLarva   <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,min,na.rm=TRUE)
  tblMaxEpisodeDurationPerLarva   <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,max,na.rm=TRUE)
  #mean(datHuntEvent$endFrame-datHuntEvent$startFrame)
  #tblMeanEpisodeDurationPerLarva  <- tapply(datHuntEvent$endFrame-datHuntEvent$startFrame, datHuntEvent$expID,mean)
  #tblMeanEpisodeDurationPerLarva  <- replace(tblMeanEpisodeDurationPerLarva,is.na(tblMeanEpisodeDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  ### Hunt Ratio /Need to total Frames For That - if no hunt record then also no TotalRec Frames here!
  ##Total Recorded Duration Per Larva - This value is repeated at each hunt event of a larva and so taking the mean would 
  # give the individual value - Missing ExpID mean no Data 
  # Empty Record Durations 0 are given for Larva Without Hunting Episodes, So as to obtain the total number of recorded frames
  # ExpID with No Data Will Appear As NA
  tblRecDurationPerLarva <- tapply(datHuntEvent$nExpFrames, datHuntEvent$expID,mean, na.rm=TRUE) 
  tblRecDurationPerLarva <- replace(tblRecDurationPerLarva,tblRecDurationPerLarva==0,1) #Set 0s Ie Missing/Empty Tracks to 1 so to include in Ration Calculation as missing Larvae
  
  ## PREY COUNTS is averaged Per Experiment###
  tblAvailablePreyCountPerLarva <-tapply(datHuntEvent$PreyCount, datHuntEvent$expID,mean, na.rm=TRUE)  
  tblAvailablePreyCountPerLarva <- replace(tblAvailablePreyCountPerLarva,is.nan(tblAvailablePreyCountPerLarva),NA)
  
  tblAvailableInitialPreyCountPerLarva <-tapply(datHuntEvent$InitPreyCount, datHuntEvent$expID,mean, na.rm=TRUE)  
  tblAvailableInitialPreyCountPerLarva <- replace(tblAvailableInitialPreyCountPerLarva,is.nan(tblAvailableInitialPreyCountPerLarva),NA)
  
  stopifnot( NROW(tblAvailablePreyCountPerLarva[ is.nan(tblAvailablePreyCountPerLarva) ]) == 0 )
  ##Calc Sum Of Hunt Ratios of each Episode, to obtain Hunt Ratio of experiment
  tblHuntRatioPerLarva   <- tapply( (datHuntEvent$endFrame-datHuntEvent$startFrame)/datHuntEvent$nExpFrames, datHuntEvent$expID,sum,na.rm=TRUE) 
  
  
  
  ##Number of Appearances in Hunt Events - (Minimum Is 1 - so as to obtain LarvaID-ExpID Link)  /  ##
  tblHunts <-table(datHuntEvent$larvaID,datHuntEvent$expID,datHuntEvent$dataSetID)
  
  ##Exp ID LookUp Table So Ican locate the same larva Across Empty->Live Condition
  ## 
  datLT <- data.frame(cbind(larvaID=levels(datHuntEvent$larvaID)[datHuntEvent$larvaID],expID=levels(datHuntEvent$expID)[datHuntEvent$expID],dataSetID=levels(datHuntEvent$dataSetID)[datHuntEvent$dataSetID]))
  udatLT <- unique(datLT)
  
  
  ## COmbine Results Into One Handy Structure containing precalculated results ##
  
  nLarva <- nrow(tblHuntDurationPerLarva)
  nAvailInitialPreyCount <- nrow(tblAvailableInitialPreyCountPerLarva[is.na(tblAvailableInitialPreyCountPerLarva) == FALSE])
  nPreyCountSamples <- nrow(tblAvailablePreyCountPerLarva[is.na(tblAvailablePreyCountPerLarva) == FALSE]) 
  lGroupHuntStats <- list(nLarva                        = nLarva,
                          vHPreyCount                   = tblAvailablePreyCountPerLarva,
                          vHInitialPreyCount            = tblAvailableInitialPreyCountPerLarva,
                          nPreyCountSamples             = nPreyCountSamples, ##COunt of Samples Used to Calc Mean
                          meanPreyCount                 = mean(tblAvailablePreyCountPerLarva,na.rm=TRUE),
                          sdPreyCount                   = sd(tblAvailablePreyCountPerLarva,na.rm=TRUE),
                          sePreyCount                   = sd(tblAvailablePreyCountPerLarva,na.rm=TRUE)/sqrt(nPreyCountSamples),
                          minPreyCount                  = min(tblAvailablePreyCountPerLarva,na.rm=TRUE),
                          medianPreyCount               = median(tblAvailablePreyCountPerLarva,na.rm=TRUE),
                          maxPreyCount                  = max(tblAvailablePreyCountPerLarva,na.rm=TRUE),
                          ninitPreyCountSamples         = nAvailInitialPreyCount,
                          initPreyCount                 = mean(tblAvailableInitialPreyCountPerLarva,na.rm=TRUE), ##Mean Count At Event 1
                          initsePreyCount               = sd(tblAvailableInitialPreyCountPerLarva,na.rm=TRUE)/sqrt(nAvailInitialPreyCount), ##Mean Count At Event 1
                          totalFrames                   = sum(tblRecDurationPerLarva),
                          vDataSetID                    = levels(datHuntEvent$dataSetID),
                          vExpID                        = levels(datHuntEvent$expID),
                          vLarvaID                      = levels(datHuntEvent$larvaID),
                          vIDLookupTable                = udatLT,
                          vHDurationPerLarva            = tblHuntDurationPerLarva,
                          vmeanHEpisodeDurationPerLarva	= tblMeanEpisodeDurationPerLarva,
                          vHLarvaEventCount             = tblHuntsCounts,
                          groupHuntEvents               = sum(tblHuntsCounts),
                          meanHuntingEventsPerLarva     = mean(tblHuntsCounts),
                          stdHuntingEventsPerLarva      = sd(tblHuntsCounts),
                          seHuntingEventsPerLarva       = sd(tblHuntsCounts)/sqrt(nLarva),
                          medHuntingEventsPerLarva      = median(tblHuntsCounts),
                          maxHuntingEventsPerLarva      = max(tblHuntsCounts),
                          minHuntingEventsPerLarva      = min(tblHuntsCounts),
                          totalHuntFrames               = sum(tblHuntDurationPerLarva),
                          vLarvaHRatio                  = tblHuntRatioPerLarva,
                          meanDuration                  = mean(tblHuntDurationPerLarva),
                          sdDuration                    = sd(tblHuntDurationPerLarva),
                          seDuration                    = sd(tblHuntDurationPerLarva)/sqrt(nLarva),
                          medDuration                   = median(tblHuntDurationPerLarva),
                          meanEpisodeDuration           = mean(tblMeanEpisodeDurationPerLarva,na.rm=TRUE),
                          sdEpisodeDuration             = sd(tblMeanEpisodeDurationPerLarva,na.rm=TRUE), ##Std Dev Of Means
                          seEpisodeDuration             = sd(tblMeanEpisodeDurationPerLarva,na.rm=TRUE)/sqrt(nLarva), ##Std Err Of Means Of Episode Means
                          medEpisodeDuration            = median(tblMedianEpisodeDurationPerLarva,na.rm=TRUE), #Median OF Medians
                          maxEpisodeDuration            = max(tblMaxEpisodeDurationPerLarva,na.rm=TRUE),
                          minEpisodeDuration            = min(tblMinEpisodeDurationPerLarva,na.rm=TRUE),
                          groupHuntRatio                = sum(tblHuntDurationPerLarva)/sum(tblRecDurationPerLarva), ###Need Actual Frame Numbers Here
                          meanHuntRatioOfGroup          = mean(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          stdHuntRatioOfGroup           = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          seHuntRatioOfGroup            = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva)/sqrt(nLarva),
                          medHuntRatioPerLarva          = median(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          maxHuntRatioPerLarva          = max(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          minHuntRatioPerLarva          = min(tblHuntDurationPerLarva/tblRecDurationPerLarva)
  )
  
  
} ### END oF CALC Stat ##



####### Calc Hunt Statistics Giving list of Hunting Events of a group - 
## Note:  Aadded Nabla Of Prey Count
calcHuntStat3 <- function(datHuntEvent)
{
  
  message(paste("##V3 Calculate Hunting Statitistics for Group ",unique(datHuntEvent$groupID), " ##" ) )
  
  stopifnot(length(unique(datHuntEvent$groupID) ) ==1 ) ##Only One Condition Should Be analysed at a time - OtherWise LarvaID may mix results between conditions 
  
  ##This Method Produces The Vector WIth The zero Values for Non Hunting Larvae - but mean And SD are correct with aggregate Method
  tblHuntDurationPerLarva  <- tapply(datHuntEvent$endFrame-datHuntEvent$startFrame, datHuntEvent$expID,sum)
  tblHuntDurationPerLarva  <- replace(tblHuntDurationPerLarva,is.na(tblHuntDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  
  ###But  Mean Hunting Duration Per Larva is also correct if using aggregate  ##
  #tblHuntDurationsPerLarva <- aggregate(x=datHuntEvent$endFrame-datHuntEvent$startFrame, by=list(expID=datHuntEvent$expID), FUN=sum)
  
  ## Hunting Duration Per Episode ##
  
  ##Replace 0s with NA - Ignore 0 Duration Episodes
  datHuntEventNonZeroEpi <- datHuntEvent[datHuntEvent$endFrame-datHuntEvent$startFrame > 0,] 
  ##Calc number of hunting events stat per Experiment
  tblHuntsCounts<-table(datHuntEventNonZeroEpi$expID)
  ##Given An Episode Occured
  tblMeanEpisodeDurationPerLarva  <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,mean)
  #tblMeanEpisodeDurationPerLarva  <- replace(tblMeanEpisodeDurationPerLarva,is.na(tblMeanEpisodeDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  tblMedianEpisodeDurationPerLarva  <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,median)
  ##Min Should Return Smallest Episode Length - ie collects from where episodes did happen to get the shortest one
  tblMinEpisodeDurationPerLarva   <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,min,na.rm=TRUE)
  tblMaxEpisodeDurationPerLarva   <- tapply(datHuntEventNonZeroEpi$endFrame-datHuntEventNonZeroEpi$startFrame, datHuntEventNonZeroEpi$expID,max,na.rm=TRUE)
  #mean(datHuntEvent$endFrame-datHuntEvent$startFrame)
  #tblMeanEpisodeDurationPerLarva  <- tapply(datHuntEvent$endFrame-datHuntEvent$startFrame, datHuntEvent$expID,mean)
  #tblMeanEpisodeDurationPerLarva  <- replace(tblMeanEpisodeDurationPerLarva,is.na(tblMeanEpisodeDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  ### Hunt Ratio /Need to total Frames For That - if no hunt record then also no TotalRec Frames here!
  ##Total Recorded Duration Per Larva - This value is repeated at each hunt event of a larva and so taking the mean would 
  # give the individual value - Missing ExpID mean no Data 
  # Empty Record Durations 0 are given for Larva Without Hunting Episodes, So as to obtain the total number of recorded frames
  # ExpID with No Data Will Appear As NA
  tblRecDurationPerLarva <- tapply(datHuntEvent$nExpFrames, datHuntEvent$expID,mean, na.rm=TRUE) 
  tblRecDurationPerLarva <- replace(tblRecDurationPerLarva,tblRecDurationPerLarva==0,1) #Set 0s Ie Missing/Empty Tracks to 1 so to include in Ration Calculation as missing Larvae
  
  ## PREY COUNTS is averaged Per Experiment###
  tblPreyCountReductionPerLarvaHunt <-tapply(datHuntEvent$InitPreyCount-datHuntEvent$PreyCount, datHuntEvent$expID,mean, na.rm=TRUE)  
  tblPreyCountReductionPerLarvaHunt <- replace(tblPreyCountReductionPerLarvaHunt,is.nan(tblPreyCountReductionPerLarvaHunt),NA)
  
  tblAvailableInitialPreyCountPerLarva <-tapply(datHuntEvent$InitPreyCount, datHuntEvent$expID,mean, na.rm=TRUE)  
  tblAvailableInitialPreyCountPerLarva <- replace(tblAvailableInitialPreyCountPerLarva,is.nan(tblAvailableInitialPreyCountPerLarva),NA)
  
  stopifnot( NROW(tblPreyCountReductionPerLarvaHunt[ is.nan(tblPreyCountReductionPerLarvaHunt) ]) == 0 )
  ##Calc Sum Of Hunt Ratios of each Episode, to obtain Hunt Ratio of experiment
  tblHuntRatioPerLarva   <- tapply( (datHuntEvent$endFrame-datHuntEvent$startFrame)/datHuntEvent$nExpFrames, datHuntEvent$expID,sum,na.rm=TRUE) 
  
  
  ##Number of Appearances in Hunt Events - (Minimum Is 1 - so as to obtain LarvaID-ExpID Link)  /  ##
  tblHunts <-table(datHuntEvent$larvaID,datHuntEvent$expID,datHuntEvent$dataSetID)
  
  ##Exp ID LookUp Table So Ican locate the same larva Across Empty->Live Condition
  ## 
  datLT <- data.frame(cbind(larvaID=levels(datHuntEvent$larvaID)[datHuntEvent$larvaID],expID=levels(datHuntEvent$expID)[datHuntEvent$expID],dataSetID=levels(datHuntEvent$dataSetID)[datHuntEvent$dataSetID]))
  udatLT <- unique(datLT)
  
  
  ## COmbine Results Into One Handy Structure containing precalculated results ##
  
  nLarva <- nrow(tblHuntDurationPerLarva)
  nAvailInitialPreyCount <- nrow(tblAvailableInitialPreyCountPerLarva[is.na(tblAvailableInitialPreyCountPerLarva) == FALSE])
  nPreyCountSamples <- nrow(tblPreyCountReductionPerLarvaHunt[is.na(tblPreyCountReductionPerLarvaHunt) == FALSE]) 
  lGroupHuntStats <- list(nLarva                        = nLarva,
                          vHNablaPreyCount              = tblPreyCountReductionPerLarvaHunt, ##Mean Prey Reduction
                          vHInitialPreyCount            = tblAvailableInitialPreyCountPerLarva,
                          nPreyCountSamples             = nPreyCountSamples, ##COunt of Samples Used to Calc Mean
                          meanNablaPreyCount            = mean(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE),
                          sdNablaPreyCount              = sd(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE),
                          seNablaPreyCount              = sd(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE)/sqrt(nPreyCountSamples),
                          minNablaPreyCount             = min(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE),
                          medianNablaPreyCount          = median(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE),
                          maxNablaPreyCount             = max(tblPreyCountReductionPerLarvaHunt,na.rm=TRUE),
                          ninitPreyCountSamples         = nAvailInitialPreyCount,
                          initPreyCount                 = mean(tblAvailableInitialPreyCountPerLarva,na.rm=TRUE), ##Mean Count At Event 1
                          initsePreyCount               = sd(tblAvailableInitialPreyCountPerLarva,na.rm=TRUE)/sqrt(nAvailInitialPreyCount), ##Mean Count At Event 1
                          totalFrames                   = sum(tblRecDurationPerLarva),
                          vDataSetID                    = levels(datHuntEvent$dataSetID),
                          vExpID                        = levels(datHuntEvent$expID),
                          vLarvaID                      = levels(datHuntEvent$larvaID),
                          vIDLookupTable                = udatLT,
                          vHDurationPerLarva            = tblHuntDurationPerLarva,
                          vmeanHEpisodeDurationPerLarva	= tblMeanEpisodeDurationPerLarva,
                          vHLarvaEventCount             = tblHuntsCounts,
                          groupHuntEvents               = sum(tblHuntsCounts),
                          meanHuntingEventsPerLarva     = mean(tblHuntsCounts),
                          stdHuntingEventsPerLarva      = sd(tblHuntsCounts),
                          seHuntingEventsPerLarva       = sd(tblHuntsCounts)/sqrt(nLarva),
                          medHuntingEventsPerLarva      = median(tblHuntsCounts),
                          maxHuntingEventsPerLarva      = max(tblHuntsCounts),
                          minHuntingEventsPerLarva      = min(tblHuntsCounts),
                          totalHuntFrames               = sum(tblHuntDurationPerLarva),
                          vLarvaHRatio                  = tblHuntRatioPerLarva,
                          meanDuration                  = mean(tblHuntDurationPerLarva),
                          sdDuration                    = sd(tblHuntDurationPerLarva),
                          seDuration                    = sd(tblHuntDurationPerLarva)/sqrt(nLarva),
                          medDuration                   = median(tblHuntDurationPerLarva),
                          meanEpisodeDuration           = mean(tblMeanEpisodeDurationPerLarva,na.rm=TRUE),
                          sdEpisodeDuration             = sd(tblMeanEpisodeDurationPerLarva,na.rm=TRUE), ##Std Dev Of Means
                          seEpisodeDuration             = sd(tblMeanEpisodeDurationPerLarva,na.rm=TRUE)/sqrt(nLarva), ##Std Err Of Means Of Episode Means
                          medEpisodeDuration            = median(tblMedianEpisodeDurationPerLarva,na.rm=TRUE), #Median OF Medians
                          maxEpisodeDuration            = max(tblMaxEpisodeDurationPerLarva,na.rm=TRUE),
                          minEpisodeDuration            = min(tblMinEpisodeDurationPerLarva,na.rm=TRUE),
                          groupHuntRatio                = sum(tblHuntDurationPerLarva)/sum(tblRecDurationPerLarva), ###Need Actual Frame Numbers Here
                          meanHuntRatioOfGroup          = mean(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          stdHuntRatioOfGroup           = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          seHuntRatioOfGroup            = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva)/sqrt(nLarva),
                          medHuntRatioPerLarva          = median(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          maxHuntRatioPerLarva          = max(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          minHuntRatioPerLarva          = min(tblHuntDurationPerLarva/tblRecDurationPerLarva)
  )
  
  
} ### END oF CALC Stat ##




######################################
################################### STAT HUNT EVENTS ###################
# 
# calcHuntStat <- function(datAllFrames,vexpID)
# {
#   message("##Start Hunting Analysis For Group##")
#   lGroupHunting    <- list()
#   lHuntingEvents    <- list()
#   #lHuntingDuration <-list()
#   #vHuntStartFrames <- list()  
#   #vHuntEndFrames <- list()  
#   
#   try(rm("vHuntStartFrames","vHuntEndFrames","vHuntDeltaFrames","lHuntingDuration"),silent=TRUE)
#   
#   
#   nHuntingEventsForLarva   <- 0
#   meanHuntingEvents <- 0;
#   stdDevHuntingEvents <- 0;
#   idx = 0;
#   idxHuntRec = 0;
#   nLarva  	<- length(vexpID)
#   
#   
#   for (i in vexpID)
#   {
#     idx = idx+1;
#     datLarvaFrames <- datAllFrames[datAllFrames$expID == i,]
#     
#     #Used to Identify Experimet Set From Which Data COmes from - Plotted AS different colour
#     DataSetID             <- ifelse(any(names(datLarvaFrames) == "dataSet"),unique(datLarvaFrames$dataSet),0 )
#     stopifnot(DataSetID >= 0 )
#     DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID ) 
#     
#     vEventID = unique((datLarvaFrames$eventID))
#     
#     ##Filter To Keep Only data in range ###
#     datLarvaFrames <- datLarvaFrames[datLarvaFrames$REyeAngle > -35 & datLarvaFrames$REyeAngle <  35 & datLarvaFrames$LEyeAngle > -35 & datLarvaFrames$LEyeAngle <  35,]
#     
#     nTotalHuntFrames = 0;
#     nHuntingEventsForLarva <- 0;
#     meanHuntDuration <- sdHuntDuration  <- medHuntDuration <- maxHuntDuration<-minHuntDuration <- 0
#     
#     lHuntingDuration <-list() ##Empty List Of Event Durations
#     
#     
#     ##Note Hunt Frames Can Come From Separate Events ##
#     ## So we need to loop over them ##
#     for (k in vEventID)
#     {
#       
#       ##Select Hunt Frames of this Event ## 
#       ## Note Criterion - Left Eye Angle >  Min & Right Eye Angle <  Min (Turned Inwards) & Vergence Angle > 40
#       ## Warning - This criterion Is repeated in the plot Function of plotTrackScatterAndDensities - Warning 
#       datHuntFrames    <- datLarvaFrames[datLarvaFrames$eventID == k & 
#                                            datLarvaFrames$REyeAngle < -G_THRESHUNTANGLE &
#                                            datLarvaFrames$LEyeAngle > G_THRESHUNTANGLE & 
#                                            abs(datLarvaFrames$LEyeAngle-datLarvaFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE,]
#       nTotalHuntFrames <- nTotalHuntFrames + length(datHuntFrames$frameN)
#       ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
#       ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
#       ## Detects End Events as +ve differences
#       vHuntDeltaFrames <- diff(datHuntFrames$frameN,lag=1,differences=1)
#       
#       #Find Start Hunt Event the point leading a big gap in frameN, ie reverse diff s[(1+1):10] - s[(1):9]
#       ##ShifHuntDelteFrames
#       vsHuntDeltaFrames     <- (1:(NROW(datHuntFrames$frameN)+2))
#       vHuntStartFrames      <- list()
#       vHuntEndFrames        <- list()
#       lHuntingDuration[[k]] <- 0
#       minHuntDuration       <- 0
#       
#       
#       # Debug #
#       #message(paste(i,".",k,". Detected Hunt #Frames=", NROW(datHuntFrames$frameN) ))
#       
#       
#       ##Consider Only if there are at least 100 frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
#       if  (NROW(datHuntFrames$frameN) > G_MINEPISODEDURATION)
#       {
#         ##Add (Imaginary) Edge Frame  Numbers  (preceding hunting and after last hunting Frame)
#         vsHuntDeltaFrames[1] <- datHuntFrames$frameN[1]-G_MINGAPBETWEENEPISODES
#         vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- datHuntFrames$frameN[NROW(datHuntFrames$frameN)]+G_MINGAPBETWEENEPISODES
#         ##Copy into Shifted right (lagged) Position
#         vsHuntDeltaFrames[2:(NROW(datHuntFrames$frameN)+1)] <- datHuntFrames$frameN
#         ##Do Rev Diff - ie Nabla taking s_n = x_n-x_{n-1} # Detect Start Events as +ve diff
#         vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] <- datHuntFrames$frameN - vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN))]
#         ##\TODO:  Perhaps Throw away Hunts That Began Before Recording Started or End After It Ends - remove Edge Points ##
#         
#         ##Find Hunt Event STARTs ## Needs A gap of At least #MINGAP Between Events To Be a new Episode
#         ##attempting to Remove delay introduced by filtering 
#         vHuntStartFrames <- datHuntFrames$frameN[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
#         ##Find Hunt Event ENDs - Needs To Be atleast G_MINGAPBETWEENEPISODES Away from next start
#         ##Shift (expanded) DFrame vector to pickup on adjacent Hunt-Frame number signifying the end of the Hunt Episode
#         vsHuntDeltaFrames[2:NROW(vsHuntDeltaFrames)] <- vsHuntDeltaFrames[1:NROW(vsHuntDeltaFrames)-1]
#         ##Pickout End frames by selected frameN locations shifted by 1+
#         vHuntEndFrames  <- datHuntFrames$frameN[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
#         ###
#         
#         #message(paste("Found i=",length(vHuntEndFrames)," Hunt End points and k=",length(vHuntStartFrames), "Start Points"))
#         stopifnot(length(vHuntEndFrames) == length(vHuntStartFrames))
#         lHuntingDuration[[k]]   <- vHuntEndFrames - vHuntStartFrames
#         message(paste(" Event Hunt Duration is  ",sum(lHuntingDuration[[k]])," in ",length(vHuntStartFrames)," episodes for expID:",i," eventID:",k))
#         
#         # idxHuntRec = idxHuntRec + 1;
#         # lHuntingEvents[[idxHuntRec]] <- list(expID =rep(i,length(vHuntEndFrames)),
#         #                                      eventID = rep(k,length(vHuntEndFrames)),
#         #                                      dataSetID = rep(DataSetID,length(vHuntEndFrames)),
#         #                                      fileIdx = rep(unique(datHuntFrames$fileIdx),length(vHuntEndFrames)),
#         #                                      startFrame=unlist(vHuntStartFrames),
#         #                                      endFrame=unlist(vHuntEndFrames),
#         #                                      huntScore=rep(0,length(vHuntEndFrames))) ##0 To Be used for Manual Labelling Of Hunting Success
#         # 
#         idxHuntRec = idxHuntRec + 1;
#         lHuntingEvents[[idxHuntRec]] <- data.frame(expID =i,
#                                                    eventID = k,
#                                                    dataSetID = DataSetID,
#                                                    fileIdx = unique(datHuntFrames$fileIdx),
#                                                    startFrame=unlist(vHuntStartFrames),
#                                                    endFrame=unlist(vHuntEndFrames),
#                                                    huntScore=0) ##0 To Be used for Manual Labelling Of Hunting Success
#         
#         if( lHuntingDuration[[k]] < 0) 
#           stop(paste("Negative Hunt Duration detected expID:",i," eventID:",k," Possible Duplicate tracker file" ) ) ##Catch Error / Can be caused by Duplicate Tracked Video - Check )
#         
#         minHuntDuration  <- min(lHuntingDuration[[k]])
#         
#         nHuntingEventsForLarva  <- nHuntingEventsForLarva + length(vHuntStartFrames) ##increment Event Count
#       }else{
#         ##Hunting Episode Insufficient ##	
#         # vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN)+1)] <- 0  ##Copy into Shifted Position
#         warning(paste("Ignoring low Hunt Frame Number",NROW(datHuntFrames$frameN)," of Larva:",i,"event",k) )
#         if (is.na(k) | k < 1 )
#         {
#           stop(paste("Invalid Event Idx:",k) )
#         }
#         minHuntDuration       <-0
#         lHuntingDuration[[k]] <- 0 ##Otherwise Unused Event IDs end up NULL
#       }
#       
#       
#       #### ERROR Hanlders /  DEbug ###
#       
#       if (is.na(minHuntDuration)) 
#       {  
#         message(paste("Warning - Min Episode Duration is  ",minHuntDuration," for expID:",i," eventID:",k))
#         warning(paste("Warning - Min Episode Duration is  ",minHuntDuration," for expID:",i," eventID:",k))
#         # message("** Probably dublicate track files present **")
#       }else if (minHuntDuration < 0) 
#       {  
#         
#         message(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for expID:",i," eventID:",k))
#         warning(paste("Warning - Min Episode Duration is -Ve:",minHuntDuration," for expID:",i," eventID:",k))
#         message("** Probably dublicate track files present **")
#       }
#       
#       if( any( is.null( lHuntingDuration ) ) )  
#         stop(paste("Null Hunting duration for expID:",i," eventID:",k,". Possible Missing Event IDs",  as.character( unlist(which(sapply(lHuntingDuration, is.null) == TRUE)) )) )
#       
#     }##For each Event Of this Larva  ##########  
#     #####
#     vHuntingDuration = unlist(lHuntingDuration) ##Vector With All Event Durations
#     if (length(vHuntingDuration) > 0)
#     {
#       meanHuntDuration <- mean(vHuntingDuration)
#       sdHuntDuration <-   sd(vHuntingDuration)
#       medHuntDuration  <- median(vHuntingDuration)
#       maxHuntDuration  <- max(vHuntingDuration)
#       minHuntDuration  <- min(vHuntingDuration)
#       ##Debug
#     }
#     else
#     {
#       meanHuntDuration = medHuntDuration = sdHuntDuration = maxHuntDuration = minHuntDuration = 0
#     }
#     
#     if (is.nan(meanHuntDuration) | meanHuntDuration < 0)
#     {
#       stop(paste("Invalid Hunt Duration for expID:",i," eventID:",k) )
#     }
#     
#     if (meanHuntDuration == 0)
#     {
#       warning(paste("ZERO Hunt Duration for expID:",i," All events" ) )
#       message(paste("ZERO Hunt Duration for expID:",i," All events") )
#     }
#     
#     ##Total Duration of Hunting For This Larva ##
#     ##+ Number of Hunting Events Counted as points where eyevergence events are more than 300 frames apart  ##
#     lGroupHunting[[idx]] = list(expID=i,
#                                 dataSetID=DataSetID,
#                                 count=nHuntingEventsForLarva,
#                                 huntframes=nTotalHuntFrames,
#                                 vHuntDurations=vHuntingDuration, ## Hunt Duration of each Episode  (All data)    	                        
#                                 meanHDuration=meanHuntDuration, ##Mean Duration for this Larva (per Episode)
#                                 sdHDuration=sdHuntDuration,
#                                 seHDuration=sdHuntDuration/sqrt(nLarva),
#                                 medHDuration=medHuntDuration,
#                                 maxHDuration = maxHuntDuration,
#                                 minHDuration = minHuntDuration,
#                                 totalframes=length(datLarvaFrames$frameN)
#     )
#     
#     
#     #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
#   } ### For Each Larva In Group ###
#   ##### Note That vHuntStartFrames is invalid beyond this point #
#   
#   datGroupHunting  = do.call(rbind,lGroupHunting)
#   datHuntingEvents = do.call(rbind,lHuntingEvents)
#   #  datGroupHunting = as.data.frame(lGroupHunting)
#   ntotalRecFrames               <- unlist(datGroupHunting[,"totalframes"])
#   vHuntEpisodeDuration          <- as.vector(do.call(rbind,lapply(lGroupHunting,"[[","vHuntDurations"))) ## Of All Episodes
#   
#   vHuntRatio                    <- unlist(datGroupHunting[,"huntframes"])/ntotalRecFrames 
#   ##Where LArva Did not Hunt Then Set Ratio to 0 SO this shows up in data
#   vHuntRatio[is.nan(vHuntRatio)== TRUE] <- 0 
#   #vHuntRatio <- vHuntRatio[is.nan(vHuntRatio)== FALSE]  #Remove Nan Entries So The mean Is not Affected / Produced When Total Frames Is 0 for some LArva ID
#   
#   
#   vHuntingDuration          <- unlist(datGroupHunting[,"huntframes"])##Vector with Total for each Larva
#   vmeanLarvaEpisodeHuntingDuration <- unlist(datGroupHunting[,"meanHDuration"]) ## Vector with Mean Episode Duration of eacg Larva
#   vnHuntingEvents           <- unlist(datGroupHunting[,"count"])
#   vDataSetID                <- unlist(datGroupHunting[,"dataSetID"])
#   nHuntingDuration          <- sum(vHuntingDuration)
#   nHuntingEvents            <- sum(vnHuntingEvents);
#   ntotalFrames              <- sum(unlist(datGroupHunting[,"totalframes"]))
#   meanHuntRatio      <- mean(vHuntRatio) ##Group
#   sdHuntRatio        <- sd(vHuntRatio) 
#   medHuntRatio       <- median(vHuntRatio)
#   maxHuntRatio       <- max(vHuntRatio)
#   minHuntRatio       <- min(vHuntRatio)
#   
#   stopifnot(any(meanHuntRatio >= 0 ) )
#   ##Debug
#   message("Hunting Event For This Group:",nHuntingEvents)
#   
#   ### Hunt Duration Per Larva ###    
#   meanLarvaHuntDuration <- mean(vHuntingDuration)
#   sdLarvaHuntDuration   <- sd(vHuntingDuration)
#   medLarvaHuntDuration  <- median(vHuntingDuration)
#   maxLarvaHuntDuration  <- max(vHuntingDuration)
#   minLarvaHuntDuration  <- min(vHuntingDuration)
#   
#   ### Hunt Duration Per Hunting Episode ###
#   meanHuntEpisodeDuration <- mean(vHuntEpisodeDuration)
#   sdHuntEpisodeDuration   <- sd(vHuntEpisodeDuration)
#   medHuntEpisodeDuration  <- median(vHuntEpisodeDuration)
#   minHuntEpisodeDuration  <- min(vHuntEpisodeDuration)
#   maxHuntEpisodeDuration  <- max(vHuntEpisodeDuration)
#   
#   
#   
#   ### Number of Hunting Events Per Larva ##
#   meanHuntingEvents  <- mean(vnHuntingEvents);
#   sdDevHuntingEvents <- sd(vnHuntingEvents);
#   medHuntingEvents   <- median(vnHuntingEvents);
#   maxHuntingEvents   <- max(vnHuntingEvents);
#   minHuntingEvents   <- min(vnHuntingEvents);
#   if (nHuntingEvents > 1)
#   {
#     
#   }
#   else
#   {
#     warning("No Hunting Events Detected in Group!")
#     message("-* No Hunting Events Detected in Group! *- ")
#     #
#     #meanHuntRatio       = 0
#     #sdHuntRatio       = sd(unlist(datGroupHunting[,"huntframes"])/unlist(datGroupHunting[,"totalframes"])) 
#     
#     meanHuntingEvents = 0;
#     sdDevHuntingEvents = 0;
#     
#   }
#   
#   
#   
#   ###Return Results As List of Variable
#   lGroupHuntStats <- list(nLarva=nLarva,
#                           totalFrames         =  ntotalFrames,
#                           vDataSetID          = vDataSetID,
#                           vexpID            = vexpID,
#                           vHDurationPerLarva = vHuntingDuration,
#                           vmeanHEpisodeDurationPerLarva	= vmeanLarvaEpisodeHuntingDuration,
#                           vHLarvaEventCount   = vnHuntingEvents,
#                           vLarvaHRatio        = vHuntRatio,
#                           huntFrames=nHuntingDuration,
#                           meanDuration=meanLarvaHuntDuration,
#                           sdDuration=sdLarvaHuntDuration,
#                           seDuration=sdHuntDuration/sqrt(nLarva),
#                           medDuration             =medLarvaHuntDuration,
#                           maxDuration             = maxLarvaHuntDuration,
#                           minDuration             = minLarvaHuntDuration,
#                           meanEpisodeDuration=meanHuntEpisodeDuration,
#                           sdEpisodeDuration=sdHuntEpisodeDuration,
#                           seEpisodeDuration=sdHuntDuration/sqrt(nLarva),
#                           medEpisodeDuration=medHuntEpisodeDuration,
#                           maxEpisodeDuration = maxHuntEpisodeDuration,
#                           minEpisodeDuration = minHuntEpisodeDuration,
#                           groupHuntRatio=nHuntingDuration/ntotalFrames,
#                           meanHuntRatioPerLarva=meanHuntRatio,
#                           stdHuntRatioPerLarva=sdHuntRatio,
#                           seHuntRatioPerLarva=sdHuntRatio/sqrt(nLarva),
#                           medHuntRatioPerLarva=medHuntRatio,
#                           maxHuntRatioPerLarva=maxHuntRatio,
#                           minHuntRatioPerLarva=minHuntRatio,
#                           groupHuntEvents=nHuntingEvents,
#                           meanHuntingEventsPerLarva=meanHuntingEvents,
#                           stdHuntingEventsPerLarva=sdDevHuntingEvents,
#                           seHuntingEventsPerLarva=sdDevHuntingEvents/sqrt(nLarva),
#                           seHuntingEventsPerLarva=sdDevHuntingEvents/sqrt(nLarva),
#                           medHuntingEventsPerLarva=medHuntingEvents,
#                           maxHuntingEventsPerLarva=maxHuntingEvents,
#                           minHuntingEventsPerLarva=minHuntingEvents,
#                           vHuntingEventsList= lHuntingEvents,
#                           vdHuntingEventsList=datHuntingEvents
#   )
#   
#   message(paste("Number of Hunting Events detected:",nHuntingEvents, " Mean:", meanHuntingEvents, " SD:",sdDevHuntingEvents));
#   message(paste("Hunting Episode Duration :",meanLarvaHuntDuration, " Episode Mean:", meanHuntEpisodeDuration, " SD:",sdHuntEpisodeDuration));
#   
#   message(paste("Logged ",NROW(datHuntingEvents)," hunt events") )
#   
#   return (lGroupHuntStats)
# }
# 









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


