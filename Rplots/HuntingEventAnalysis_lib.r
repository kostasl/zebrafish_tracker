# ################ Hunt Event Detection on 1st pass on  Tracked Data #####
# ## \todo Code to *Detect* Prey Hunting Events - Locate Videos And Start Frames - ###
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
  z=regexpr('.*_tracks',x,T)
  
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
      z <- x
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
    lHuntStat[[i]] <- calcHuntStat3(datHuntEvent)
    
  }
  
  datHuntStat = do.call(rbind,lHuntStat)
  #datHuntStat <- rbindlist(lapply(lHuntStat,alloc.col))
  
  return(datHuntStat)
  
}


#{r combine-dispersion-with-hunt events, echo=FALSE,cache=TRUE,results=FALSE,warning=FALSE}
### Combines information across Detected HuntEvents and their outcome, the dispersion of the larval trajectory during which these were
## initiated  ,and the position at which the hunt-event was initiated
mergeDispersionOntoHuntEvents <- function(datDispersion, datAllFrames, datHuntLabelledEventsSBMerged_fixed)
{
  message(paste(" Loading Hunt Event List to Analyse... "))
  
  start.time <- Sys.time()
  ## Attach Frame Number to Dispersion Data
  datDispersion <- cbind(datDispersion,frameN=datAllFrames[datDispersion[,"frameRow"],'frameN'],
                         posX=datAllFrames[datDispersion[,"frameRow"],'posX'],
                         posY=datAllFrames[datDispersion[,"frameRow"],'posY'])
  
  ## Extract HUnt Event Foraging State
  ### 1. Get Dispersion Measure of each hunt event
  ## MERGE Dispersion With HuntEvent Records
  # Need to convert to char  GroupID factor( datHuntLabelledEventsSBMerged_fixed$groupID) 
  datDispersion$groupID <- levels( datDispersion$groupID)[datDispersion$groupID]
  datHEventDispersion <- merge(datDispersion,datHuntLabelledEventsSBMerged_fixed,
                               all.y=TRUE,
                               by.x=c("frameN","expID","eventID","groupID"),
                               by.y=c("startFrame","expID","eventID","groupID") )
  ## Add Position In Arena
  saveRDS(datHEventDispersion,file=paste0(strDataStore,"/huntEvent_mergedwith_Dispersion",tsec_timeWindow,"sec.rds") )
  message("Saved to:",paste0(strDataStore,"/huntEvent_mergedwith_Dispersion",tsec_timeWindow,"sec.rds") )
  
  nmergemissingEvents <-  NROW(datHEventDispersion[!(datHEventDispersion$frameRow %in% datDispersion$frameRow),] )
  if (nmergemissingEvents > 0)
    warning("** ",nmergemissingEvents, " Hunt Events failed to merge ")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(data.frame(datHEventDispersion ) )
}






######################################################### EXTRACT HUNTING EVENTS ##################################
##Focus on extracting and Identifying Hunting events - Eye Vergence
##Return list of HuntEvents / With start and End Frame / and Video FileName
## TODO : Add Hysterisis in the ON/Off of Eye Vergence Event Detection
detectHuntEvents <- function(datAllGroupFrames,vexpID,ptestCond,vdatasetID)
{  
  message("## Extract Hunting Events For Group ##")
  lGroupHunting    <- list()
  lHuntingEvents    <- list()
  #lHuntingDuration <-list()
  #vHuntStartFrames <- list()  
  #vHuntEndFrames <- list()  
  if (exists("vHuntStartFrames"))
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
    datLarvaFramesRaw <- datAllGroupFrames[which(datAllGroupFrames$expID == i & datAllGroupFrames$testCond == ptestCond),]
    
    #Used to Identify Experimet Set From Which Data COmes from - Plotted AS different colour
    DataSetID             <- ifelse(any(names(datLarvaFramesRaw) == "dataSet"),unique(datLarvaFramesRaw$dataSet),0 )
    stopifnot(DataSetID >= 0 )
    stopifnot(any(is.numeric(DataSetID ) ) )
    DataSetID             <- ifelse(any(is.na(DataSetID)),G_DATASETPALLETSIZE+1,DataSetID ) 
    
    
    vEventID = unique((datLarvaFramesRaw$eventID)) ## +1 To Account for 0 event
    groupID <- unique((datLarvaFramesRaw$group))
    larvaID = unique((datLarvaFramesRaw$larvaID))
    testCond = unique(datLarvaFramesRaw$testCond)
    
    ##Filter To Keep Only data in range ###
    datLarvaFrames <- datLarvaFramesRaw[which(datLarvaFramesRaw$REyeAngle > -G_THRESHCLIPEYEDATA & datLarvaFramesRaw$REyeAngle <  G_THRESHCLIPEYEDATA &
                                                datLarvaFramesRaw$LEyeAngle > -G_THRESHCLIPEYEDATA & datLarvaFramesRaw$LEyeAngle <  G_THRESHCLIPEYEDATA),]
    
    nTotalHuntFrames      <- 0
    nTotalRecordedFrames  <- NROW(datLarvaFrames) ##Will be calculated as diff in frame N for each Event (Each Event FrameN starts from 0 in tracker)
    
    nHuntingEventsForLarva <- 0;
    meanHuntDuration <- sdHuntDuration  <- medHuntDuration <- maxHuntDuration<-minHuntDuration <- 0
    
    lHuntingDuration  <- list() ##Empty List Of Event Durations
    lHuntingIntervals <- list() ##Empty List Of Numbr of Frames Between Hunting Episodes Within An Event
    
    ##Note Hunt Frames Can Come From Separate Events ##
    ## So we need to loop over them ##
    nInitialPrey <- NA ##Count THe Total Prey Count On 1st Event
    nFinalPrey  <- NA
    for (k in vEventID)
    {
        k <- max(k,1) ##Event ID min is 1 / Do not Allow 0
        
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
        
        ##On 1st Event Save The Initial mean Prey Count
        if (k==1 & NROW(datEventFrames) > 0)
        {
          nInitialPrey <- mean(datEventFrames$PreyCount[1:min(NROW(datEventFrames),PREY_COUNT_FRAMEWINDOW)],na.rm = TRUE)
          stopifnot(!is.na(nInitialPrey))
        }
       
         
        ##Count Prey On Final Event
        if ( (k == max(vEventID)) & NROW(datEventFrames) > 0)
        {
          nFinalPrey <- mean(datEventFrames$PreyCount[max(1,NROW(datEventFrames)-PREY_COUNT_FRAMEWINDOW):NROW(datEventFrames)],na.rm = TRUE)
          stopifnot(!is.na(nFinalPrey))
        }
        
        ##Add to Total Number of frames - Using Actual Frame Number of steps Instead of simpler NRow(datLarvaFrames) 
        ##(ideally these should match tho) but tracker can loose object for a few frames and skip
        #nTotalRecordedFrames <- nTotalRecordedFrames + sum(diff(datLarvaFrames[datLarvaFrames$eventID == k,] ))
        
        ##Detect Vergence Consecutive Blocks and find Hunting Event Initiation ##
        ## This COmputes s_n = x_n+1 - x_n # ie SIgnals the end of a consecutive Block
        ## Detects End Events as +ve differences
        vHuntDeltaFrames <- diff(datHuntFrames$frameN,lag=1,differences=1)
        
        #Find Start Hunt Event the point leading a big gap in frameN, ie reverse diff s[(1+1):10] - s[(1):9]
        ##ShifHuntDelteFrames
        vsHuntDeltaFrames      <- (1:(NROW(datHuntFrames$frameN)+2))
        vHuntStartFrames       <- list()
        vHuntEndFrames         <- list()
        lHuntingDuration[[k]]  <- 0
        lHuntingIntervals[[k]] <- 0
        minHuntDuration        <- 0
        
        
        # Debug #
        #message(paste(i,".",k,". Detected Hunt #Frames=", NROW(datHuntFrames$frameN) ))
  
        ##Consider Only if there are at a min of Hunting frames in this Event  - Surragate to proper Hunt Length Filter frames Durations Minimum Duration  / Here MisUsing The Name EpisodeDuration
        if  (NROW(datHuntFrames$frameN) > G_MINEPISODEDURATION)
        {
          ## Add (Imaginary) Edge Frame  Numbers  (preceding hunting and after last hunting Frame)
          vsHuntDeltaFrames[1] <- datHuntFrames$frameN[1]-G_MINGAPBETWEENEPISODES
          vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- datHuntFrames$frameN[NROW(datHuntFrames$frameN)]+G_MINGAPBETWEENEPISODES
          ## Copy into Shifted right (lagged) Position
          vsHuntDeltaFrames[2:(NROW(datHuntFrames$frameN)+1)] <- datHuntFrames$frameN
          ## Do Rev Diff - ie Nabla taking s_n = x_n-x_{n-1} # Detect Start Events as +ve diff
          vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] <- datHuntFrames$frameN - vsHuntDeltaFrames[1:(NROW(datHuntFrames$frameN))]
          ##\TODO:  Perhaps Throw away Hunts That Began Before Recording Started or End After It Ends - remove Edge Points ##
          ##Find Hunt Event STARTs ## Needs A gap of At least #MINGAP Between Events To Be a new Episode
          ##attempting to Remove delay introduced by filtering 
          vHuntStartEyeVergence <- datHuntFrames[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES,"LEyeAngle"]-datHuntFrames[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES,"REyeAngle"]
          vHuntStartFrames      <- datHuntFrames[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES,"frameN" ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
          vHuntStartFrameRowIDs <- rownames(datHuntFrames[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES,])
          ##Find Hunt Event ENDs - Needs To Be atleast G_MINGAPBETWEENEPISODES Away from next start
          ##Shift (expanded) DFrame vector to pickup on adjacent Hunt-Frame number signifying the end of the Hunt Episode
          vsHuntDeltaFrames[2:NROW(vsHuntDeltaFrames)] <- vsHuntDeltaFrames[1:NROW(vsHuntDeltaFrames)-1]
          #Add ENd Edge of DS Vector 
          vsHuntDeltaFrames[NROW(vsHuntDeltaFrames)] <- G_MINGAPBETWEENEPISODES
          ##Pickout End frames by selected frameN locations shifted by 1+ 
          ##End Frames Selected By FrameDIff > Than MinGapBetween Episodes
          ## But it can happen that a single episode is shorter than G_MINGAPBETWEENEPISODES if G_MINEPISODEDURATION < G_MINGAPBETWEENEPISODES
          vHuntEndFrames  <- datHuntFrames$frameN[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
          vHuntEndFrameRowIDs <- rownames( datHuntFrames[vsHuntDeltaFrames[2:(NROW(vsHuntDeltaFrames)-1)+1] >= G_MINGAPBETWEENEPISODES, ])#datHuntFrames$frameN[vsHuntDeltaFrames[1:NROW(datHuntFrames$frameN)] >= G_MINGAPBETWEENEPISODES ]-as.integer(nFrWidth/3) ##Note Ignores hunting event at the very beginning of recording frames for event  (<300 frames away from start)
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

          ##Check If More than Hunt Episode occurred in this Event - And Calculate Interval Between them
          meanInterval         <- NA
          shiftHuntStartFrames <- NA
          if (length(vHuntStartFrames) > 1  )
          {
            ##Calc diff between end of event to next event start : Shift start frame vector to right and shrink End Frames vector from end- 
            shiftHuntStartFrames   <- vHuntStartFrames[2:NROW(vHuntStartFrames)]
            shiftHuntEndFrames     <- vHuntEndFrames[1:NROW(vHuntEndFrames)-1] 
            lHuntingIntervals[[k]] <- shiftHuntStartFrames-shiftHuntEndFrames
            meanInterval           <- mean(lHuntingIntervals[[k]])
            shiftHuntStartFrames  <-  c(shiftHuntStartFrames,NA) ##Add NA On End Which We use for indicating Next Hunt Frame When there is no Event
          }else
            lHuntingIntervals[[k]] <- NA

          message(paste(" Event Hunt Duration :",sum(lHuntingDuration[[k]]), " in ",length(vHuntStartFrames)," t interval between: ", meanInterval, " episodes for expID:",i," eventID:",k))
            
          idxHuntRec = idxHuntRec + 1;
          
          muEpiPreyCount <- mean(datHuntFrames$PreyCount,na.rm = TRUE) ##Mean Prey Count Across Hunt Frames
          
          stopifnot(!is.na(muEpiPreyCount))
          
          lHuntingEvents[[idxHuntRec]] <- data.frame(expID               = factor(i,levels=vexpID),
                                                     eventID             = k,
                                                     dataSetID           = factor(DataSetID,levels=vdatasetID),
                                                     larvaID             = factor(larvaID,levels=seq(1:10)) , ##Identifies Larva Between Empty And Live Test Conditions
                                                     groupID             = groupID,
                                                     testCond            = testCond, ##Latest Additions / Previous ExpID minor decimal indicated condition
                                                     fileIdx             = unique(datHuntFrames$fileIdx),
                                                     eyeVergence         = unlist(vHuntStartEyeVergence),
                                                     startFrame          = unlist(vHuntStartFrames),
                                                     endFrame            = unlist(vHuntEndFrames),
                                                     startFrameRowID     = unlist(vHuntStartFrameRowIDs),
                                                     endFrameRowID       = unlist(vHuntEndFrameRowIDs),
                                                     nextHuntFrame       = unlist(shiftHuntStartFrames), ## Indicates when next Event Starts / Last Hunt Episode will have an NA as next frame
                                                     nExpFrames          = nTotalRecordedFrames,
                                                     InitPreyCount      = nInitialPrey, ##Mean Prey Count Across Hunt PREY_COUNT_FRAMEWINDOW Frames of 1st Event
                                                     FinalPreyCount     = c(rep(NA,length(vHuntEndFrames)-1),nFinalPrey), ##Mean Prey COunt On Last Event's PREY_COUNT_FRAMEWINDOW frames-Add only To Last Record
                                                     PreyCount          =  muEpiPreyCount,
                                                     huntScore          = 0,
                                                     markTracked        = NA, ## Flags that this event has been been closelly Retracked for Analysis
                                                     stringsAsFactors = FALSE) ##0 To Be used for Manual Labelling Of Hunting Success
          
          if( lHuntingDuration[[k]] < 0) 
          {
            stop(paste("Negative Hunt Duration detected expID:",i," eventID:",k," Possible Duplicate tracker file" ) ) ##Catch Error / Can be caused by Duplicate Tracked Video - Check )
          }
          
          minHuntDuration  <- min(lHuntingDuration[[k]])
          
          nHuntingEventsForLarva  <- nHuntingEventsForLarva + length(vHuntStartFrames) ##increment Event Count
        }###(NROW(datHuntFrames$frameN) > G_MINEPISODEDURATION)
        else
        {
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
                                                 testCond   = testCond, ##Latest Additions 
                                                 fileIdx    = 0,
                                                 eyeVergence = 0,
                                                 startFrame = 0,
                                                 endFrame   = 0,
                                                 startFrameRowID     = 0,
                                                 endFrameRowID       = 0,
                                                 nextHuntFrame      = 0,
                                                 nExpFrames         = nTotalRecordedFrames,
                                                 InitPreyCount      = nInitialPrey,
                                                 FinalPreyCount     = nFinalPrey,
                                                 PreyCount          = nmeanPreyCount,
                                                 huntScore          = 0,
                                                 markTracked        = NA, ## Flags that this event has been been closelly Retracked for Analysis
                                                 stringsAsFactors = FALSE) ##0 To Be used for Manual Labelling Of Hunting Success
      
    }
    
    
    #message("Hunting Events For This Larva:",nHuntingEventsForLarva)
    
  } ### For Each Larva In Group ###
  ##### Note That vHuntStartFrames is invalid beyond this point #
  
  datHuntingEvents = do.call(rbind,lHuntingEvents)
  #datHuntingEvents <- rbindlist(lapply(lHuntingEvents,alloc.col))
  
  #stopifnot(any(is.nan(datHuntingEvents) ) ==FALSE )
  #  datGroupHunting = as.data.frame(lGroupHunting)

  nHuntingEvents            <- NROW(datHuntingEvents);
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
## Note:  Aadded Nabla Of Prey Count
calcHuntStat3 <- function(datHuntEvent)
{
  
  message(paste("##V3 Calculate Hunting Statitistics for Group ",unique(datHuntEvent$groupID),"#n",length(unique(datHuntEvent$groupID) )," ##" ) )
  if (NROW(datHuntEvent[is.na(datHuntEvent$groupID) ,] ) > 0 )
  {
    warning("calcHuntStat3: NA found in datHuntEvent GroupID - NA rows removed")
    datHuntEvent <- datHuntEvent[!is.na(datHuntEvent$groupID) ,]
  } 
  stopifnot( length(unique(datHuntEvent$groupID) ) ==1 ) ## Only One Condition Should Be analysed at a time - OtherWise LarvaID may mix results between conditions 
  ## Redo Factor On Subset Of Data - Excluding any ExpID that do not belong Here
  datHuntEvent$expID <- factor(datHuntEvent$expID) 
  ##This Method Produces The Vector WIth The zero Values for Non Hunting Larvae - but mean And SD are correct with aggregate Method
  ##Remove Factor from expID, so as to select only expIDs from this group
  tblHuntDurationPerLarva  <- tapply(datHuntEvent$endFrame-datHuntEvent$startFrame, (datHuntEvent$expID),sum)
  tblHuntDurationPerLarva  <- replace(tblHuntDurationPerLarva,is.na(tblHuntDurationPerLarva),0) ##Replace NA with 0 -Duration
  
  
  ###But  Mean Hunting Duration Per Larva is also correct if using aggregate  ##
  #tblHuntDurationsPerLarva <- aggregate(x=datHuntEvent$endFrame-datHuntEvent$startFrame, by=list(expID=datHuntEvent$expID), FUN=sum)
  
  ## Hunting Duration Per Episode ##
  
  ##Replace 0s with NA - Ignore 0 Duration Episodes
  
  datHuntEventNonZeroEpi <- datHuntEvent[datHuntEvent$endFrame-datHuntEvent$startFrame > 0 &
                                           datHuntEvent$eventID   != 0 ,]
  ##Events that include Some form of Target tracking 
  datHuntEventPursuitEpi <- datHuntEventNonZeroEpi[grepl("Success",datHuntEventNonZeroEpi$huntScore) | 
                                                     grepl("Fail",datHuntEventNonZeroEpi$huntScore) ,]
  ##Calc number of hunting events stat per Experiment (Convert to numeric so as to select only exp from selected GroupID stats)
  tblHuntsCounts<-table((datHuntEventNonZeroEpi$expID) ) 
  tblHuntCaptureAttempt <- table((datHuntEventPursuitEpi$expID) ) 
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
  
  
  ###--- Hunt Intervals ---#####
  # nextHuntFrame can be NA - which will produce wantings for inf when substracting
  nIntervalSamples <- length(datHuntEventNonZeroEpi[!is.na(datHuntEventNonZeroEpi$nextHuntFrame),])
  with(datHuntEventNonZeroEpi[ !is.na(datHuntEventNonZeroEpi$nextHuntFrame),], { 
    tblMeanEpisodeIntervalPerLarva  <<- tapply(nextHuntFrame-endFrame, expID,mean,na.rm=TRUE)
    tblMedianEpisodeIntervalPerLarva <<- tapply(nextHuntFrame-endFrame, expID,median,na.rm=TRUE)
    tblMinEpisodeIntervalPerLarva   <<- tapply(nextHuntFrame-endFrame, expID,min,na.rm=TRUE)
    tblMaxEpisodeIntervalPerLarva   <<- tapply(nextHuntFrame-endFrame, expID,max,na.rm=TRUE)
  })
  #stopifnot((! is.finite(tblMinEpisodeIntervalPerLarva)))
  #### #### ## # # # # 
  
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
  
  # Prey Reduction Per Experiment ##
  ## \todo This does not seems to work Correctly - Better Use the continuously recorded Prey Number
  tblPreyCountReductionPerLarva <-  tapply(datHuntEvent$FinalPreyCount, datHuntEvent$expID,sum, na.rm=TRUE)   
  #tblPreyCountReductionPerLarva <- replace(tblPreyCountReductionPerLarva,is.nan(tblPreyCountReductionPerLarva),NA)
  
  ##Calc Sum Of Hunt Ratios of each Episode, to obtain Hunt Ratio of experiment
  tblHuntRatioPerLarva   <- tapply( (datHuntEvent$endFrame-datHuntEvent$startFrame)/datHuntEvent$nExpFrames, datHuntEvent$expID,sum,na.rm=TRUE) 
  
  
  ##Number of Appearances in Hunt Events - (Minimum Is 1 - so as to obtain LarvaID-ExpID Link)  /  ##
  tblHunts <-table(datHuntEvent$larvaID,datHuntEvent$expID,datHuntEvent$dataSetID)
  
  ##Exp ID LookUp Table So Ican locate the same larva Across Empty->Live Condition
  ## 
  datLT <- data.frame(cbind(larvaID=levels(datHuntEvent$larvaID)[datHuntEvent$larvaID],expID=levels(datHuntEvent$expID)[datHuntEvent$expID],dataSetID=levels(datHuntEvent$dataSetID)[datHuntEvent$dataSetID]))
  udatLT <- unique(datLT)
  
  
  ## COmbine Results Into One Handy Structure containing precalculated results ##
  
  nLarva <- NROW(tblHuntDurationPerLarva)
  nAvailInitialPreyCount <- NROW(tblAvailableInitialPreyCountPerLarva[is.na(tblAvailableInitialPreyCountPerLarva) == FALSE])
  nPreyCountSamples <- NROW(tblPreyCountReductionPerLarvaHunt[is.na(tblPreyCountReductionPerLarvaHunt) == FALSE]) 

  if ( is.na(nPreyCountSamples) | is.nan(nPreyCountSamples) ) 
  {
    stop(paste("nPreyCountSamples not counted correctly : ",nPreyCountSamples) )
  }
  
  lGroupHuntStats <- list(nLarva                        = nLarva,
                          vHNablaPreyCount              = tblPreyCountReductionPerLarvaHunt, ##Mean Prey Reduction
                          vHInitialPreyCount            = tblAvailableInitialPreyCountPerLarva,
                          vHPreyReductionPerLarva       = tblPreyCountReductionPerLarva,
                          nPreyCountSamples             = nPreyCountSamples, ##COunt of Samples Used to Calc Mean
                          meanPreyReductionPerLarva     = mean(tblPreyCountReductionPerLarva,na.rm=TRUE),
                          sePreyReductionPerLarva       = sd(tblPreyCountReductionPerLarva,na.rm=TRUE)/sqrt(nLarva),
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
                          vHLarvaCaptureEventCount      = tblHuntCaptureAttempt,
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
                          seEpisodeDuration             = sd(tblMeanEpisodeDurationPerLarva,na.rm=TRUE)/sqrt(nIntervalSamples), ##Std Err Of Means Of Episode Means
                          medEpisodeDuration            = median(tblMedianEpisodeDurationPerLarva,na.rm=TRUE), #Median OF Medians
                          maxEpisodeDuration            = max(tblMaxEpisodeDurationPerLarva,na.rm=TRUE),
                          minEpisodeDuration            = min(tblMinEpisodeDurationPerLarva,na.rm=TRUE),
                          vmeanHuntInterval             = tblMeanEpisodeIntervalPerLarva,
                          meanHuntInterval              = mean(tblMeanEpisodeIntervalPerLarva,na.rm=TRUE),
                          seHuntInterval                = sd(tblMeanEpisodeIntervalPerLarva,na.rm=TRUE)/sqrt(nIntervalSamples), ##Consider that is should be Mean of Hunt Events so use tblHunt counts
                          minHuntInterval               = min(tblMinEpisodeIntervalPerLarva,na.rm=TRUE),
                          maxHuntInterval               = max(tblMaxEpisodeIntervalPerLarva,na.rm=TRUE),
                          groupHuntRatio                = sum(tblHuntDurationPerLarva)/sum(tblRecDurationPerLarva), ###Need Actual Frame Numbers Here
                          meanHuntRatioOfGroup          = mean(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          stdHuntRatioOfGroup           = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          seHuntRatioOfGroup            = sd(tblHuntDurationPerLarva/tblRecDurationPerLarva)/sqrt(nLarva),
                          medHuntRatioPerLarva          = median(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          maxHuntRatioPerLarva          = max(tblHuntDurationPerLarva/tblRecDurationPerLarva),
                          minHuntRatioPerLarva          = min(tblHuntDurationPerLarva/tblRecDurationPerLarva)
  )
  
  
} ### END oF CALC Stat ##

## Make Hunt Stat For dataFrames of HuntEvents with Mixed GroupID
##Filters out Events labelled as nOn Hunting and produces the hunting Statistics on the rest
## Spliting the statistics of Each Group
makeHuntStat <- function(datHuntEvent)
{
  lHuntStat <- list()
  #groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
  strCondTags <- unique(datHuntEvent$groupID)
  datHuntEvent$huntScore <- convertToScoreLabel( datHuntEvent$huntScore)
  #datHuntEvent$expID <- factor(datHuntEvent$expID) ##Maybe Factor Could resolve issue to count 0 Hunt Events
  for (i in strCondTags)
  {
    message(paste("#### ProcessGroup ",i," ###############"))
    ##ExPORT
    #load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    
    ##Filter Hunt Events ##
    datHuntEventFilt <- datHuntEvent[datHuntEvent$groupID == i,]
    #'These are prefiltered by Score Labels but also selected based on Score in calcHuntStat 
    datHuntEventFiltH <- datHuntEventFilt[datHuntEventFilt$huntScore != "NA" &
#                                           datHuntEventFilt$huntScore != "Not_HuntMode/Delete" &
#                                           datHuntEventFilt$huntScore != "Out_Of_Range" &
#                                           datHuntEventFilt$huntScore != "Duplicate/Overlapping" &
                                           datHuntEventFilt$huntScore != "Near-Hunt State" |
                                           datHuntEventFilt$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events
   
    # Some Larvae only produced Near-Hunt States- And THus no actual HuntEvents / We need to count their rate as 0 however, and not remove them from record
    ##The missing ExpID once Filtering For NearHunt State
    missingExpID <- rownames((table(datHuntEventFilt$expID)))[!(rownames((table(datHuntEventFilt$expID))) %in% rownames((table(datHuntEventFiltH$expID))))] 
    ##Add An Empty Event For Each of the Missing Larvae after filtering non Hunt Events
    for (expID in missingExpID )
    {
      rec = head(datHuntEventFilt[datHuntEventFilt$expID == expID,])##Pick One Event Of this Larva
      rec$eventID = 0
      rec$huntScore = 0
      datHuntEventFiltH <- rbind(datHuntEventFiltH,rec)##Append Empty Place Holder Event
    }
   
   lHuntStat[[i]] <- calcHuntStat3(datHuntEventFiltH)
  }
  
  datHuntStat = do.call(rbind,lHuntStat)#
  
  return(datHuntStat)
}

##After Labelling The the nextHuntFrame Field would need adjustments - as we change frame durations.
##Durations Are Difficult as the begining of a hunt event is not clearly defined (Tirverdi et al ) used a monocular critirion eye vergence criterion
##
fixNextHuntEventFrame <- function(datHuntEvents)
{
  
}




### COmbine HuntEvent DataFrames from datafiles Found in a directory - and save into a merged Data file
### Loads Each HuntEvent RData File Found in A given directory and merges the records,
###which can then be identified by groupID
mergeHuntEventRecords <- function(strSrcDir,strExt = "*.RData")
{
  
  ldatHunt <- list()
  #  groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
  #  strCondTags <- names(groupsrcdatList)
  ldatFiles <- getFileSet("",strSrcDir,strExt)
  
  i <- 1
  for (f in ldatFiles)
  {
    message(paste("#### Load HuntEvent File ",f," ###############"))
    load(file=paste(f,"",sep="" )) ##Save With Dataset Idx Identifier
    
    #datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13),labels=vHuntEventLabels )##Set To NoTHuntMode
    
    ldatHunt[[i]] <- datHuntEvent
    i <- i +1
  }
  
  datHuntEvents = do.call(rbind,ldatHunt)#
  
  nDat <- NROW(unique(datHuntEvents$dataSetID) )
  first <- levels(datHuntEvents$dataSetID)[min(unique(as.numeric(datHuntEvents$dataSetID)) )]
  last <- levels(datHuntEvents$dataSetID)[max(unique(as.numeric(datHuntEvents$dataSetID) ) )]
  
  saveRDS(datHuntEvents,file=paste(strDataExportDir,"/setn",nDat,"-D",first,"-",last,"-HuntEvents-Merged.rds",sep="")  )
}

##
## Returns dataframe with  Experiment ID that link the spontaneous and evoked test conditions coming from the same larva - 
## \Note: This Data organization has now been superseeded - Exp ID now remains the same between test conditions with the addition Of testCond Field Used to distinguish test conditions the 
##
getSpontaneousEvokedExperimentPairs <- function(datHuntStat)
{
  lExpPairs <- list()
  i <- 1
  strCondTags <-  rownames(datHuntStat) #c("LE","LL","NE","NL","DE","DL")
  
  ### Plot Connected Larva Event Counts - To Show Individual Behaviour In Spontaneous Vs Evoked Activity
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    vRegL <- datHuntStat[,"vIDLookupTable"][[gL]]
    vRegE <- datHuntStat[,"vIDLookupTable"][[gE]]
    
    ## Fix Missing LarvaID: REMOVE FActor Field / Set NAs Which Are really ID 5
    vRegL[,"larvaID"] <- as.numeric(vRegL[,"larvaID"])
    vRegE[,"larvaID"] <- as.numeric(vRegE[,"larvaID"])
    vRegL[is.na(vRegL$larvaID),"larvaID"] <- 5 
    vRegE[is.na(vRegE$larvaID),"larvaID"] <- 5
    
    ## Drop Factors To Skip Errors 
    vRegL[,"dataSetID"] <- as.numeric(vRegL[,"dataSetID"])
    vRegE[,"dataSetID"] <- as.numeric(vRegE[,"dataSetID"])
    
    for (k in 1:NROW(vRegE) )
    {
      e <- vRegE[k,]
      i <- i + 1
      EvokedExp <- (vRegL[vRegL$dataSetID == e$dataSetID & vRegL$larvaID == e$larvaID,])
      if (NROW(EvokedExp) == 0)
        next() ##Skip If Matched LiveTest Larva Is not Found
      dexppair <- data.frame(groupID.E = gL, groupID.S = gE, expID.E = as.character(EvokedExp$expID), expID.S = as.character(e$expID) )
      lExpPairs[[i]] <- dexppair
    }
  }
  datExpPairs <- do.call(rbind,lExpPairs)
  
  return(data.frame( datExpPairs) )
}####


## Handle Export of  Detected Hunt Events to File ###
writeHuntEventToFile <- function(datHuntEvent,dataSetsToProcess,groupsrcdatListPerDataSet)
{
  
  if (NROW(datHuntEvent) > 0 )
  {
    message("Writing Hunting Data...")
    ## Recover The File names from fileIdxs
    ## Run It for Each DataSet Idx - Not Sure Of Better Way extracting the relevant filelist
    datHuntEvent$filenames = "." ##Create Field
    tcond =  unique(datHuntEvent$testCond) #Test Condition
    rgroup =  unique(datHuntEvent$groupID) #Test Condition
    for (d in idxDataSet)
    {
      
      ##Get Files Used for This DataSet, and this Condition
      
      filelist <- getVideofilePath(unlist(groupsrcdatListPerDataSet[d][[1]][[rgroup]][[1]]),strVideoFilePath)
      ## Override 
      #filelist <- groupsrcdatListPerDataSet[d][[1]][[tcond]][[1]] 
      
      ##Set File Name
      datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0,]$filenames <- filelist[ datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0 ,]$fileIdx ]
    }
    
    strDataFileName <- paste("setn",NROW(dataSetsToProcess),"-D",dataSetsToProcess[1],"-",dataSetsToProcess[NROW(dataSetsToProcess)],"-HuntEvents-",rgroup,"-",tcond,sep="") ##To Which To Save After Loading
    
    write.csv(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".csv",sep="" ) , row.names=FALSE ) 
    ###Save Hunt Event Data Frame
    
    message(paste(" Exporting to:",strDataFileName))
    ##ExPORT 
    datHuntEvent$groupID = i
    save(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    saveRDS(datHuntEvent,file=paste(strDataExportDir,"/",strDataFileName,".rds",sep="" )) ##So It can be loaded into a custom Named Structure
    
  }else{
    message("No Hunting Event to write!")
  }
  
  
} ##Write Hunting Data




