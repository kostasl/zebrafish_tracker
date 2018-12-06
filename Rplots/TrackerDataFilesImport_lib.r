## Initial Import From CSV to Data.frames -

#### Filter FUNCTIONS Using Window to estimate value at point i using  k surrounding values on a window centred at i###
##  Adjust window Size to Vector Length, in case small data samples are provided
## Replace NaN with NA, (this can happen if  calculating mean over NA entries) ,tproc[is.nan(tproc) ] = NA; / If ,na.rm=TRUE Is removed then At LEast A NA is returned By Default
medianf <- function(t,k) {n=length(t);tproc=rep(NA,n); k=min(k,n); for(i in (k/2):n) tproc[i]=median(t[max(1,i-k/2): min(n, i+k/2) ]);  return(tproc)}
##Note It Returns Mean of Vector values centrered at x , with a window width k
meanf <- function(t,k) {n=length(t);tproc=t;k=min(k,n); for(i in (k/2):n) tproc[i]=mean(t[max(1,i-k/2): min(n, i+k/2) ]);  return(tproc)} #tproc=rep(NA,n);

##Make Angles appear as offsets from 0 angle, within range -180 to 180 # Used onr Diff Tail angles
wrapAngle <- function(x)
{
  res <- x %% 360
  if (is.na(x))
    return(0)
  
#  if (abs(res) > abs(360-res))
#    res <- res-360
  
  return(res)
}

##Calculates Change In Degree Angle assumes Circle Of 0-360 Degrees
##Acounts For Circular Nature In THe ANgles going -360 To +360
##Polar Diff   
diffPolar <- function(X)
{
  X <-  X %% 360 
  Y <- rep(0,NROW(X))
  # X[X < 0] <- X[X < 0] + 360
  
  if (NROW(X) < 2 )
  {
    warning("Empty Vector to DiffPolar")
    return (Y)
  }
  
  
  for (i in 2:NROW(X) )
  {
    ##Skip NA
    if (is.na(X[i]) | is.na(X[i-1]) )
      next

    Y[i]<- X[i]-X[i-1]
    
    if (Y[i] > 180) ##Replace with Closest Distance Around Circle
      Y[i] <- X[i]-X[i-1] - 360
    
    if (Y[i] < -180)
      Y[i] <- X[i] - X[i-1] + 360
  }
  
  return(Y)
}


##Fixes Lost Tracking / and Out of Range values by Filling In Gaps with Last known good Value 
clipEyeRange <- function(vEyeAngle,lMin,lMax)
{
  if (NROW(vEyeAngle) < 2)
    return(NA)
  
  ##Min Idx TO Start fROM IS 2 
  idxStart <- max(  min(which(!is.na(vEyeAngle) )   ),2)
  ##Check for Errors, Occuring when all values are NA (usually a very short vector)!
  if (is.infinite(idxStart))
    return(vEyeAngle)
  
  ##Check  Value Before the start one- As it will be propagated forward when values are missing
  if (is.na(vEyeAngle[idxStart-1])) 
    vEyeAngle[idxStart-1] <- lMin

  if (vEyeAngle[idxStart-1] > lMax  ) 
    vEyeAngle[idxStart-1] <- lMax

  if (vEyeAngle[idxStart-1] < lMin  ) 
    vEyeAngle[idxStart-1] <- lMin
  
  
  for (e in idxStart:NROW(vEyeAngle))
  {
    
      
    if (is.na(vEyeAngle[e])) 
      vEyeAngle[e] <- vEyeAngle[e-1]

    if (vEyeAngle[e] == 180)
      vEyeAngle[e] <- vEyeAngle[e-1]
    
    #print(vEyeAngle[e])
    if (vEyeAngle[e] > lMax)
      vEyeAngle[e] <- vEyeAngle[e-1]
    
    if (vEyeAngle[e] < lMin)
      vEyeAngle[e] <- vEyeAngle[e-1]
    # ## This part Was used to remove all sudden Spikes - But simpler solution of Clipping Works Fine
    # ##Check Diff
    # dd <- vEyeAngle[e-1]-vEyeAngle[e]
    # ##If Eyes Move More than 3 degrees per frame
    # if (abs(dd) > 20 )
    # {
    #   ##Interpolate // Find Next Value Close to this one
    #   for (ff in (e):NROW(vEyeAngle))
    #   {
    #     dd2 <- vEyeAngle[e-1]-vEyeAngle[ff]
    #     if (abs(dd2) <= 45 )
    #       break; # Found Next Matching Value - Exit Seach Loop
    # 
    #     vEyeAngle[ff] <- vEyeAngle[e-1] #seq(vEyeAngle[e],vEyeAngle[ff],by=-dd2/(ff-e) )
    #   }
    # }
    ##vEyeAngle[e] <- vEyeAngle[e-1]-vEyeAngle[e-2]dd*0.01
  }
  
  return(vEyeAngle)
}



### Find Peak in 1D vector / Courtesy of https://github.com/stas-g/findPeaks
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  ##KL Custom - Remove Peaks Closer than m Apart 
  pks <- pks[diff(pks) > m]
  pks
}

colTraj <- function(x){
  ids = unique(x);
  rcol = rfc(length(ids)) ##Colour Function Defined in Main File
  
  #z_scl <- (x - min(x, na.rm=T))/(max(x, na.rm=T) - min(x, na.rm=T))
  #return(r[z_scl*length(r)])
  return(rcol[which(ids %in% x)])
}

##Return A List of Filenames for All Exp sets in strsrc list for a specified condition in  strCondDir
getFileSet <- function(strCondDir,strsrc,strCondR = "*.csv")
{
  
  lfileList = list()
  for (i in strsrc)
  {
    lfileList = append(lfileList, list.files(path=paste(i,strCondDir,sep = "/"), pattern=strCondR,full.names = TRUE))

    if (length(unlist(lfileList)) < 1)
      warning(paste("No files found while loading from :",paste(i,strCondDir,sep = "/")))
    
  }
  
  return(unlist(lfileList))
}

##Load Files To Tracker Data and Filter Them Out##
## 
importTrackerFilesToFrame <- function(listSrcFiles,strNameFieldFUN) {
  datProcessed <- list();
  ##CHANGE HASH/ID to select between datasets/groups ##
  strCondTags = names(listSrcFiles);
  
  
  procDatIdx = 1;
  groupDatIdx = 1;
  for (i in strCondTags)
  {
    message(paste("#### Load Data Files Of Group ",i," ###############"))
    
    TrackerData <- list();	
    subsetDat = listSrcFiles[[i]];
    temp <-  unlist(subsetDat[1])
    TrackerData[[i]] = lapply(temp, read.delim)
    nDat = length(TrackerData[[i]])
    
    
    groupDatIdx = 0;
    procDatFrames = 0;
    ## FOR EACH DATA FIle IN Group - Filter Data And combine into Single DataFrame For Group ##
    for (j in 1:nDat)
    {
      message(paste(j,". Filtering Data :",  basename( temp[[j]] ) ) )
      procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$frameN);
      message(paste("Found #Rec:",  length(TrackerData[[i]][[j]]$frameN) ))
      
      ## Extract fields values from filename using function name provided##
      lNameDat <- do.call(strNameFieldFUN, list(temp[[j]]) )
      
      
      
      ## Save standart expected fields for experiments ##
      expID <- lNameDat$expID  #as.numeric(brokenname[[1]][length(brokenname[[1]])-3]);
      eventID <-lNameDat$eventID #   as.numeric(brokenname[[1]][length(brokenname[[1]])-2]);
      trackID = lNameDat$trackID #as.integer( gsub("[^0-9]","",brokenname[[1]][length(brokenname[[1]])])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
      ##Extract Larva ID - Identifies larva in group across food condition - ie which larva in Empty group is the same one in the fed group
      #NOTE: Only Available In files names of more Recent Experiments
      larvaID <-lNameDat$larvaID##as.integer( gsub("[^0-9]","",brokenname[[1]][length(brokenname[[1]])-4]) )

      if(!is.numeric(larvaID)  ) ##Check As it Could Be missing
      {
        larvaID <- NA
        warning(paste("No LarvaID In Filename ",temp[[j]] ) )
      }
       
      #Filter Out Empty Files - ones with less than 300 frames ( ~1 sec of data )
      if (!is.numeric(expID) | !is.numeric(eventID) | is.na(expID) | is.na(eventID)  ) 
      {
        #expID <- j
        #message(paste("Auto Set To expID:",expID))
        stop(paste("Could not extract Larva ID and event ID from File Name ",temp[[j]]))
      }

      stopifnot(!is.na(expID))
      stopifnot(!is.na(eventID))
      
      ##FILTER Out NA values - Set to 0
      #message(NROW(TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$EyeLDeg),]))
      if (any(is.na(TrackerData[[i]][[j]]$EyeLDeg)) )
      {
        ##Filter Out NA Values
        TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$EyeLDeg),]$EyeLDeg <- 1000
        message(paste("**NA Values EyeLDeg of procDatIdx: ",procDatIdx, " will be replaced by 0 in file:" , temp[[j]]))
        warning(paste(" EyeLDeg  NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      if (any(is.na(TrackerData[[i]][[j]]$EyeRDeg)) )
      {
        ##Filter Out NA Values
        TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$EyeRDeg),]$EyeRDeg <- 1000
        message(paste("**NA Values in EyeRDeg procDatIdx: ",procDatIdx, " will be replaced by 0 in file:" , temp[[j]]))
        warning(paste("EyeRDeg NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      
      
      if ( length(TrackerData[[i]][[j]]$frameN) > 1  )
      {       
        Nn <- length(TrackerData[[i]][[j]]$EyeLDeg)
        datProcessed[[procDatIdx]] = data.frame(LEyeAngle= (TrackerData[[i]][[j]]$EyeLDeg),
                                                REyeAngle= (TrackerData[[i]][[j]]$EyeRDeg),##medianf(x,nFrWidth)
                                                posX = TrackerData[[i]][[j]]$Centroid_X,
                                                posY =TrackerData[[i]][[j]]$Centroid_Y,
                                                BodyAngle = TrackerData[[i]][[j]]$AngleDeg,
                                                ThetaSpine_0 = TrackerData[[i]][[j]]$ThetaSpine_0, ## Angle of 1st Spine Tail Seg on Global Coordinates
                                                DThetaSpine_1 = TrackerData[[i]][[j]]$DThetaSpine_1,#sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_1,5),wrapAngle) , ##Relative Angle Diff Between Next 2nd Tail Seg And The 1st one / Call wrap to wrap DAngle between -180 to 180
                                                DThetaSpine_2 = TrackerData[[i]][[j]]$DThetaSpine_2, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_2,5),wrapAngle), ##Consecutive Angle diffs
                                                DThetaSpine_3 = TrackerData[[i]][[j]]$DThetaSpine_3, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_3,5),wrapAngle),
                                                DThetaSpine_4 = TrackerData[[i]][[j]]$DThetaSpine_4, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_4,5),wrapAngle),
                                                DThetaSpine_5 = TrackerData[[i]][[j]]$DThetaSpine_5, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_5,5),wrapAngle),
                                                DThetaSpine_6 = TrackerData[[i]][[j]]$DThetaSpine_6, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_6,5),wrapAngle),
                                                DThetaSpine_7 = TrackerData[[i]][[j]]$DThetaSpine_7, #sapply(medianf(TrackerData[[i]][[j]]$DThetaSpine_7,5),wrapAngle),
                                                frameN=TrackerData[[i]][[j]]$frameN,
                                                fileIdx=rep(j,Nn),
                                                #expID=rep(expID,Nn), ##Now attached via cbind to the lNameDat
                                                #eventID=rep(eventID,Nn), ##From Filename - Sequence # of Event captured during recording 
                                                #larvaID=rep(larvaID,Nn), ##As defined in the filename 
                                                #trackID=rep(trackID,Nn), ##The ID given to the pointtrack from the tracker 
                                                group=rep(i,Nn),
                                                trackletID= TrackerData[[i]][[j]]$fishID,
                                                PreyCount=meanf(TrackerData[[i]][[j]]$RotiferCount,nFrWidth*8),
                                                countEyeErrors=TrackerData[[i]][[j]]$nFailedEyeDetectionCount,
                                                TailFitError=TrackerData[[i]][[j]]$lastTailFitError,
                                                templateScore=TrackerData[[i]][[j]]$templateScore
                                                );
        groupDatIdx = groupDatIdx + 1; ##Count Of Files Containing Data
        
      }   
      else
      {
        ###No Records So add Empty Row - Such that Event And Larva Are on Record
        datProcessed[[procDatIdx]] = data.frame(LEyeAngle= 180,
                                                REyeAngle= 180,
                                                posX = 0,
                                                posY = 0,
                                                BodyAngle = 0,
                                                ThetaSpine_0 = 0, ## Angle of 1st Spine Tail Seg
                                                DThetaSpine_1 = 0, ##Relative Angle Diff Between Next 2nd Tail Seg And The 1st one
                                                DThetaSpine_2 = 0, ##Consecutive Angle diffs
                                                DThetaSpine_3 = 0,
                                                DThetaSpine_4 = 0,
                                                DThetaSpine_5 = 0,
                                                DThetaSpine_6 = 0,
                                                DThetaSpine_7 = 0,
                                                frameN=0,
                                                fileIdx=j,
                                                #expID=expID, ##Now attached via cbind to the lNameDat
                                                #eventID=eventID,
                                                #larvaID=larvaID,
                                                #trackID=trackID,
                                                group=i,
                                                trackletID=0,
                                                PreyCount=0,
                                                countEyeErrors=0,
                                                TailFitError=0,
                                                templateScore=0.0
                                                );
        
        message(paste("No Data for ΕχpID",expID,"event ",eventID," larva ",larvaID))
        
      }
      
      ## Attach The FileName extracted Data, to the data frame
      datProcessed[[procDatIdx]] <- cbind(lNameDat,datProcessed[[procDatIdx]])
      
      
      ## Report NA Values ##
      if (any(is.na(datProcessed[[procDatIdx]] ) ))
      {
        message(paste("**NA Values from TrackData ",i,j," Still Present in datProcessed procDatIdx: ",procDatIdx, " in file:" , temp[[j]]))
        warning(paste("NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      #stopifnot(is.numeric(datProcessed[[procDatIdx-1]]$frameN ))
      
      procDatIdx = procDatIdx+1; ##INcreased Count Of Processed Files    
      
    } ##For Each File In Group ##
    datAllFrames = do.call(rbind,datProcessed);
    #datAllFrames <- rbindlist(datProcessed ) 
    
    #xxList = list(LEyeAngle=unlist(ttEyeAnglesAll[,1]),REyeAngle=unlist(ttEyeAnglesAll[,2]),frameN=unlist(ttEyeAnglesAll[,3]),fileIdx=unlist(ttEyeAnglesAll[,4]))
    message(paste("Non empty Data files Count :",  groupDatIdx, " total Frames :",procDatFrames));
    
    message(paste("###### Finished Loading Data Files Of Group ",i," ###############"))
    
  } ##For Each Group Tag ##
  message("#### Importing Fish Tracks Complete -Return Data.frame###")
  message(paste("Total Usable Data files Count :",  procDatIdx, " total Frames :",procDatFrames));
  
  
  return (datAllFrames)
}##END OF IMPORT FUNCTION

##################################
###TODO: Finds Matching food tracker file for dataframe
## Import food data, and merges, using frameNumber, experiment, EventID as keys/
mergeFoodTrackerFilesToFrame <- function(listSrcFoodFiles,datHuntEventFrames) {
  
  datProcessed <- list();
  ##CHANGE HASH/ID to select between datasets/groups ##
  strCondTags = names(listSrcFoodFiles);
  
  
  procDatIdx = 0;
  groupDatIdx = 0;
  for (i in strCondTags)
  {
    message(paste("#### Load Prey Track Files Of Group ",i," ###############"))
    
    TrackerData <- list();	
    subsetDat = listSrcFoodFiles[[i]];
    temp <-  unlist(subsetDat[1])
    TrackerData[[i]] = lapply(temp, read.delim)
    nDat = length(TrackerData[[i]])

    groupDatIdx = 0;
    procDatFrames = 0;
    ## FOR EACH DATA FIle IN Group - Filter Data And combine into Single DataFrame For Group ##
    for (j in 1:nDat)
    {

      message(paste(j,". Filtering Prey Data :",  temp[[j]]))
      procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$FrameN);
      message(paste("Found Prey #Rec:",  length(TrackerData[[i]][[j]]$FrameN) ))
      
      ##Extract Experiment ID
      brokenname = strsplit(temp[[j]],"_")
      expID =  as.numeric(brokenname[[1]][length(brokenname[[1]])-3]);
      eventID = as.numeric(brokenname[[1]][length(brokenname[[1]])-2]);
      trackID = as.integer( gsub("[^0-9]","",brokenname[[1]][length(brokenname[[1]])])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
      

      ##Extract Larva ID - Identifies larva in group across food condition - ie which larva in Empty group is the same one in the fed group
      #NOTE: Only Available In files names of more Recent Experiments
      larvaID <- as.integer( gsub("[^0-9]","",brokenname[[1]][length(brokenname[[1]])-4]) )
      if(!is.numeric(larvaID)  ) ##Check As it Could Be missing
      {
        larvaID <- NA
        warning(paste("No LarvaID In Filename ",temp[[j]] ) )
      }
      
      #Filter Out Empty Files - ones with less than 300 frames ( ~1 sec of data )
      if (!is.numeric(expID) | !is.numeric(eventID) | is.na(expID) | is.na(eventID)  ) 
      {
        #expID <- j
        #message(paste("Auto Set To expID:",expID))
        stop(paste("Could not extract Larva ID and event ID from File Name ",temp[[j]]))
      }
      
      stopifnot(!is.na(expID))
      stopifnot(!is.na(eventID))
      
      
      ## Separate Data For Each Prey ID recorded
      vTrackPreyID <- unique(TrackerData[[i]][[j]]$FoodID)
      for (p in vTrackPreyID)
      {
        procDatIdx = procDatIdx+1; ##INcreased Count Of Processed Prey Track Data
       
        message(paste("#Prey ID:",p ) )
         
        datPreyTracks <- TrackerData[[i]][[j]][TrackerData[[i]][[j]]$FoodID == p,]
        
        ##FILTER Out NA values - Set to 0
        #message(NROW(TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$EyeLDeg),]))
        if (any(is.na(datPreyTracks$Centroid_X)) )
        {
          ##Filter Out NA Values
          datPreyTracks[is.na(datPreyTracks$Centroid_X),]$Centroid_X <- 1000
          message(paste("**NA Values Centroid_X of procDatIdx: ",procDatIdx, " will be set to 1000 in file:" , temp[[j]]))
          warning(paste(" Food Position X  NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
        }
        if (any(is.na(datPreyTracks$Centroid_Y)) )
        {
          ##Filter Out NA Values
          datPreyTracks[is.na(datPreyTracks$Centroid_Y),]$Centroid_Y <- 1000
          message(paste("**NA Values in Centroid_Y procDatIdx: ",procDatIdx, " will be setto 1000 in file:" , temp[[j]]))
          warning(paste("Food Position Y NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
        }
        
        
        if ( length(datPreyTracks$FrameN) > 1  )
        {       
          Nn <- length(datPreyTracks$Centroid_X)
          datProcessed[[procDatIdx]] = data.frame(Prey_X= medianf(datPreyTracks$Centroid_X,nFrWidth),
                                                  Prey_Y= medianf(datPreyTracks$Centroid_Y,nFrWidth),
                                                  Prey_Radius= medianf(datPreyTracks$Radius,nFrWidth),
                                                  frameN=as.numeric(datPreyTracks$FrameN),
                                                  inactiveFrames = datPreyTracks$InactiveFrames,
                                                  ROI = datPreyTracks$ROI,
                                                  fileIdx=rep(j,Nn),
                                                  expID=rep(expID,Nn),
                                                  eventID=rep(eventID,Nn),
                                                  larvaID=rep(larvaID,Nn),
                                                  trackID=rep(trackID,Nn),
                                                  PreyID =datPreyTracks$FoodID,
                                                  group=rep(i,Nn)
          );
          
          groupDatIdx = groupDatIdx + 1; ##Count Of Files Containing Data
          
        }   
        else
        {
          ###No Records So add Empty Row - Such that Event And Larva Are on Record
          datProcessed[[procDatIdx]] = data.frame(Prey_X= 180,
                                                  Prey_Y= 180,
                                                  Prey_Radius= 0,
                                                  frameN=0,
                                                  inactiveFrames=0,
                                                  ROI =0,
                                                  fileIdx=j,
                                                  expID=expID,
                                                  eventID=eventID,
                                                  larvaID=larvaID,
                                                  trackID=trackID,
                                                  PreyID = NA,
                                                  group=i
          );
          message(paste("No  Prey Track Data for ΕχpID",expID,"event ",eventID," larva ",larvaID, " TrackNo",trackID))
          
        }
      
      }#### For Each Prey Id In Food File
      
      ## Report NA Values ##
      if (any(is.na(datProcessed[[procDatIdx]] ) ))
      {
        message(paste("**NA Values from TrackData ",i,j," Still Present in datProcessed procDatIdx: ",procDatIdx, " in file:" , temp[[j]]))
        warning(paste("NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      #stopifnot(is.numeric(datProcessed[[procDatIdx-1]]$frameN ))
      

      
    } ##For Each File In Group ##
    datAllFrames = do.call(rbind,datProcessed);
    #datAllFrames <- rbindlist(datProcessed ) 
    
    #xxList = list(LEyeAngle=unlist(ttEyeAnglesAll[,1]),REyeAngle=unlist(ttEyeAnglesAll[,2]),frameN=unlist(ttEyeAnglesAll[,3]),fileIdx=unlist(ttEyeAnglesAll[,4]))
    message(paste("Non empty Data files Count :",  groupDatIdx, " total Frames :",procDatFrames));
    
    message(paste("###### Finished Loading Data Files Of Group ",i," ###############"))
    
  } ##For Each Group Tag ##
  message("#### Loading Food Data Complete -Return Data.frame###")
  message(paste("Total Usable Data files Count :",  procDatIdx, " total Frames :",procDatFrames));
  
  datAllFrames$frameN       <- as.numeric(datAllFrames$frameN) ##Convert to number so all leading zeros on frameN are ignored
  datHuntEventFrames$frameN <- as.numeric(datHuntEventFrames$frameN)
  # Merge two data frames by ID
  datMergedFrames <- merge(x=datHuntEventFrames,y=datAllFrames,by=intersect(names(datHuntEventFrames), names(datAllFrames)),all.x = TRUE )  ## Works like inner join, can set all.x so its a left outer join
  
  ##Return The Food & Fish Merged Frames
  return (datMergedFrames)
  
  
  
}


#/// Returns a list of name value pairs extracted from TrackerFile name used for the Hunting  Assay
extractFileNameParams_huntingExp <- function(strFileName)
{
  ##Extract Experiment ID
  basename <- basename(strFileName)
  brokenname = unlist(strsplit(basename,"_"))
  expID <-  as.numeric(brokenname[4]);
  eventID <- as.numeric(brokenname[5]);
  larvaID <- as.numeric(gsub("[^0-9]","",brokenname[3]) );
  trackID <- as.integer( gsub("[^0-9]","",brokenname[length(brokenname)])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
  fps     <- as.integer( gsub("[^0-9]","",brokenname[1])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
 # timeMin <- as.integer( gsub("[^0-9]","",brokenname[5])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
  
  return(list(expID=expID,eventID=eventID,larvaID=larvaID,fps=fps) )
}

#/// Returns a list of name value pairs extracted from TrackerFile name used for the PreyCount Feeding Assay
extractFileNameParams_preycountExp <- function(strFileName)
{
  ##Extract Experiment ID
  basename <- basename(strFileName)
  brokenname = unlist(strsplit(basename,"_"))
  expID <-  as.numeric(brokenname[6]);
  eventID <- as.numeric(brokenname[7]);
  larvaID <- as.numeric(gsub("[^0-9]","",brokenname[4]) );
  trackID <- as.integer( gsub("[^0-9]","",brokenname[length(brokenname)])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
  fps     <- as.integer( gsub("[^0-9]","",brokenname[2])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
  timeMin <- as.integer( gsub("[^0-9]","",brokenname[5])  ) ##Extract the Track Sequence In The filename Given Automatically By the tracker , when a file already exists
  
  return(list(expID=expID,eventID=eventID,larvaID=larvaID,trackID=trackID,time=timeMin,fps=fps) )
}
