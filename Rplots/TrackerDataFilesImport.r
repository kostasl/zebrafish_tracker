

#### Filter FUNCTIONS Using Window to estimate value at point i using  k surrounding values on a window centred at i###
medianf <- function(t,k) {n=length(t);tproc=rep(NA,n); for(i in (k/2):n) tproc[i]=median(t[max(1,i-k/2): min(n, i+k/2) ],na.rm=TRUE); return(tproc)}
##Note It Returns Mean of Vector values centrered at x , with a window width k
meanf <- function(t,k) {n=length(t);tproc=rep(NA,n); for(i in (k/2):n) tproc[i]=mean(t[max(1,i-k/2): min(n, i+k/2) ],na.rm=TRUE); return(tproc)}

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
  ##KL Custom - Combine Peaks Closer than m Apart
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
importTrackerFilesToFrame <- function(listSrcFiles) {
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
      message(paste(j,". Filtering Data :",  temp[[j]]))
      procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$frameN);
      message(paste("Found #Rec:",  length(TrackerData[[i]][[j]]$frameN) ))
      
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
        datProcessed[[procDatIdx]] = data.frame(LEyeAngle= medianf(TrackerData[[i]][[j]]$EyeLDeg,nFrWidth),
                                                REyeAngle= medianf(TrackerData[[i]][[j]]$EyeRDeg,nFrWidth),
                                                posX = TrackerData[[i]][[j]]$Centroid_X,
                                                posY =TrackerData[[i]][[j]]$Centroid_Y,
                                                BodyAngle = TrackerData[[i]][[j]]$AngleDeg,
                                                ThetaSpine_0 = TrackerData[[i]][[j]]$ThetaSpine_0, ## Angle of 1st Spine Tail Seg
                                                DThetaSpine_1 = TrackerData[[i]][[j]]$DThetaSpine_1, ##Relative Angle Diff Between Next 2nd Tail Seg And The 1st one
                                                DThetaSpine_2 = TrackerData[[i]][[j]]$DThetaSpine_2, ##Consecutive Angle diffs
                                                DThetaSpine_3 = TrackerData[[i]][[j]]$DThetaSpine_3,
                                                DThetaSpine_4 = TrackerData[[i]][[j]]$DThetaSpine_4,
                                                DThetaSpine_5 = TrackerData[[i]][[j]]$DThetaSpine_5,
                                                DThetaSpine_6 = TrackerData[[i]][[j]]$DThetaSpine_6,
                                                DThetaSpine_7 = TrackerData[[i]][[j]]$DThetaSpine_7,
                                                frameN=TrackerData[[i]][[j]]$frameN,
                                                fileIdx=rep(j,Nn),
                                                expID=rep(expID,Nn),
                                                eventID=rep(eventID,Nn),
                                                larvaID=rep(larvaID,Nn),
                                                trackID=rep(trackID,Nn),
                                                group=rep(i,Nn),
                                                PreyCount=meanf(TrackerData[[i]][[j]]$RotiferCount,nFrWidth*8),
                                                countEyeErrors=TrackerData[[i]][[j]]$nFailedEyeDetectionCount,
                                                TailFitError=TrackerData[[i]][[j]]$lastTailFitError
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
                                                expID=expID,
                                                eventID=eventID,
                                                larvaID=larvaID,
                                                trackID=trackID,
                                                group=i,
                                                PreyCount=0,
                                                countEyeErrors=0,
                                                TailFitError=0
                                                );
        message(paste("No Data for ΕχpID",expID,"event ",eventID," larva ",larvaID))
        
      }
      
      
      
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
    message(paste("#### Load Food Track Files Of Group ",i," ###############"))
    
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
      procDatIdx = procDatIdx+1; ##INcreased Count Of Processed Files    
      
      message(paste(j,". Filtering Data :",  temp[[j]]))
      procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$frameN);
      message(paste("Found #Rec:",  length(TrackerData[[i]][[j]]$frameN) ))
      
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
      
      ##FILTER Out NA values - Set to 0
      #message(NROW(TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$EyeLDeg),]))
      if (any(is.na(TrackerData[[i]][[j]]$Centroid_X)) )
      {
        ##Filter Out NA Values
        TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$Centroid_X),]$Centroid_X <- 1000
        message(paste("**NA Values Centroid_X of procDatIdx: ",procDatIdx, " will be set to 1000 in file:" , temp[[j]]))
        warning(paste(" Food Position X  NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      if (any(is.na(TrackerData[[i]][[j]]$Centroid_Y)) )
      {
        ##Filter Out NA Values
        TrackerData[[i]][[j]][is.na(TrackerData[[i]][[j]]$Centroid_Y),]$Centroid_Y <- 1000
        message(paste("**NA Values in Centroid_Y procDatIdx: ",procDatIdx, " will be setto 1000 in file:" , temp[[j]]))
        warning(paste("Food Position Y NA Values in procDatIdx:",procDatIdx,  temp[[j]]))
      }
      
      
      if ( length(TrackerData[[i]][[j]]$frameN) > 1  )
      {       
        Nn <- length(TrackerData[[i]][[j]]$Centroid_X)
        datProcessed[[procDatIdx]] = data.frame(Prey_X= medianf(TrackerData[[i]][[j]]$Centroid_X,nFrWidth),
                                                Prey_Y= medianf(TrackerData[[i]][[j]]$Centroid_Y,nFrWidth),
                                                frameN=TrackerData[[i]][[j]]$frameN,
                                                ROI = TrackerData[[i]][[j]]$ROI,
                                                fileIdx=rep(j,Nn),
                                                expID=rep(expID,Nn),
                                                eventID=rep(eventID,Nn),
                                                larvaID=rep(larvaID,Nn),
                                                trackID=rep(trackID,Nn),
                                                PreyID =TrackerData[[i]][[j]]$foodID,
                                                group=rep(i,Nn)
        );
        
        groupDatIdx = groupDatIdx + 1; ##Count Of Files Containing Data
        
      }   
      else
      {
        ###No Records So add Empty Row - Such that Event And Larva Are on Record
        datProcessed[[procDatIdx]] = data.frame(Prey_X= 180,
                                                Prey_Y= 180,
                                                frameN=0,
                                                fileIdx=j,
                                                expID=expID,
                                                eventID=eventID,
                                                larvaID=larvaID,
                                                trackID=trackID,
                                                PreyID = NA,
                                                group=i
        );
        message(paste("No Data for ΕχpID",expID,"event ",eventID," larva ",larvaID, " TrackNo",trackID))
        
      }
      
      
      
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
  
  
  # Merge two data frames by ID
  datMergedFrames <- merge(x=datHuntEventFrames,y=datAllFrames,by=intersect(names(datHuntEventFrames), names(datAllFrames)),all.x = TRUE )  ## Works like inner join, can set all.x so its a left outer join
  
  ##Return The Food & Fish Merged Frames
  return (datMergedFrames)
  
  
  
}