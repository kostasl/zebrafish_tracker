

#### FUNCTIONS ###
medianf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=median(t[max(1,i-k):i]); return(tproc)}
meanf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=mean(t[max(1,i-k):i]); return(tproc)}


colTraj <- function(x){
  ids = unique(x);
  rcol = rf(length(ids))
  
  #z_scl <- (x - min(x, na.rm=T))/(max(x, na.rm=T) - min(x, na.rm=T))
  #return(r[z_scl*length(r)])
  return(rcol[which(ids %in% x)])
}

##Return A List of Filenames for All Exp sets in strsrc list for a specified condition in  strCondDir
getFileSet <- function(strCondDir,strsrc)
{
  strCondR = "*.csv"
  lfileList = list()
  for (i in strsrc)
  {
    lfileList = append(lfileList, list.files(path=paste(i,strCondDir,sep = "/"), pattern=strCondR,full.names = TRUE))

    if (length(unlist(lfileList)) < 1)
      stop(paste("No files found while loading from :",paste(i,strCondDir,sep = "/")))
    
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
    
    
    groupDatIdx = 1;
    procDatFrames = 0;
    ## FOR EACH DATA FIle IN Group - Filter Data And combine into Single DataFrame For Group ##
    for (j in 1:nDat)
    {
      message(paste(j,". Filtering Data :",  temp[[j]]))
      procDatFrames = procDatFrames + length(TrackerData[[i]][[j]]$frameN);
      message(paste("Found #Rec:",  length(TrackerData[[i]][[j]]$frameN) ))
      
      ##Extract Larva ID
      brokenname = strsplit(temp[[j]],"_")
      larvaID =  as.numeric(brokenname[[1]][length(brokenname[[1]])-3]);
      eventID = as.numeric(brokenname[[1]][length(brokenname[[1]])-2]);
      
      #Filter Out Empty Files - ones with less than 300 frames ( ~1 sec of data )
      if (!(is.numeric(larvaID) || !is.numeric(eventID) || !is.na(larvaID) || !is.na(eventID)  ) ) 
      {
        stop(paste("Could not extract Larva ID and event ID from File Name ",temp[[j]]))
      }
      
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
        datProcessed[[procDatIdx]] = data.frame(LEyeAngle= medianf(TrackerData[[i]][[j]]$EyeLDeg,nFrWidth),
                                                REyeAngle= medianf(TrackerData[[i]][[j]]$EyeRDeg,nFrWidth),
                                                posX = TrackerData[[i]][[j]]$Centroid_X,
                                                posY =TrackerData[[i]][[j]]$Centroid_Y,
                                                frameN=TrackerData[[i]][[j]]$frameN,
                                                fileIdx=rep(j,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                                larvaID=rep(larvaID,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                                eventID=rep(eventID,length(TrackerData[[i]][[j]]$EyeLDeg)),
                                                trackID=rep(eventID,length(TrackerData[[i]][[j]]$fishID)),
                                                group=rep(i,length(TrackerData[[i]][[j]]$EyeLDeg)  )
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
                                                frameN=0,
                                                fileIdx=j,
                                                larvaID=larvaID,
                                                eventID=eventID,
                                                trackID=0,
                                                group=i );
        message(paste("No Data for Larva",larvaID,"event ",eventID))
        
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
    #xxList = list(LEyeAngle=unlist(ttEyeAnglesAll[,1]),REyeAngle=unlist(ttEyeAnglesAll[,2]),frameN=unlist(ttEyeAnglesAll[,3]),fileIdx=unlist(ttEyeAnglesAll[,4]))
    message(paste("Non empty Data files Count :",  groupDatIdx, " total Frames :",procDatFrames));
    
    message(paste("###### Finished Loading Data Files Of Group ",i," ###############"))
    
  } ##For Each Group Tag ##
  message("#### Loading Complete -Return Data.frame###")
  message(paste("Total Usable Data files Count :",  procDatIdx, " total Frames :",procDatFrames));
  
  
  return (datAllFrames)
}##END OF IMPORT FUNCTION
