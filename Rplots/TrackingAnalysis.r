## Processing Of Tracker Eye TrackerData, Using R 
# Kostasl Nov 2017
# 10-12-17 : Fixed Error on counting number of hunting episodes
#            Converted to Adding A Row for empty data files such that larva That were tested but produced no data count towards the mean / The case that a fish is invalid should be actually handled when testing in empty conditions and reject a fish that does nothing
#            Otherwise, a non appearing fish counts towards the group mean since its tested for a fixed amount of time (10mins each fish)
# 14-12-17 :  Added MotionTrajectory Analysis - PathLengths/Speed/Ratio #frames(Speed>0) over All ANalysed Event Frames Of A Larva 
#            
# Consider What the Hunt Ratio Is On a No Show Larva? Currently Set To 0 - 
#TODO: Add Colour Marker of Hunting On Trajectories
# ## Pio Eykolo Na diaspasoume to Atomo Para mia prokatalipsi ##
library(tools)
library(RColorBrewer);
library("MASS");
#library(hexbin)
rm("temp","subsetDat","TrackerData","frameNAll");

####################

source("labelHuntEvents.r")
########

## GLOBAL VARS ###


## Required Variables - Locations 
# Home Desktop
setwd("/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData")
strVideoFilePath <- "/media/extStore/ExpData/zebrafish_preycapturesetup/" #Home Desk ubuntu
strVideoFilePath         <- "/media/kostasl/extStore/ExpData/zebrafish_preycapturesetup/AnalysisSet/"
strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData"

## Office PC
#setwd("/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/")
#strVideoFilePath  <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/expDataKostas/AnalysisSetAlpha/" 
#strTrackerPath    <- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_9_2_GCC_64bit-Release/"
#strTrackeroutPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo18Jan/"

## Laptop
#setwd("~/Dropbox/Calculations/zebrafishtrackerData/")
#strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
#strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
#strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/"


G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100

nFrWidth                 <- 50 ## Sliding Window Filter Width

rf <- colorRampPalette(rev(brewer.pal(11,'Dark2')));
r <- c(rf(30),"#FF0000");


strDataSetDirectories <- list("./Tracked12-10-17/", ##Dataset 1
                              "./Tracked26-10-17/",
                              "./Tracked02-11-17/",##MDataset 3 -NOTE: Does not Larva ID on File Name 
                              "./Tracked08-11-17/", #4 350fps - Missing a condition WTDryFed3Roti - So removed One Set Larva of Data from other conditions to balance the dataset
                              "./Tracked16-11-17/",#5 400fps - Strict Timer Dataset
                              "./Tracked30-11-17/",#6 420fps
                              "./Tracked07-12-17/",#7
                              "./Tracked14-12-17/",#8
                              "./Tracked21-12-17/",
                              "./Tracked11-01-18/",
                              "./TrackedC18-01-18/")##Dataset n 
strCondR  <- "*.csv"; 

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rf(G_DATASETPALLETSIZE),"#FF0000");

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
lastDataSet = NROW(strDataSetDirectories)
firstDataSet = lastDataSet
source("runimportTrackerDataFiles.r")

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
lastDataSet = NROW(strDataSetDirectories)
firstDataSet = lastDataSet

dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)

datAllSets <- list()
groupsrcdatListPerDataSet <- list() ##Holds File List Per Data Set - Used to Cross Ref FileIdx 
for (i in dataSetsToProcess )
{
  strDataFileName <- paste("setn1_Dataset_",strsplit(strDataSetDirectories[[i]],"/")[[1]][[2]],".RData",sep="") ##To Which To Save After Loading
  
  message(paste("...Loading ",strDataFileName) )
  load(strDataFileName,verbose=TRUE)
  
  groupsrcdatListPerDataSet[[i]] <- groupsrcdatList 
  stopifnot(datAllFrames$dataSet ==  i) ##Identify DataSet
  datAllSets[[i]] <- datAllFrames
}

datAllFrames = do.call(rbind,datAllSets);
##Done LOADING Required DataSets


######################################   START DATA Processing SCRIPT   #####################
#######################################################################################
lHuntStat <- list();
lMotionStat <- list();

####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")

source("TrajectoryAnalysis.r")

source("labelHuntEvents.r")

### TRAJECTORIES Indicating Hunting  - With distinct colour for each larva ####
source("plotTrackScatterAndDensities.r")
##########

strCondTags <- names(groupsrcdatList)
#### Process Files  #####
### Calculates Statistics/Makes Plots 
## Extracts Hunting events and saves them on output csv files / one for each conditiongroup - showing start-end frames and file
for (i in strCondTags)
{
  message(paste("#### ProcessGroup ",i," ###############"))
  subsetDat = groupsrcdatList[[i]]
  strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);

  ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
  #vexpID = unique(filtereddatAllFrames$expID)
  ##Select Larvaof this Group
  
  datAllGroupFrames <- datAllFrames[which(datAllFrames$group == i),]
  #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
  vexpID = unique(datAllGroupFrames$expID)
  idxDataSet <- unique(datAllGroupFrames$dataSet)
  
#  lHuntStat[[i]] = calcHuntStat(datAllGroupFrames,vexpID)
  ##Combine Hunting Events across fish in this Condition In One
  #datHuntEvent = do.call(rbind,lHuntStat[[i]]$vHuntingEventsList )
  
  ## Extract Hunting Events From Data
  datHuntEvent = detectHuntEvents(datAllGroupFrames,vexpID,dataSetsToProcess)
  #lMotionStat[[i]] <- calcMotionStat(datAllGroupFrames,vexpID,dataSetsToProcess)
  
  lHuntStat[[i]] <- calcHuntStat2(datHuntEvent)
  
  stopifnot(length(lHuntStat[[i]]$vHLarvaEventCount) > 0)
  ##Reconstruct DataSet File List - So As to link fileIdx To Files
  #filelist <- getFileSet("LiveFed/Empty/",strDataSetDirectories[[idxDataSet]])
  #filelist <-  

  
  if (NROW(datHuntEvent) > 0 )
  {
    message("Writing Hunting Data...")
    ## Recover The File names from fileIdxs
    ## Run It for Each DataSet Idx - Not Sure Of Better Way extracting the relevant filelist
    datHuntEvent$filenames = "." ##Creat Field
    for (d in idxDataSet)
    {
      
      ##Get Files Used for This DataSet, and this Condition
      filelist <- getVideofilePath(unlist(groupsrcdatListPerDataSet[d][[1]][[i]][[1]]),strVideoFilePath)
      
      ##Set File Name
      
      datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0,]$filenames <- filelist[ datHuntEvent[datHuntEvent$dataSet == d & datHuntEvent$fileIdx != 0,]$fileIdx ]
    }
  
    #write.csv(datHuntEvent,file=paste("out/HuntEvents",i,".csv",sep="-" ),row.names=FALSE ) 
    ###Save Hunt Event Data Frame
    strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,".RData",sep="-") ##To Which To Save After Loading
    message(paste(" Exporting to:",strDataFileName))
    ##ExPORT 
    datHuntEvent$groupID = i
    save(datHuntEvent,file=strDataFileName) ##Save With Dataset Idx Identifier
    
  }else{
    message("No Hunting Event to write!")
  }

  #plotGroupMotion(datAllGroupFrames,lHuntStat[[i]],vexpID)
  #######################################################################
  ###  EYE - PLOT Scatter and Eye Densities #####
  strCond = i;
  #source("EyeScatterAndDensities.r")
  #####
} ##Process Stat For Each Condition / Write Hunting Events




## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)
datMotionStat = do.call(rbind,lMotionStat)


source("plotHuntStat.r") 

#################   MOTION ######################
##Motion Plots  - Path Length ##
par(bg="white")
sampleSize = sum(unlist(datMotionStat[,"nLarva"],use.names = FALSE))
totalFrames = sum(unlist(datMotionStat[,"totalFrames"],use.names = FALSE))
moveFrames = sum(unlist(datMotionStat[,"totalMotionFrames"],use.names = FALSE))

FPS = 420;
strsub = paste("#n=", sampleSize, " #F:",totalFrames,
               "(",format(totalFrames/FPS/60,digits =3),"min)",
               " F_M/F:",format( moveFrames/totalFrames,use.names = FALSE ,digits =3) ,
               " #Events:",sum(unlist(datMotionStat[,"totalEventCount"])),
               collapse=NULL)


strPlotName = "plots/meanPathLengthLarva.pdf"
vDat <- datMotionStat[,"vPathLengths"]
vDatSetID <- datMotionStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]

datmean <- unlist(datMotionStat[,"meanPathLength"],use.names = FALSE)
datse <- unlist(datMotionStat[,"sePathLength"],use.names = FALSE)
strtitle <- "Mean Path Length Per Larva"

##Fish With No Hunting Events #
#datHuntStat[,"vexpID"]$NE[datHuntStat[,"vHLarvaEventCount"]$NE == 0]
##*All:unlist(datHuntStat[,"vexpID"],,use.names=FALSE)[unlist(datHuntStat[,"vHLarvaEventCount"],use.names=FALSE) == 0]
pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))

xbarcenters <- boxplotPerCondition(datMotionStat, datmean,datse,strtitle,strsub,strPlotName,ylim)
vpch = c(0:25,32:127)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vpch[vDatSetID[[g]]],col=rDataset[vDatSetID[[g]]] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags)
dev.off();

#plot(rep(xbarcenters[1],NROW(datMotionStat[[1,"vPathLengths"]]) ),datMotionStat[[1,"vPathLengths"]] )

##Motion Plots - Speed ##
strPlotName = "plots/meanSpeedLarva.pdf"
vDat <- datMotionStat[,"vSpeed"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanSpeed"],use.names = FALSE)
datse <- unlist(datMotionStat[,"seSpeed"],use.names = FALSE)
strtitle <- "Movement Speed Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}
dev.off();

##Motion Ratio -  ##
strPlotName = "plots/moveRatioGroup.pdf"
vDat <- datMotionStat[,"vMovementRatio"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanMotionRatio"])
datse <- unlist(datMotionStat[,"seSpeed"],use.names = FALSE)
strtitle <- "Movement Ratio Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}


dev.off();

##Motion Sinuosity -  ##
strPlotName = "plots/pathSinuosityGroup.pdf"
vDat <- datMotionStat[,"vSinuosity"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanSinuosity"])
datse <- unlist(datMotionStat[,"seSinuosity"],use.names = FALSE)
strtitle <- "Path Sinuosity Ratio Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}

dev.off()


# #Calculate EyeVergence Index
# #(Vectors) Are (RightEyeAngle,LeftEyeAngle)
# vRestCentre = c(-10,10);
# vVerg = c(-30,30)-vRestCentre
# 
# Vidx <- function(x) {CosA=((x-vRestCentre)%*%vVerg)/sqrt(sum((x-vRestCentre)^2))*sqrt(sum((vVerg)^2)); return(CosA)}
# 
# ##Obtain Change In Eye Vergence 
# vDEyeVergence = meanf(diff(LeyeAngleAll-ReyeAngleAll),30);
# #Pick + VE Changes 
# huntFiles = fileIdxAll[vDEyeVergence > 50]
# huntFrames = fileNAll[vDEyeVergence > 50]

#temp[[huntFiles]]
#X11();
#plot(ReyeAll,LeyeAll,cex=.1,xlim=c(-50,20),ylim=c(-40,60),col='blue')

#
#hL=hist(ttL[[1]],breaks=seq(-200,200,length=500),xlim=c(-50,100),mfg=c(1,2,3, 3))
#X11();
#hR=hist(ttR[[1]],breaks=seq(-200,200,length=500),xlim=c(-50,100),mfg=c(1,3,3, 3))


#Mean Filtered TrackerData
#meanf <- function(t,k) {tproc=rep(NA,length(t)); for(i in 1:length(t)) tproc[i]=mean(t[max(1,i-k):i]); return(tproc)}
#meanttL=meanf(TrackerData[1]$EyeLDeg,200);
#meanttR=meanf(TrackerData[1]$EyeRDeg,200);
#plot(meanttR,meanttL,cex=.1,xlim=c(-50,20),ylim=c(-40,60));



###HISTOGRAM 2D
##Try HexBin ##
# Create hexbin object and plot
#eyeDm <- data.frame(ttL[[1]],ttR[[1]]) 
#h <- hexbin(eyeDm)
#plot(h, colramp=rf)


##### OPTION 5: The Hard Way (DIY) #######
# http://stackoverflow.com/questions/18089752/r-generate-2d-histogram-from-raw-data
#nbins <- 25
#x.bin <- seq(floor(min(df[,1])), ceiling(max(df[,1])), length=nbins)
#y.bin <- seq(floor(min(df[,2])), ceiling(max(df[,2])), length=nbins)

#freq <-  as.data.frame(table(findInterval(df[,1], x.bin),findInterval(df[,2], y.bin)))
#freq[,1] <- as.numeric(freq[,1])
#freq[,2] <- as.numeric(freq[,2])

#freq2D <- diag(nbins)*0
#freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

# Normal
#image(x.bin, y.bin, freq2D, col=r)

# Log
#image(x.bin, y.bin, log(freq2D), col=r)

##### Addendum: 2D Histogram + 1D on sides (from Computational ActSci w R) #######
#http://books.google.ca/books?id=YWcLBAAAQBAJ&pg=PA60&lpg=PA60&dq=kde2d+log&source=bl&ots=7AB-RAoMqY&sig=gFaHSoQCoGMXrR9BTaLOdCs198U&hl=en&sa=X&ei=8mQDVPqtMsi4ggSRnILQDw&redir_esc=y#v=onepage&q=kde2d%20log&f=false


