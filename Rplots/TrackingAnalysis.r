## Processing Of Tracker Eye TrackerData, Using R 
# Kostasl Nov 2017
# 10-12-17 : Fixed Error on counting number of hunting episodes
#            Converted to Adding A Row for empty data files such that larva That were tested but produced no data count towards the mean / The case that a fish is invalid should be actually handled when testing in empty conditions and reject a fish that does nothing
#            Otherwise, a non appearing fish counts towards the group mean since its tested for a fixed amount of time (10mins each fish)
# 14-12-17 :  Added MotionTrajectory Analysis - PathLengths/Speed/Ratio #frames(Speed>0) over All ANalysed Event Frames Of A Larva 
#            
# Consider What the Hunt Ratio Is On a No Show Larva? Currently Set To 0 - 
#TODO: Add Colour Marker of Hunting On Trajectories
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
#strVideoFilePath <- "/media/extStore/ExpData/zebrafish_preycapturesetup/" #Home Desk ubuntu
#strVideoFilePath         <- "/media/kostasl/extStore/ExpData/zebrafish_preycapturesetup/AnalysisSet/"
#strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
#strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData"

## Office PC
#strVideoFilePath  <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/expDataKostas/AnalysisSetAlpha/" 
#strTrackerPath    <- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_9_2_GCC_64bit-Release/"
#strTrackeroutPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/"

## Laptop
setwd("~/Dropbox/Calculations/zebrafishtrackerData/")
strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/"


G_THRESHUNTANGLE         <- 20 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100

nFrWidth                 <- 50

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')));
r <- c(rf(30),"#FF0000");


strDataSetDirectories <- list("./Tracked12-10-17/", ##Dataset 1
                              "./Tracked26-10-17/",
                              "./Tracked02-11-17/",##Dataset 3
                              "./Tracked08-11-17/", #4
                              "./Tracked-T116-11-17/",#5
                              "./Tracked30-11-17/",#6
                              "./Tracked07-12-17/",#7
                              "./Tracked14-12-17/",#8
                              "./Tracked21-12-17/")##Dataset n
strCondR  <- "*.csv"; 

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rf(G_DATASETPALLETSIZE),"#FF0000");

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
#lastDataSet = NROW(strDataSetDirectories)
#firstDataSet = 1
#source("runimportTrackerDataFiles.r")
###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
lastDataSet = NROW(strDataSetDirectories)
firstDataSet = 5

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
source("HuntingEventAnalysis.r")
source("TrajectoryAnalysis.r")
source("labelHuntEvents.r")
########################################
##
### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")
#################################

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
  #vlarvaID = unique(filtereddatAllFrames$larvaID)
  ##Select Larvaof this Group
  
  datAllGroupFrames <- datAllFrames[datAllFrames$group == i,]
  #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
  vlarvaID = unique(datAllGroupFrames$larvaID)
  idxDataSet <- unique(datAllGroupFrames$dataSet)
  
  lHuntStat[[i]] = calcHuntStat(datAllGroupFrames,vlarvaID)
  lMotionStat[[i]] <- calcMotionStat(datAllGroupFrames,vlarvaID)
  
  ##Reconstruct DataSet File List - So As to link fileIdx To Files
  #filelist <- getFileSet("LiveFed/Empty/",strDataSetDirectories[[idxDataSet]])
  #filelist <-  

  ##Combine Hunting Events across fish in this Condition In One
  datHuntEvent = do.call(rbind,lHuntStat[[i]]$vHuntingEventsList )
  
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
      
      datHuntEvent[datHuntEvent$dataSet == d,]$filenames <- filelist[ datHuntEvent[datHuntEvent$dataSet == d,]$fileIdx ]
    }
  
    write.csv(datHuntEvent,file=paste("out/HuntEvents",i,".csv",sep="-" ),row.names=FALSE ) 
    ###Save Hunt Event Data Frame
    strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,".RData",sep="-") ##To Which To Save After Loading
    message(paste(" Exporting to:",strDataFileName))
    ##ExPORT 
    datHuntEvent$GroupID = i
    save(datHuntEvent,file=strDataFileName) ##Save With Dataset Idx Identifier
    
  }else{
    message("No Hunting Event to write!")
  }

  #plotGroupMotion(datAllGroupFrames,lHuntStat[[i]],vlarvaID)
  #######################################################################
  ###  EYE - PLOT Scatter and Eye Densities #####
  strCond = i;
  #source("EyeScatterAndDensities.r")
  #####
} ##Process Stat For Each Condition / Write Hunting Events



par(bg="black")
## Hunt Statistics Summary - Combine Rows ##
datHuntStat = do.call(rbind,lHuntStat)
datMotionStat = do.call(rbind,lMotionStat)

##Motion Plots  - Path Length ##
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
datmean <- unlist(datMotionStat[,"meanPathLength"],use.names = FALSE)
datse <- unlist(datMotionStat[,"sePathLength"],use.names = FALSE)
strtitle <- "Mean Path Length Per Larva"

##Fish With No Hunting Events #
#datHuntStat[,"vLarvaID"]$NE[datHuntStat[,"vHLarvaEventCount"]$NE == 0]
##*All:unlist(datHuntStat[,"vLarvaID"],,use.names=FALSE)[unlist(datHuntStat[,"vHLarvaEventCount"],use.names=FALSE) == 0]
pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

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


dev.off();

#plot(rep(xbarcenters[1],NROW(datMotionStat[[1,"vPathLengths"]]) ),datMotionStat[[1,"vPathLengths"]] )

##Motion Plots - Speed ##
strPlotName = "plots/meanSpeedLarva.pdf"
vDat <- datMotionStat[,"vSpeed"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanSpeed"],use.names = FALSE)
datse <- unlist(datMotionStat[,"seSpeed"],use.names = FALSE)
strtitle <- "Movement Speed Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

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

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

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

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}


dev.off()

################-HUNTING-###################
## Bar Plot Mean Hunting Events Per Animal #
## Common Subtitle info            ########
sampleSize = sum(unlist(datHuntStat[,"nLarva"],use.names = FALSE))
totalFrames = sum(unlist(datHuntStat[,"totalFrames"],use.names = FALSE))
#huntFrames = sum(unlist(datMotionStat[,"huntframes"],use.names = FALSE))

FPS = 420;
strsub = paste("#n=", sampleSize, " #F:",totalFrames,
               "(",format(totalFrames/FPS/60,digits =3),"min)",
               " F_H/F:",format(sum(unlist(datHuntStat[,"huntFrames"] ,use.names = FALSE) )/totalFrames,digits =3) ,
               " #Hunts:",sum(unlist(datHuntStat[,"groupHuntEvents"] ,use.names = FALSE) ),
               collapse=NULL)
##### Done Subtitle ##


#X11()

strPlotName = "plots/meanHuntingEventsPerLarva.pdf"
vDat <- datHuntStat[,"vHLarvaEventCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]

datmean <- unlist(datHuntStat[,"meanHuntingEventsPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntingEventsPerLarva"],use.names = FALSE)
strtitle <- "Hunting Events Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #vDatSetID <- ,-1,vDatSetID[[g]] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}
dev.off()
#######################################
######## EPISODE DURATION ############
#X11()

strPlotName = "plots/meanEpisodeDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanEpisodeDuration"],use.names = FALSE) #Of the Group
#vDat    <- datHuntStat[,"vmeanHLarvaDuration"] # Mean Episode Duration of Each LArva
vDat    <- datHuntStat[,"vmeanHEpisodeDurationPerLarva"]
vDatSetID <- datHuntStat[,"vDataSetID"]
datse   <- unlist(datHuntStat[,"seEpisodeDuration"],use.names = FALSE)
strtitle <- "Mean Duration of each Hunting Episode"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}
dev.off()
#######################################

######## HUNTING DURATION PER LARVA ############
#X11()
strPlotName = "plots/meanHuntDurationPerLarva.pdf"
datmean <- unlist(datHuntStat[,"meanDuration"],use.names = FALSE)
datse   <- unlist(datHuntStat[,"seDuration"],use.names = FALSE)
vDat    <- datHuntStat[,"vHDurationPerLarva"] #Total H Duration Per Larva
vDatSetID <- datHuntStat[,"vDataSetID"]
strtitle <- "Duration of Hunting per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}
dev.off()
#######################################
############### HUNT RATIOn ###############
#X11()

strPlotName = "plots/meanHuntRatioPerLarva.pdf"
datmean <- unlist(datHuntStat[,"meanHuntRatioPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntRatioPerLarva"],use.names = FALSE)
vDat    <- datHuntStat[,"vLarvaHRatio"]
vDatSetID <- datHuntStat[,"vDataSetID"]
strtitle <- "Ratio of Time spent Hunting Over all Frames"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=vDatSetID[[g]],col=rDataset[vDatSetID[[g]]] )
}
dev.off()
#############


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


