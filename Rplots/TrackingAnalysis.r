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


####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")

source("TrajectoryAnalysis.r")

source("labelHuntEvents.r")
########

## GLOBAL VARS ###


## Required Variables - Locations 
# Home Desktop
#setwd("/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData")
#strVideoFilePath  <- "/media/extStore/ExpData/zebrapreyCap/AnalysisSet/"
#strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
#strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData"
#strTrackInputPath <- "./" ##Same As Working Dir

##Emily ##
#setwd("/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData")
#strVideoFilePath         <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/expDataEmily/"
#strTrackeroutPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/"
#strDataSetDirectories <- list("/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackedEm/")##Dataset n 
#strDataSetDirectories <- list("/media/extStore/ExpData/zebrapreyCap/TrackedEmilyDat/")##Dataset n 


## Office PC
setwd("/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/")
strVideoFilePath  <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/expDataKostas/AnalysisSetAlpha/" 
strTrackerPath    <- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_9_2_GCC_64bit-Release/"
strTrackeroutPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo01Feb/"
strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from 


## Laptop
#setwd("~/Dropbox/Calculations/zebrafishtrackerData/")
#strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
#strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
#strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/" ##Where to stre the Tracker output csv files when labelling events
#strTrackInputPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData"##Where to source the Tracker csv files from 

G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100

nFrWidth                 <- 50 ## Sliding Window Filter Width


strDataSetDirectories <- paste(strTrackInputPath, list(
                              "/Tracked12-10-17/", ##Dataset 1
                              "/Tracked26-10-17/",
                              "/Tracked02-11-17/",##MDataset 3 -NOTE: Does not Larva ID on File Name 
                              "Tracked08-11-17/", #4 350fps - Missing a condition WTDryFed3Roti - So removed One Set Larva of Data from other conditions to balance the dataset
                              "/Tracked16-11-17/",#5 400fps - Strict Timer Dataset
                              "/Tracked30-11-17/",#6 420fps
                              "/Tracked07-12-17/",#7
                              "/Tracked14-12-17/",#8
                              "Tracked21-12-17/",
                              "/Tracked11-01-18/",
                              "/Tracked18-01-18/",
                              "/Tracked25-01-18/",
                              "/Tracked01-02-18/",
                              "/Tracked08-02-18/",
                              "/Tracked15-02-18/"##Dataset n 
                              ),sep="/")
##Add Source Directory

strCondR  <- "*.csv"; 

#display.brewer.all() to see avaulable options
rf <- colorRampPalette(rev(brewer.pal(11,'Set3')));
r <- c(rf(30),"#FF0000");

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rf(G_DATASETPALLETSIZE),"#FF00AA");

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
  
  lastDataSet = NROW(strDataSetDirectories)-2
  firstDataSet = NROW(strDataSetDirectories)-2
  source("runimportTrackerDataFiles.r")

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  firstDataSet = NROW(strDataSetDirectories)-10
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##oad Frames and HuntStats
  source("loadAllDataSets.r")

  ## Calculates HuntEvents And Hunt Statistics On Loaded Data ##
  source("processLoadedData.r")


### Make Eye Phase Space Density Plots ##########
for (i in strCondTags)
{
  message(paste("#### Eye ProcessGroup ",i," ###############"))
  subsetDat = groupsrcdatList[[i]]
  strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);
  
  ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
  #vexpID = unique(filtereddatAllFrames$expID)
  ##Select Larvaof this Group
  
  datAllGroupFrames <- datAllFrames[which(datAllFrames$group == i),]
  #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
  vexpID = unique(datAllGroupFrames$expID)
  #plotGroupMotion(datAllGroupFrames,lHuntStat[[i]],vexpID)
  #######################################################################
  ###  EYE - PLOT Scatter and Eye Densities #####
  strCond = i;
  source("EyeScatterAndDensities.r")
  #####
}


######## Make Hunt Statistics ##
source("plotHuntStat.r") 




#### LABEL MANUALLY THE HUNT EVENTS WITH THE HELP OF THE TRACKER ###
gc <- "DL"
strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
message(paste(" Loading Hunt Events: ",strDataFileName))
##ExPORT 
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
labelHuntEvents(datHuntEvent,strVideoFilePath,strTrackerPath )
save(datHuntEvent,file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
##########################




