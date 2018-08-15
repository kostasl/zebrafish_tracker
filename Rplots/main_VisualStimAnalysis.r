## Main Function TO IMport and Process Tracker CSV files for the Visual Stimulation Control Experiment (Rectangular Arenas) ##
#            
# Consider What the Hunt Ratio Is On a No Show Larva? Currently Set To 0 - 
#TODO: Add Colour Marker of Hunting On Trajectories
# ## Pio Eykolo Na diaspasoume to Atomo Para mia prokatalipsi ##
library(tools)
library(RColorBrewer);
library("MASS");
#library(data.table) ##Required for rBindList
#library(hexbin)
rm("temp","subsetDat","TrackerData","frameNAll");

####################




## GLOBAL VARS ###


## Required Variables - Locations 
# Home Desktop
setwd("/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData")
strVideoFilePath  <- "/media/extStore/ExpData/zebrapreyCap/AnalysisVidSet_VSControl/" ##Same As Working Dir
strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/" ##When  Labelling Hunt Events
strTrackInputPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/"
strDatDir        <- "./dat/VSControlTrialA" ##Where Are the Imported RData Stored
strDataExportDir <- "./out/"
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
strTrackeroutPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/" ##For Labelling
#strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from 
strTrackInputPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/"
strDatDir        <- "./dat/VSControlTrialB" ##Where Are the Imported RData Stored
strDataExportDir <- "./out/"

## Laptop
setwd("~/Dropbox/Calculations/zebrafishtrackerData/")
strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/" ##Where to stre the Tracker output csv files when labelling events
strTrackInputPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData"##Where to source the Tracker csv files from 
strDataExportDir <- "./out/"

####################

### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")

source("TrajectoryAnalysis.r")

source("DataLabelling/labelHuntEvents_lib.r")
########


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
G_APPROXFPS              <- 420
G_THRESHUNTANGLE         <- 22 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 45 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 45 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments

nFrWidth                 <- 50 ## Sliding Window Filter Width


##Add Source Directory
strDataSetDirectories <- paste(strTrackInputPath, list(
  "VSControlTrial450fps_03-05-18R", ##Dataset 1
  "VSControlTrial450fps_10-05-18", ##Dataset 2
  "VSControlTrial450fps_17-05-18"
),sep="/")


strCondR  <- "*.csv"; 

#display.brewer.all() to see avaulable options
rfc <- colorRampPalette(rev(brewer.pal(11,'Set3')));
r <- c(rfc(30),"#FF0000");

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rfc(G_DATASETPALLETSIZE),"#FF00AA");

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
 #source("TrackerDataFilesImport.r")

  lastDataSet = NROW(strDataSetDirectories)-2
  firstDataSet = NROW(strDataSetDirectories)-2
  source("importVSTrialTrackerDataFiles.r")

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  firstDataSet = NROW(strDataSetDirectories)
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##oad Frames and HuntStats
  source("loadAllDataSets.r")

  ##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
  ## Calculates HuntEvents And Hunt Statistics On Loaded Data ##
  groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]]
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  source("processLoadedData.r")


  source("plotTrackScatterAndDensities.r")
### Make Eye Phase Space Density Plots ##########
for (i in strCondTags)
{
  message(paste("####  Processing  Group ",i," Eye/Trajectory ###############"))
  subsetDat = groupsrcdatList[[i]]
  strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);
  
  ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
  #vexpID = unique(filtereddatAllFrames$expID)
  ##Select Larvaof this Group
  
  datAllGroupFrames <- datAllFrames[which(datAllFrames$group == i),]
  #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
  vexpID = unique(datAllGroupFrames$expID)
  plotGroupMotion(datAllGroupFrames,lHuntStat[[i]],vexpID)
  #######################################################################
  ###  EYE - PLOT Scatter and Eye Densities #####
  strCond = i;
  #source("EyeScatterAndDensities.r")
  #####
}


######## Make Hunt Statistics ##
source("plotHuntStat.r") 

source("plotMotionStat.r")
