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
#library(data.table) ##Required for rBindList
#library(hexbin)
rm("temp","subsetDat","TrackerData","frameNAll");

####################




## GLOBAL VARS ###


## Required Variables - Locations 
# Home Desktop
setwd("/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData")
strVideoFilePath  <- "/media/extStore/ExpData/zebrapreyCap/AnalysisSet/"
strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/"
strTrackInputPath <- "/media/extStore//kostasl/Dropbox/Calculations/zebrafishtrackerData/" ##Same As Working Dir
strDatDir        <- "./dat/TrackedSessionA" ##Where Are the Imported RData Stored

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
strTrackeroutPath <- "/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/"
strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from 
strDatDir        <- "./dat/TrackedSessionA" ##Where Are the Imported RData Stored

## Laptop
setwd("~/Dropbox/Calculations/zebrafishtrackerData/")
strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_UpTo21Dec/" ##Where to stre the Tracker output csv files when labelling events
strTrackInputPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData"##Where to source the Tracker csv files from 


####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")

source("TrajectoryAnalysis.r")

source("DataLabelling/labelHuntEvents_lib.r")
########


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
G_APPROXFPS              <- 420
G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments

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
                              "/Tracked15-02-18/",
                              "/Tracked22-02-18/"##Dataset n 
                              ),sep="/")
##Add Source Directory


## These Were For VS Control ###


strDataSetDirectories <- paste(strTrackInputPath, list(
  "VSControlTrial450fps_03-05-18", ##Dataset 1
  "VSControlTrial450fps_10-05-18" ##Dataset 2
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
  
  lastDataSet = NROW(strDataSetDirectories)-1
  firstDataSet = 1
  source("runimportTrackerDataFiles.r")

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  firstDataSet = NROW(strDataSetDirectories)-11
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##oad Frames and HuntStats
  source("loadAllDataSets.r")

  ##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
  ## Calculates HuntEvents And Hunt Statistics On Loaded Data ##
  groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]]
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
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


###
  source("plotMotionStat.r")

#### LABEL MANUALLY THE HUNT EVENTS WITH THE HELP OF THE TRACKER ###
gc <- "DL"
firstDataSet = NROW(strDataSetDirectories)-11
lastDataSet = NROW(strDataSetDirectories)
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)

strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
message(paste(" Loading Hunt Events: ",strDataFileName))
##ExPORT 
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEvent <- labelHuntEvents(datHuntEvent,strDataFileName,strVideoFilePath,strTrackerPath,strTrackeroutPath )
##Saving is done in labelHuntEvent on Every loop - But repeated here
save(datHuntEvent,file=paste(strDataFileName,"-backup.RData",sep="" )) ##Save With Dataset Idx Identifier
message(paste("Saved Backup :",strDataFileName,"-backup.RData",sep="") )
##########################


##How to Summarize Success / Fail Scores :
gc <- "LL"
strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12),labels=vHuntEventLabels )##Set To NoTHuntMode
message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
tblLLStat <- table(datHuntEvent$huntScore)

nFailLL <- tblLLStat[[4]]+tblLLStat[[5]]+tblLLStat[[10]]+tblLLStat[[11]]
nSuccessLL <- tblLLStat[[3]]+tblLLStat[[12]]


gc <- "NL"
strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12),labels=vHuntEventLabels )##Set To NoTHuntMode
message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
tblNLStat <- table(datHuntEvent$huntScore)

nFailNL <- tblNLStat[[4]]+tblNLStat[[5]]+tblNLStat[[10]]+tblNLStat[[11]]
nSuccessNL <- tblNLStat[[3]]+tblNLStat[[12]]


gc <- "DL"
strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12),labels=vHuntEventLabels )##Set To NoTHuntMode
message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
tblDLStat <- table(datHuntEvent$huntScore)

nFailDL <- tblDLStat[[4]]+tblDLStat[[5]]+tblDLStat[[10]]+tblDLStat[[11]]
nSuccessDL <- tblDLStat[[3]]+tblDLStat[[12]]
###


######## CALC Stat On Hunt Events ######
## Re-process Hunt Stat On Modified Events
lHuntStat <- list()
groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
strCondTags <- names(groupsrcdatList)
for (i in strCondTags)
{
  message(paste("#### ProcessGroup ",i," ###############"))
  strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
  message(paste(" Loading Hunt Events: ",strDataFileName))
  ##ExPORT 
  load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
  
  datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12),labels=vHuntEventLabels )##Set To NoTHuntMode
  ##Filter Hunt Events ##
  datHuntEventFilt <- datHuntEvent[datHuntEvent$huntScore != which(levels(huntLabels) == "NA") &
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Not_HuntMode/Delete") &
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Out_Of_Range") & 
                                   datHuntEvent$huntScore != which(levels(huntLabels) == "Duplicate/Overlapping") |
                                   datHuntEvent$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events
  
  
  lHuntStat[[i]] <- calcHuntStat3(datHuntEventFilt)
}

datHuntStat = do.call(rbind,lHuntStat)#
################