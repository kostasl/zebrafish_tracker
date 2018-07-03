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
setwd("/home/klagogia1/workspace/zebrafishtrack/Rplots")
strVideoFilePath  <- "/media/extStore/ExpData/zebrapreyCap/AnalysisSet/"
strTrackerPath    <- "/home/klagogia1/workspace/build-zebraprey_track-Release/" 
strTrackeroutPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/"
strTrackInputPath <- "/media/extStore//kostasl/Dropbox/Calculations/zebrafishtrackerData/" ##Same As Working Dir
strDatDir        <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
strDataExportDir <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
strPlotExportPath <- "/media/extStore/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots" ##Where to source the Tracker csv files from 


## Office PC
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
strVideoFilePath  <- "/media/LinuxDat/expDataKostas/AnalysisSetAlpha/" 
strTrackerPath    <- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_11_1_GCC_64bit-Release/"
strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/"
#strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from
strTrackInputPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/" ##Where to source the Tracker csv files from
strDatDir        <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
strDataExportDir <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
strPlotExportPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots" ##Where to source the Tracker csv files from 

## Laptop
setwd("~/workspace/zebrafishtrack/Rplots")
strVideoFilePath  <- "/media/kostasl/FLASHDATA/AnalysisSet"
strTrackerPath <-  "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release"
strTrackeroutPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackerOnHuntEvents_UpTo22Feb/"
strTrackInputPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData"##Where to source the Tracker csv files from 
strDatDir        <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
strDataExportDir <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
strPlotExportPath <- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots"
####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis.r")

source("TrajectoryAnalysis.r")

source("labelHuntEvents.r")
########


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
G_APPROXFPS              <- 420
G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
G_MIN_BOUTSPEED          <- 0.05 ##px/frame - Need to be above to be considered A Motion Bout
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments

nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags


strDataSetDirectories <- paste(strTrackInputPath, list(
                              "/Tracked12-10-17/", ##Dataset 1
                              "/Tracked26-10-17/",
                              "/Tracked02-11-17/",##MDataset 3 -NOTE: Does not Larva ID on File Name 
                              "Tracked08-11-17/", #4 350fps - Missing a condition WTDryFed3Roti - So removed One Set Larva of Data from other conditions to balance the dataset
                              "/Tracked16-11-17/",#5 400fps - Strict Timer Dataset
                              "/Tracked30-11-17/",#6 420fps ## Most Simular - Start ANalysis From here
                              "/Tracked07-12-17/",#7
                              "/Tracked14-12-17/",#8
                              "Tracked21-12-17/", # 9
                              "/Tracked11-01-18/",#10
                              "/Tracked18-01-18/",#11
                              "/Tracked25-01-18/",#12
                              "/Tracked01-02-18/",#13
                              "/Tracked08-02-18/",#14
                              "/Tracked15-02-18/",#15
                              "/Tracked22-02-18/",#16
                              "/Tracked_07-06-18/",##Dataset 17 
                              "/Tracked14-06-18/"##Dataset n 
                              ),sep="/")
##Add Source Directory



strCondR  <- "*.csv"; 

#display.brewer.all() to see avaulable options
rfc <- colorRampPalette(rev(brewer.pal(11,'Set3')));
r <- c(rfc(30),"#FF0000");

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rfc(G_DATASETPALLETSIZE),"#FF00AA");

#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
##NOTE: Assumes Files Begin with "Auto" and end with "track"
  lastDataSet = NROW(strDataSetDirectories)
  firstDataSet = lastDataSet 
  source("runimportTrackerDataFiles.r") 

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  ##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
  firstDataSet = NROW(strDataSetDirectories)-1
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##oad Frames and HuntStats
  source("loadAllDataSets.r")

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


######## Plot Hunt Statistics using datHuntStat##
source("plotHuntStat.r") 

###
source("plotMotionStat.r")

  source("labelHuntEvents.r") ##Library 
  source("DataLabelling/main_LabellingBlind.r")


  
  ######## CALC Stat On Hunt Events ######
## Re-process Hunt Stat On Modified Events
source("HuntingEventAnalysis.r")
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
  
  datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13),labels=vHuntEventLabels )##Set To NoTHuntMode
  ##Filter Hunt Events ##
  datHuntEventFilt <- datHuntEvent[datHuntEvent$huntScore != "NA" &
                                   datHuntEvent$huntScore != "Not_HuntMode/Delete" &
                                   datHuntEvent$huntScore != "Out_Of_Range" & 
                                   datHuntEvent$huntScore != "Duplicate/Overlapping" &
                                   datHuntEvent$huntScore != "Near-Hunt State" |
                                   datHuntEvent$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events
  
  
  lHuntStat[[i]] <- calcHuntStat3(datHuntEventFilt)
}

datHuntStat = do.call(rbind,lHuntStat)#
################



# 
# #### LABEL MANUALLY THE HUNT EVENTS WITH THE HELP OF THE TRACKER ###
# gc <- "LL"
# firstDataSet = 5
# lastDataSet = 16 #NROW(strDataSetDirectories)
# dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
# 
# #strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"-D-",firstDataSet,"-",lastDataSet,"-","HuntEvents-",gc,sep="") ##To Which To Save After Loading
# message(paste(" Loading Hunt Events: ",strDataFileName))
# ##ExPORT 
# load(file=paste(strDatDir,"/HuntEvents/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# TargetLabel = which(vHuntEventLabels == "UnLabelled")-1;
# 
# datHuntEvent <- labelHuntEvents(datHuntEvent,strDataFileName,strVideoFilePath,strTrackerPath,strTrackeroutPath, TargetLabel)
# ##Saving is done in labelHuntEvent on Every loop - But repeated here
# save(datHuntEvent,file=paste(strDatDir,"/",strDataFileName,"-backup.RData",sep="" )) ##Save With Dataset Idx Identifier
# message(paste("Saved Backup :",strDatDir,"/",strDataFileName,"-backup.RData",sep="") )
# 
# ##########################
# ####

