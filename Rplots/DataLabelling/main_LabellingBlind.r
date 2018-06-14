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
library("MASS");
#library(data.table) ##Required for rBindList
#library(hexbin)


####################




## GLOBAL VARS ###


## Office PC
setwd("/home/sabina/zfishLabel/")
strVideoFilePath  <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/expDataKostas/AnalysisSetAlpha/" 
#strVideoFilePath  <- "/media/sabina/zfishDataAlpha/AnalysisSetAlpha/" 
strTrackerPath    <- "/home/sabina/zfishLabel/build-zebraprey_track-Desktop_Qt_5_9_2_GCC_64bit-Release/"
strTrackeroutPath <- "/home/sabina/zfishLabel/out/"
strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from 
strDatDir        <- "./dat/" ##Where Are the Imported RData Stored
strDataExportDir <- "./out/"



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

strCondR  <- "*.csv"; 


### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)

#source("HuntingEventAnalysis.r")
#source("TrajectoryAnalysis.r")
source("labelHuntEvents.r")



firstDataSet = NROW(strDataSetDirectories)-11
lastDataSet = NROW(strDataSetDirectories)
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
# 
# ###COMBINE LISTS ##
# lHuntEvents <- list()
# gc <- "LL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# 
# lHuntEvents[[gc]] <- datHuntEvent
# 
# gc <- "NL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# datHuntEvent$markTracked <-NA
# lHuntEvents[[gc]] <- datHuntEvent
# 
# gc <- "DL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# datHuntEvent$markTracked <-NA
# lHuntEvents[[gc]] <- datHuntEvent
# 
# datHuntEventAllGroup = do.call(rbind,lHuntEvents)
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","SB","ALL",sep="-") ##To Which To Save After Loading
# save(datHuntEventAllGroup,file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# 

message(paste(" Loading Hunt Event List to Validate... "))
strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","KL","ALL",sep="-") ##To Which To Save After Loading
load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEventAllGroupToValidate <- datHuntEventAllGroup

groupsList <- unique(datHuntEventAllGroupToValidate$groupID)

##Load The List To process
strProcDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","SB","ALL",sep="-") ##To Which To Save After Loading
message(paste(" Loading Hunt Event List to Process... "))
load(file=paste(strDatDir,"/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#datHuntEventAllGroupToLabel <- datHuntEventAllGroup


##Select Randomly From THe Already Labelled Set ##
##Main Sample Loop
Keyc <- 'n'
while (Keyc != 'q')
{
  Keyc <- readline(prompt="### Press q to exit, 'n' for next, or type event number you wish to label  :")
  
  if (Keyc == 'q')
    break
  
  TargetLabel = which(vHuntEventLabels == "UnLabelled")-1;
  gc <- sample(groupsList,1)
  idx <- NA
  TargetLabels <- vHuntEventLabels
  
  if (Keyc == 'n')
  {
    datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore != "UnLabelled" & datHuntEventAllGroupToValidate$eventID != 0 ,]
    datHuntEventPool <- datHuntEventPool[datHuntEventPool$groupID == gc & datHuntEventPool$huntScore != TargetLabel,]
    expID <- sample(datHuntEventPool$expID,1)
    datHuntEventPool <- datHuntEventPool[datHuntEventPool$expID == expID ,]
    eventID <- sample(datHuntEventPool$eventID,1)
    ###
    TargetLabels <- vHuntEventLabels[vHuntEventLabels=="UnLabelled"]
  }
  
  if (!is.na(as.numeric(Keyc) ) )
  {
    idx <- as.numeric(Keyc)
    datHuntEventPool <- datHuntEventAllGroupToLabel[idx,]
    expID <- datHuntEventPool$expID
    eventID <- datHuntEventPool$eventID
    TargetLabels <- vHuntEventLabels
  }
  ##ExPORT 

  datHuntEventAllGroupToLabel <- labelHuntEvents(datHuntEventAllGroupToLabel,
                                                 strProcDataFileName,strVideoFilePath,
                                                 strTrackerPath,strTrackeroutPath,
                                                 TargetLabels,expID,eventID,idx)
  
  ##Saving is done in labelHuntEvent on Every loop - But repeated here
  save(datHuntEventAllGroupToLabel,file=paste(strDatDir,"/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
  strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","SB","ALL",sep="-") ##To Which To Save After Loading
  save(datHuntEventAllGroupToLabel,file=paste(strDatDir,"/",strProcDataFileName,"-backup.RData",sep="" )) ##Save With Dataset Idx Identifier
  message(paste("Saved Backup :",strDatDir,"/",strProcDataFileName,"-SB-backup.RData",sep="") )
  
  
}

##########################
####

# 
# 
# ########################################## SUMMARY OF LABELLING #####################
# ##How to Summarize Success / Fail Scores :
# gc <- "LL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# datHuntEvent$huntScore <- convertToScoreLabel( datHuntEvent$huntScore)##Set To NoTHuntMode
# message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
# tblLLStat <- table(datHuntEvent$huntScore)
# write.csv(tblLLStat,file="tbLLHuntLabelStat.csv")
# 
# nFailLL <- tblLLStat[[4]]+tblLLStat[[5]]+tblLLStat[[10]]+tblLLStat[[11]]
# nSuccessLL <- tblLLStat[[3]]+tblLLStat[[12]]
# 
# 
# gc <- "NL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# datHuntEvent$huntScore <- convertToScoreLabel( datHuntEvent$huntScore)##Set To NoTHuntMode
# message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
# tblNLStat <- table(datHuntEvent$huntScore)
# write.csv(tblNLStat,file="tbNLHuntLabelStat.csv")
# 
# nFailNL <- tblNLStat[[4]]+tblNLStat[[5]]+tblNLStat[[10]]+tblNLStat[[11]]
# nSuccessNL <- tblNLStat[[3]]+tblNLStat[[12]]
# 
# 
# gc <- "DL"
# strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents",gc,sep="-") ##To Which To Save After Loading
# load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
# datHuntEvent$huntScore <- convertToScoreLabel( datHuntEvent$huntScore)##Set To NoTHuntMode
# message(paste(NROW(datHuntEvent[datHuntEvent$huntScore != "UnLabelled",]),"/",NROW(datHuntEvent), " Data has already been labelled" ) )
# tblDLStat <- table(datHuntEvent$huntScore)
# write.csv(tblDLStat,file="tbDLHuntLabelStat.csv")
# 
# nFailDL <- tblDLStat[[4]]+tblDLStat[[5]]+tblDLStat[[10]]+tblDLStat[[11]]
# nSuccessDL <- tblDLStat[[3]]+tblDLStat[[12]]
# 
# message(paste("Rates:",nSuccessLL/nFailLL,nSuccessNL/nFailNL,nSuccessDL/nFailDL,sep="  "))
# ###
# 
# 
# ######## CALC Stat On Hunt Events ######
# ## Re-process Hunt Stat On Modified Events
# source("HuntingEventAnalysis.r")
# lHuntStat <- list()
# groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
# strCondTags <- names(groupsrcdatList)
# for (i in strCondTags)
# {
#   message(paste("#### ProcessGroup ",i," ###############"))
#   strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
#   message(paste(" Loading Hunt Events: ",strDataFileName))
#   ##ExPORT 
#   load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#   
#   datHuntEvent$huntScore <- factor(x=datHuntEvent$huntScore,levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13),labels=vHuntEventLabels )##Set To NoTHuntMode
#   ##Filter Hunt Events ##
#   datHuntEventFilt <- datHuntEvent[datHuntEvent$huntScore != "NA" &
#                                    datHuntEvent$huntScore != "Not_HuntMode/Delete" &
#                                    datHuntEvent$huntScore != "Out_Of_Range" & 
#                                    datHuntEvent$huntScore != "Duplicate/Overlapping" &
#                                    datHuntEvent$huntScore != "Near-Hunt State" |
#                                    datHuntEvent$eventID   == 0 , ] ##Keep THose EventID 0 so as to identify All experiments - even those with no events
#   
#   
#   lHuntStat[[i]] <- calcHuntStat3(datHuntEventFilt)
# }
# 
# datHuntStat = do.call(rbind,lHuntStat)#
# ################