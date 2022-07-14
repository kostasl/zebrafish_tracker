## Processing Of Tracker Eye TrackerData, Using R 
# Kostasl Nov 2017
# 10-12-17 : Fixed Error on counting number of hunting episodes
#            Converted to Adding A Row for empty data files such that larva That were tested but produced no data count towards the mean / The case that a fish is invalid should be actually handled when testing in empty conditions and reject a fish that does nothing
#            Otherwise, a non appearing fish counts towards the group mean since its tested for a fixed amount of time (10mins each fish)
# 14-12-17 :  Added MotionTrajectory Analysis - PathLengths/Speed/Ratio #frames(Speed>0) over All ANalysed Event Frames Of A Larva 
#            
## NOTES : The processing sequence of data:
##    *   Hunt Events are extracted from All tracking data using eye vergence theeshold
##    *   Hunting Events Bouts are analysed and saved in huntEpisodeAnalysis_MotionBoutDataXXX file - 
##    *   Then each hunt event is checked and Labelled According to outcome - The start-End frame adjusted 
###   *   Succesful hunt episodes are retracked in supervised manner we can analyse in finer detail -
##    *  The bouts from each retracked  hunt episode are extracted / and labelled according to seauence 
##    *   Motion bouts were then validated and noted the position of mouth and prey just prior to capture - Appended to huntEpisodeAnalysis_MotionBoutDataXXX_Validated
##    *   From there a new data list that compunes capture info with First turn behaviour at the onset of hunting is created :
##        Confusingly named it is called the firstboutpoints, with latest filename :  huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated
## 
### -- Clarifying the extracted data :
### 1st Pass Tracked Data is in : "datAllFramesFix1_Ds-5-19.RData"
##  2nd Pass Retracked (Successfull)  hunt events are in : datAllHuntEventAnalysisFrames_setC.RData, with a register datTrackedEventsRegister: setn_huntEventsTrackAnalysis_Register_ToValidate.rds
##  The Bout data After analysis of Retracked Events are in : huntEpisodeAnalysis_MotionBoutData_ToValidate.rds 
##  Data is then Labelled and register saved in :  "setn15-HuntEvents-SB-Updated-Merged3" loaded from getLabelledHuntEventsSet() function - LIst of labels : vHuntEventLabels
##  and once validated by manually adding position mouth and prey position these become: huntEpisodeAnalysis_MotionBoutData_peakSpeedExtended_Validated.rds - A manually recorded file of the labels also exists zebrafishtrackeData/HuntEvents_Retracked/boutValidate/manualRecs.txt
##  From there the Initial Turn and Capture Behaviour are linked in the lFirstBoutPoints : huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated.rds


#TODO: Add Colour Marker of Hunting On Trajectories
# ## Pio Eykolo Na diaspasoume to Atomo Para mia prokatalipsi ##
library(tools)
library(RColorBrewer);
library("MASS");
library(extrafont) ##For Font Embdedding in PDF, Run import_fonts() after install
library(here)
library(boot)
#library(data.table) ##Required for rBindList
#library(hexbin)
rm("temp","subsetDat","TrackerData","frameNAll");

setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
#setwd(here())
source("config_lib.R")

setEnvFileLocations("OFFICE") #HOME,OFFICE,#LAPTOP

source("HuntEpisodeAnalysis/HuntEpisodeAnalysis_lib.r")
source("TrajectoryAnalysis.r")
source("DataLabelling/labelHuntEvents_lib.r")

########   Directory to Source Tracker Exported CSV Files 
# Hunting Assay Experiments
# strDataSetDirectories <- paste(strTrackInputPath, list(
#                               "HB80_7dpf_LF3/", ##Dataset 2
#                               "HB70_7dpf_NF1/",
#                               "HB60_7dpf_LF2",
#                               "HB50_7dpf_NF0",
#                               "HB40_7dpf_LF1"
#                               ),sep="/")
## Hunger Exp
strDataSetDirectories <- paste(strTrackInputPath, list(
                                 "DS_7dpf/" ##Dataset 2
                                 ),sep="/")


### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rfc(G_DATASETPALLETSIZE),"#FF00AA");

strCondR  <- "*.csv"; 
#display.brewer.all() to see avaulable options


#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
##NOTE: Assumes Files Begin with "Auto" and end with "track"
## These need to be grouped in folders per GroupID
  lastDataSet   = NROW(strDataSetDirectories)
  firstDataSet  = 1 
  strFileNameFn = "extractFileNameParams_FOntogeny"
  source("runimportTrackerDataFiles.r") 

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  ##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
  firstDataSet = 1#NROW(strDataSetDirectories)#-14
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##Load All Tracked Frames into datAllFrames data frame ##
  ##Warning : Merging is memory intensive
  #source("loadAllDataSets.r")
  ## Best to Load  a merged file datAllFrames instead :
  #load(paste(strDatDir,"datAllFramesFix1_Ds-5-19.RData",sep="/"))
  #load(paste(strDatDir,"groupsrcdatListPerDataSet_Ds-5-19.RData",sep="/"))
  
  load(paste0(strDatDir,"/datAllFrames_Ds-",firstDataSet,"-",lastDataSet,".RData",sep=""))
  load(paste0(strDatDir,"/groupsrcdatListPerDataSet_Ds-",firstDataSet,"-",lastDataSet,".RData"))
  ## Load Tracklet Stat
  #load(file = paste(strDatDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datTrackletStat.RData",sep=""))
  
  ## \todo Problem Detecting Hunt Events -Across conditions
  ## Calculates HuntEvents And Hunt Statistics On Loaded Data ##
  groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]]
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  strCondTags <- names(groupsrcdatList)
  source("processLoadedData.r") ##Detects HuntEvents
 
  ##Once Processed you can Check and Validate Hunt Events Using main_LabellingBlind.r
  
  
    ### Make Eye Phase Space Density Plots ##########
  strCondTags <- unique(datAllFrames$groupID)
  for (i in strCondTags)
  {
    message(paste("#### Eye ProcessGroup ",i," ###############"))
    subsetDat = groupsrcdatList[[i]]
    strCond   <- paste(strCondR,subsetDat[2],collapse=NULL);
    
    ##Take All larva IDs recorded - Regardless of Data Produced - No Tracks Is Also Data
    #vexpID = unique(filtereddatAllFrames$expID)
    ##Select Larvaof this Group
    
    datAllGroupFrames <- datAllFrames[which(datAllFrames$groupID == i),]
    #Note:A Larva ID Corresponds to A specific Condition ex. NF1E (Same Fish Is tested in 2 conditions tho ex. NF1E, NF1L)
    vexpID = unique(datAllGroupFrames$expID)
    #plotGroupMotion(datAllGroupFrames,lHuntStat[[i]],vexpID)
    #######################################################################
    ###  EYE - PLOT Scatter and Eye Densities #####
    strCond = i;
    source("EyeScatterAndDensities.r")
    #####
  }


######## Plot Hunt Statistics using datHuntStat###
strDataFileName <-"setn1D1-1_datHuntStat.RData" #strHuntStateFilename# paste("setn14-D5-18-HuntEvents-Merged",sep="" )
source("plotHuntStat.r") 

###
source("plotMotionStat.r")

  source("DataLabelling/labelHuntEvents_lib.r")
  source("DataLabelling/main_LabellingBlind.r")


  
  ######## CALC Stat On Hunt Events ######
## Re-process Hunt Stat On Modified Events
source("HuntingEventAnalysis_lib.r")
  ##Assumes datHuntEvent contains all groups
 datHuntStat <- makeHuntStat(datHuntEvent)  

lHuntStat <- list()
groupsrcdatList <- groupsrcdatListPerDataSet[[NROW(groupsrcdatListPerDataSet)]] ##Load the groupsrcdatListPerDataSetFile
strCondTags <- names(groupsrcdatList)
for (i in strCondTags)
{
  message(paste("#### ProcessGroup ",i," ###############"))
  strDataFileName <- paste("/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
  message(paste(" Loading Hunt Events: ",strDataFileName))
  ##ExPORT
  load(file=paste(strDatDir,strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier

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
###############



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

