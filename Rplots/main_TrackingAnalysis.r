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
library(extrafont) ##For Font Embdedding in PDF, Run import_fonts() after install
#library(data.table) ##Required for rBindList
#library(hexbin)
rm("temp","subsetDat","TrackerData","frameNAll");

setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
source("config_lib.R")

setEnvFileLocations("OFFICE") #HOME,OFFICE,#LAPTOP

source("HuntingEventAnalysis.r")
source("TrajectoryAnalysis.r")
source("DataLabelling/labelHuntEvents_lib.r")
########

strDataSetDirectories <- paste(strTrackInputPath, list(
                              "/Tracked12-10-17/", ##Dataset 1
                              "/Tracked26-10-17/",
                              "/Tracked02-11-17/",##MDataset 3 -NOTE: Does not Larva ID on File Name 
                              "Tracked08-11-17/", #4 350fps - Missing a condition WTDryFed3Roti - So removed One Set Larva of Data from other conditions to balance the dataset
                              "/Tracked16-11-17/",#5 400fps - Strict Timer Dataset -  ANalysis Starts From here ***
                              "/Tracked30-11-17/",#6 420fps ## Most Simular 
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
                              "/Tracked14-06-18/",##Dataset 18
                              "/Tracked_21-08-18/"##Dataset 19
                              ),sep="/")
##Add Source Directory


### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)
rDataset <- c(rfc(G_DATASETPALLETSIZE),"#FF00AA");

strCondR  <- "*.csv"; 
#display.brewer.all() to see avaulable options


#################IMPORT TRACKER FILES # source Tracker Data Files############################### 
##Saves imported Data In Group Separeted RData Files as setn1_Dataset_...RData
##NOTE: Assumes Files Begin with "Auto" and end with "track"
  lastDataSet = NROW(strDataSetDirectories)
  firstDataSet = lastDataSet -14
  source("runimportTrackerDataFiles.r") 

###### END OF IMPORT TRACKER DATA ############


### LOAD Imported Data Sets - Starting From firstDataSet
  ##Alternatevelly Load The Complete Set From datAllFrames_Ds-5-16-.RData ##Avoids data.frame bug rbind
  firstDataSet = NROW(strDataSetDirectories)-14
  lastDataSet = NROW(strDataSetDirectories)
  dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)
  ##Load All Tracked Frames into datAllFrames data frame ##
  ##Warning : Merging is memory intensive
  #source("loadAllDataSets.r")
  ## Best to Load  a merged file datAllFrames instead :
  load(paste(strDatDir,"datAllFramesFix1_Ds-5-19.RData",sep="/"))
  load(paste(strDatDir,"groupsrcdatListPerDataSet_Ds-5-19.RData",sep="/"))
  load(file =paste(strDataExportDir,"/setn",NROW(dataSetsToProcess),"D",firstDataSet,"-",lastDataSet,"datTrackletStat.RData",sep=""))
  
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

  source("DataLabelling/labelHuntEvents_lib.r")
  source("DataLabelling/main_LabellingBlind.r")


  
  ######## CALC Stat On Hunt Events ######
## Re-process Hunt Stat On Modified Events
source("HuntingEventAnalysis_lib.r")
  ##Assumes datHuntEvent contains all groups
 datHuntStat <- makeHuntStat(datHuntEvent)  
#  
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

