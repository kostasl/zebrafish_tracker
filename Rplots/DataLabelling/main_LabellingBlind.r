## USER VALIDATION OF HUNTEVENTS ##
# Jun 2022 : New version: Exports CSV of hunt events detected per experiment.
#            Tracker loads these in table and lets user view, score and modify hunt event data
# Kostasl Nov 2017
# 10-12-17 : Fixed Error on counting number of hunting episodes
#            Converted to Adding A Row for empty data files such that larva That were tested but produced no data count towards the mean / The case that a fish is invalid should be actually handled when testing in empty conditions and reject a fish that does nothing
#            Otherwise, a non appearing fish counts towards the group mean since its tested for a fixed amount of time (10mins each fish)
# 14-12-17 :  Added MotionTrajectory Analysis - PathLengths/Speed/Ratio #frames(Speed>0) over All ANalysed Event Frames Of A Larva 
# Note : Can use findLabelledEvent from auxFunctions.r file to obtain Labelled Event ID (ex.  "505/2149") from Using Event Register 
# Consider What the Hunt Ratio Is On a No Show Larva? Currently Set To 0 - 
#


# ## Pio Eykolo Na diaspasoume to Atomo Para mia prokatalipsi ##

library(tools)
library("MASS");
#library(hexbin)

##FIX 4443 to fail with strike
####################
#setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
source("config_lib.R")

setEnvFileLocations("LAB") #HOME,OFFICE,#LAPTOP


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
G_APPROXFPS              <- 60#410
G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 45 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
G_MIN_BOUTSPEED          <- 0.05 ##px/frame - Need to be above to be considered A Motion Bout
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments
nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags
# 
# 
# strDataSetDirectories <- paste(strTrackInputPath, list(
#                               "/Tracked12-10-17/", ##Dataset 1
#                               "/Tracked26-10-17/",
#                               "/Tracked02-11-17/",##MDataset 3 -NOTE: Does not Larva ID on File Name 
#                               "Tracked08-11-17/", #4 350fps - Missing a condition WTDryFed3Roti - So removed One Set Larva of Data from other conditions to balance the dataset
#                               "/Tracked16-11-17/",#5 400fps - Strict Timer Dataset
#                               "/Tracked30-11-17/",#6 420fps
#                               "/Tracked07-12-17/",#7
#                               "/Tracked14-12-17/",#8
#                               "Tracked21-12-17/",
#                               "/Tracked11-01-18/",
#                               "/Tracked18-01-18/",
#                               "/Tracked25-01-18/",
#                               "/Tracked01-02-18/",
#                               "/Tracked08-02-18/",
#                               "/Tracked15-02-18/",
#                               "/Tracked22-02-18/",
#                               "/Tracked_07-06-18/",##Dataset 17 
#                               "/Tracked14-06-18/",##Dataset 18
#                               "/Tracked_21-08-18/"##Dataset n ##Dataset n 
#                               ),sep="/")

strDataSetDirectories <- paste(strTrackInputPath, list(
  "DS_7dpf/" ##Dataset 2
),sep="/")



##Add Source Directory
strCondR  <- "*.csv"; 

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)

#source("HuntingEventAnalysis.r")
#source("TrajectoryAnalysis.r")
source("DataLabelling/labelHuntEvents_lib.r")

##For Safe Sampling Of Vectors Of Size 1
resample <- function(x, ...) x[sample.int(length(x), ...)]


firstDataSet = NROW(strDataSetDirectories)#-13
lastDataSet = NROW(strDataSetDirectories)
dataSetsToProcess = seq(from=firstDataSet,to=lastDataSet)


#message(paste(" Loading Hunt Event List to Validate... "))
#strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","KL","ALL",sep="-") ##To Which To Save After Loading
#load(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#datHuntEventAllGroupToValidate <- datHuntEventAllGroup

##Load The List To process

#strProcDataFileName <-paste("setn-12","-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
#strProcDataFileName <- paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
#strProcDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged3") ##To Which To Save After Loading
strProcDataFileName <- paste("HB_allHuntEvents.rds") ##To Which To Save After Loading
message(paste(" Loading Hunt Event List to Process... ",strProcDataFileName))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#datHuntEventAllGroupToLabel <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
##Load the latest set

datHuntEventAllGroupToLabel  <- getLabelledHuntEventsSet()

#datHuntEventAllGroupToLabel <- datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$expID == 6,]
 ##<- datHuntEvent
groupsList <- unique(datHuntEventAllGroupToLabel$groupID)
str_FilterLabel <- "UnLabelled"

## Old Method - Individual Events Scored ##
## scoreIndividualEventsRandomly(datHuntEventAllGroupToLabel,str_FilterLabel <- "UnLabelled")

#str_FilterLabel <- "NA"

## New method for validation using updated Tracker -
# VALIDATION Export Hunt Events table for each experiment - These are then imported to tracker for validation
## Write in separate folder so as not to overwrite ongoing validation files ##
out_Hdir <- paste0(strDataExportDir,"/LabelledSet/ToValidate_blank/")
Keyc = 'n'
if (!dir.exists(out_Hdir))
{
  dir.create(out_Hdir)
} ## If validation CSV subfolder does not Exist - make dir and export CSV with hunt events for each video
  for (expID in unique(datHuntEventAllGroupToLabel$expID) )
  {
    
      for (testCod in unique(datHuntEventAllGroupToLabel$testCond))
      {
        if (Keyc == 'q')
          break
        datHuntEvents_exp <- datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$expID == expID &
                                                         datHuntEventAllGroupToLabel$testCond == testCod,]
        
        for (eventID in unique(datHuntEvents_exp$eventID) )
        {
          ##if (eventID == 0) ## Skip Default event 0 which records rotifers
          #  next
          datHuntEvents_exp <-  datHuntEvents_exp[datHuntEvents_exp$eventID == eventID,]
          
          datEventsForTracker <- cbind.data.frame(rowID=(rownames(datHuntEvents_exp)),
                                              startFrame=datHuntEvents_exp$startFrame,
                                              endFrame=datHuntEvents_exp$endFrame,
                                              label=datHuntEvents_exp$huntScore)

      if (NROW(datEventsForTracker) > 0 )
      {
        filename_csv <- paste0(strsplit(head(basename(datHuntEvents_exp$filename),1),split = ".",fixed=T)[[1]][1],"_huntevents.csv")  
        filename_csv <- paste0("HB",expID,"_",testCod,"_",eventID,"_huntevents.csv")  
        if (!file.exists(paste0(out_Hdir,filename_csv))) ##Do not Overwrite Existing - In case user has updated it 
          write.table(datEventsForTracker,paste0(out_Hdir,filename_csv),sep=",",row.names=F )
        else
          warning("Will not overwrite existing : ", filename_csv," Table may have been validated by user. ")
          next()
      }else
      {
        filename_csv <- paste0("HB",expID,"_",testCod,"_",eventID,"_huntevents.csv")  
        write.table(datEventsForTracker,paste0(out_Hdir,filename_csv),sep=",",row.names=F )
        warning("No hunt events for expID:",expID," cond:",testCod)
        #break ##No Events for this one
      }
      
       ## Run Tracker passing tbl of Hunt Events
       message(paste("\n Validating video  ExpID:",expID,"event :",eventID," ",testCod ) )
       strVideoFile <- head(datHuntEvents_exp$filename,1)
       strArgs = paste0(" --HideDataSource=1 --MeasureMode=1 --ModelBG=1 --SkipTracked=0 --PolygonROI=0 --invideofile=",strVideoFile,
                       " --outputdir=",strTrackeroutPath," --DNNModelFile=","/home/meyerlab/workspace/zebrafishtrack/tensorDNN/savedmodels/fishNet_loc",
                       " --HuntEventsFile=",paste0(out_Hdir,filename_csv)," --startpaused=1")
      
       message(paste(strTrackerPath,"/zebraprey_track",strArgs,sep=""))
      
      if (!file.exists(paste(strTrackerPath,"/zebraprey_track",sep="")) )
        stop(paste("Tracker software not found in :",strTrackerPath ))
      
      execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout="",stderr =NULL) ## stdout=FALSE stderr = FALSE
      
      ## execres contains all of the stdout - so cant be used for exit code
      if (execres != 0)
        stop(execres) ##Stop If Application Exit Status is not success
      ##Show Labels And As
      
      
      Keyc <- readline(prompt="### Press q to exit, 'n' for next, or type event number you wish to label  :")
     
      
      }##For Each Test Condition
    } ##For Each Event
  }##For Each Experiment




# lLabelSummary <- list()
# nLabelledDL <- sum(tblRes[3:13,"DL"])
# nLabelledLL <- sum(tblRes[3:13,"LL"])
# nLabelledNL <- sum(tblRes[3:13,"NL"])
# 
# 
# message(paste("HuntEvents Labelled (exclude NA) #DL:",nLabelledDL,"#LL:",nLabelledLL,"#NL:",nLabelledNL ) )
# lLabelSummary$HuntEventCount <- list(DL=nLabelledDL,LL=nLabelledLL,NL=nLabelledNL)
# lLabelSummary$Success <- list(DL=sum(tblRes[c(3,12),"DL"]),LL=sum(tblRes[c(3,12),"LL"]),NL=sum(tblRes[c(3,12),"NL"]) )
# lLabelSummary$SuccessRatio <- list(DL=lLabelSummary$Success$DL/lLabelSummary$HuntEventCount$DL,LL=lLabelSummary$Success$LL/lLabelSummary$HuntEventCount$LL,NL=lLabelSummary$Success$NL/lLabelSummary$HuntEventCount$NL )

####
##Testing on 4581/542 , 446/1529
##Find Event
#expID <- "3521"
#eventID <- 14 
#datRes <- datHuntEventAllGroupToValidate[as.character(datHuntEventAllGroupToValidate$expID) == expID ,]
#datRes[as.character(datRes$eventID) == eventID ,]

######## CALC Stat On Hunt Events ######
## Re-process Hunt Stat On Modified Events
#source("HuntingEventAnalysis.r")
################

## Events Of Super Hunters 
message("Super Hunters List:")
datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$expID == 3541 & datHuntEventAllGroupToLabel$huntScore == 2,]