## Script Assisting to Measure Fish Length from each larval experiment - (in order to check correlation with success)
# Kostasl 2 May 2020

library(tools)
library("MASS");


source("config_lib.R")

#setEnvFileLocations("HOME") #HOME,OFFICE,#LAPTOP


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
G_APPROXFPS              <- 410
G_THRESHUNTANGLE         <- 19 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 45 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
G_MIN_BOUTSPEED          <- 0.05 ##px/frame - Need to be above to be considered A Motion Bout
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments
nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags

### Set Colour Pallette Size from List Of Datasets
G_DATASETPALLETSIZE = NROW(strDataSetDirectories)

#source("HuntingEventAnalysis.r")
#source("TrajectoryAnalysis.r")
source("DataLabelling/labelHuntEvents_lib.r")

##For Safe Sampling Of Vectors Of Size 1
resample <- function(x, ...) x[sample.int(length(x), ...)]

strProcDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged3") ##To Which To Save After Loading
message(paste(" Using Hunt Event List to Process... ",strProcDataFileName))

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
strGroupID <- c(levels(datTrackedEventsRegister$groupID),"DE","LE","NE")


datHuntEventAllGroupToLabel  <- getLabelledHuntEventsSet()
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)

## The scriptlet to run the labelling process on a set of expID is found in auxFunctions.r
datFlatPxLength <- readRDS(file= paste(strDataExportDir,"/FishLength_Updated3.rds",sep=""))
message(paste(" Loading Measured fish length in pixels data ... "))

 ##<- datHuntEvent
groupsList <- c("DL","NL","LL") ##unique(datHuntEventAllGroupToLabel$groupID)
str_FilterLabel <- c("UnLabelled","Success","Fail","Fail-No Strike")
#str_FilterLabel <- "NA"
##Select Randomly From THe Already Labelled Set ##
##Main Sample Loop
Keyc <- 'n'
while (Keyc != 'q')
{
  Keyc <- readline(prompt="### Press q to exit, 'n' for next, or type event number you wish to label  :")
  
  if (Keyc == 'q')
    break
  
  TargetLabel = which(vHuntEventLabels %in% str_FilterLabel)-1; ##Convert to Number Score
  gc <- resample(groupsList,1)
  idx <- NA
  TargetLabels <- vHuntEventLabels
  ##Select All non Measured ExpID
  vExpID_ToMeasure <- datFishSuccessRate[!(datFishSuccessRate$expID %in% datFlatPxLength$expID),"expID"]
  
  
  if (Keyc == 'n')
  {
    ##Choose From THe Set Of Videos Already Labelled From Another User (Kostasl) So as to Verify The Label # Sample Only From THose ExpID that have not been already verified
    #datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore != "UnLabelled" & datHuntEventAllGroupToValidate$eventID != 0
    #                                           & (datHuntEventAllGroupToValidate$expID %in% datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$huntScore == TargetLabel,]$expID ),]
    
    datHuntEventPool <- datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$eventID != 0 &
                                                      datHuntEventAllGroupToLabel$expID %in% vExpID_ToMeasure ,]
    datHuntEventPool <- datHuntEventPool[ datHuntEventPool$huntScore %in% TargetLabel ,] #& is.na(datHuntEventPool$markTracked)
    if (NROW(datHuntEventPool)  == 0)
    {
      message( paste("Finished with Hunt Events for labels ",TargetLabels[TargetLabel+1], ". Try Again") )
      groupsList <- groupsList[which(groupsList != gc)]
      next
    }
    ##Choose Random Exp/Event from Subset
    expID <- resample(datHuntEventPool$expID,1)
    datHuntEventPool <- datHuntEventPool[datHuntEventPool$expID == expID ,]
    eventID <- resample(datHuntEventPool$eventID,1)
    ###
    TargetLabels <- vHuntEventLabels[vHuntEventLabels %in% str_FilterLabel] ##Convert to Text Label Score to Use for Filtering OUt
  }
  ##Extract If Any Numbers In Input/ Then User Picked a specific Row
  if (!is.na(as.numeric(gsub("[^0-9]","",Keyc)) ) )
  {
    message(paste("Goto Event:",Keyc ) )
    idx <- as.character(Keyc) ##Note It acts as key only as string, numeric would just bring out the respective order idx record
    datHuntEventPool <- datHuntEventAllGroupToLabel[idx,]
    expID <- datHuntEventPool$expID
    eventID <- datHuntEventPool$eventID
    TargetLabels <- vHuntEventLabels
    
    if (is.na(datHuntEventAllGroupToLabel[idx,]$expID))
    {
      message("Event Not Found")
      next
    }
  }
  ##ExPORT 

  
  datHuntEventAllGroupToLabel <- labelHuntEvents(datHuntEventAllGroupToLabel,
                                                 strProcDataFileName,strVideoFilePath,
                                                 strTrackerPath,strTrackeroutPath,
                                                 TargetLabels,expID,eventID,idx)
  
  ## Get User Measurement For Length
  KeyDat <- readline(prompt="### Type Fish Length In Pixels  (c exit):")
  
  if (KeyDat == 'c')
  {
    message("Stopping ", KeyDat)
    break;
  }
  
  groupID <- which(strGroupID == unique(datHuntEventPool[datHuntEventPool$expID == expID,]$groupID)  )
  datRow <- data.frame(KeyDat,groupID,expID)
  names(datRow) <- names(datFlatPxLength)
  datFlatPxLength <- rbind(datFlatPxLength,datRow)
  
  saveRDS(datFlatPxLength,file= paste(strDataExportDir,"/FishLength_Updated3.rds",sep="") )
  message(paste(" Saved updated fish length to FishLength_Updated3.csv. "))


}
#3951 106.508
#3941 93.10
#3421 86.03
#3751 95.52
#3811 94.59
#3731 98.47
#3871 92.97
#284 83.8153
#3541  97.63
#4612 96.0417
#3931 91.54
#4211 95.27
#4001 91.41 91.44
#3413 91.23 91.70 94.37 91.09
#3961 89.2749 89.96 90.60
#3441 92.72 93.13
#4041 106.2 107.4
#3681 97.94 100.62 98.61 100.49
#4071 100.49 100.88
#4051 98.00 100.285 96.93 98
#3511  95.52 97.80 92.913 93.98 95.00
#3991 93.338 95.078 97 95.42
#3591 110.11 111.73
#3591 110.60
#4251 88.45 94.13 92.61 91.52
#4561 95.80 95.210 94.00
#252 88.60 91.06 89 92.1
#3921 95.33 91.08 93.34 95 
#3781 94.429 90.79 94.11 94
#301 100.28 97.61 101.17
#285 93.193 92.633 94.047