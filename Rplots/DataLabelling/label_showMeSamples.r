### I used showMeSamples.r to reTrack specific Events of given Outcome along with prey tracking- 
## The Tracked files are then imported and combined in runimportHuntEventAnalysisDataFiles.r
##

source("DataLabelling/labelHuntEvents_lib.r")

message(paste(" Loading Hunt Event List to Validate... "))
#strDataFileName <- paste("setn14-D5-18-HuntEvents-Merged") ##To Which To Save After Loading
#strDataFileName <-paste("setn14-HuntEventsFixExpID-SB-Updated-Merged",sep="") ##To Which To Save After Loading
strDataFileName <-paste("setn15-HuntEvents-SB-Updated-Merged",sep="") ##To Which To Save After Loading
datHuntEventAllGroupToValidate <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
 
groupsList <- unique(datHuntEventAllGroupToValidate$groupID)

###  SHOW ME SAMPLES ## 
vHuntEventLabels
l=0
for (g in vHuntEventLabels )
{
  message(paste(l,g,sep="-"))
  l=l+1
}

TargetGroups <- "LL" #c("LL","DL","NL")
TargetLabel <- as.numeric(readline("### Key In A Number For Which Label You want to see:"))
#TargetLabel = which(vHuntEventLabels == vHuntEventLabels[Keyc])-1;
#gc <- sample(groupsList,1)
datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore == TargetLabel &
                                                     datHuntEventAllGroupToValidate$eventID != 0 & 
                                                     datHuntEventAllGroupToValidate$groupID %in% TargetGroups 
                                                   ,] ##&is.na(datHuntEventAllGroupToValidate$markTracked)

if (NROW(datHuntEventPool) == 0)
  stop("No Hunt records with that label were found")

expID <- sample(datHuntEventPool$expID,1)
datHuntEventPool <- datHuntEventPool[as.character(datHuntEventPool$expID) == expID ,]
#datHuntEventPool <- datHuntEventPool[as.character(datHuntEventPool$eventID) == eventID ,]
eventID <-resample(datHuntEventPool$eventID,1)

##ExPORT 2
datHuntEventPool <- labelHuntEvents(datHuntEventAllGroupToValidate,
                                    strDataFileName,strVideoFilePath,
                                    strTrackerPath,strTrackeroutPath,
                                    convertToScoreLabel(TargetLabel),expID,eventID,
                                    idxFilter=NA,
                                    bskipMarked = FALSE)

#579/2017
# 