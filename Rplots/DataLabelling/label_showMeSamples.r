### I used showMeSamples.r to reTrack specific Events of given Outcome along with prey tracking- 
## The Tracked files are then imported and combined in runimportHuntEventAnalysisDataFiles.r
##


source("DataLabelling/labelHuntEvents.r")

message(paste(" Loading Hunt Event List to Validate... "))
#strDataFileName <- paste("setn14-D5-18-HuntEvents-Merged") ##To Which To Save After Loading
strDataFileName <-paste("setn12-HuntEvents-SB-Updated",sep="") ##To Which To Save After Loading
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

TargetLabel <- as.numeric(readline(prompt="### Key In A Number For Which Label You want to see:"))
#TargetLabel = which(vHuntEventLabels == vHuntEventLabels[Keyc])-1;
#gc <- sample(groupsList,1)
datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore == TargetLabel &
                                                     datHuntEventAllGroupToValidate$eventID != 0 &
                                                     is.na(datHuntEventAllGroupToValidate$markTracked),]
expID <- sample(datHuntEventPool$expID,1)
datHuntEventPool <- datHuntEventPool[datHuntEventPool$expID == expID ,]
eventID <- sample(datHuntEventPool$eventID,1)

##ExPORT 2
datHuntEventPool <- labelHuntEvents(datHuntEventAllGroupToValidate,
                                    strDataFileName,strVideoFilePath,
                                    strTrackerPath,strTrackeroutPath,
                                    convertToScoreLabel(TargetLabel),expID,eventID)


# 