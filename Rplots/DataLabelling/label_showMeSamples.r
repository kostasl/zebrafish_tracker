

message(paste(" Loading Hunt Event List to Validate... "))
strDataFileName <- paste("setn",NROW(dataSetsToProcess),"HuntEvents","KL","ALL",sep="-") ##To Which To Save After Loading
load(file=paste(strDatDir,"/",strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntEventAllGroupToValidate <- datHuntEventAllGroup

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
datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore == TargetLabel & datHuntEventAllGroupToValidate$eventID != 0 ,]
expID <- sample(datHuntEventPool$expID,1)
datHuntEventPool <- datHuntEventPool[datHuntEventPool$expID == expID ,]
eventID <- sample(datHuntEventPool$eventID,1)

##ExPORT 
datHuntEventPool <- labelHuntEvents(datHuntEventPool,
                                    strProcDataFileName,strVideoFilePath,
                                    strTrackerPath,strTrackeroutPath,
                                    convertToScoreLabel(TargetLabel),expID,eventID)


# 