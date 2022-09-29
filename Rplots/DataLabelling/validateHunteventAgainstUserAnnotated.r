### Match Hunt Event Annotation between manual and automated methods
## Kostas 2022
source("config_lib.R")
source("HuntingEventAnalysis_lib.r")
setEnvFileLocations("LAPTOP") #HOME,OFFICE,#LAPTOP


load("/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/Analysis/dat//setn1_Dataset_VAL.RData")

vExpID <- c(3)


expID <- vExpID[1]
## Load Manually Labelled Data for Exp
strFileUserHuntEvents <- paste0(strDataExportDir,"ManuallyLabelled/fish",expID,"_video_mpeg_fixed_huntEvents.csv") 
datHuntEventsM <- read.csv(
  file=strFileUserHuntEvents, header = T)

## Load Automated detection
datHuntEvents <- detectHuntEvents(datAllFrames,vExpID,"LR",1)

## Make Unique pairs for each element of vector of start frames
datAllStartFramePairs <- expand.grid(manual=datHuntEventsM$startFrame,automatic=datHuntEvents$startFrame)

datAllStartFramePairs$frameDistance <- abs(datAllStartFramePairs$manual-datAllStartFramePairs$automatic)
idxSort <- order(datAllStartFramePairs$frameDistance,decreasing = FALSE)
nTopSelected <- NROW(datHuntEventsM)
## Select Closest For Each Manual Event ##
datAllStartFramePairs_top <- datAllStartFramePairs[head(idxSort, nTopSelected ),]
## Select Events That 
#datAllStartFramePairs_matched <- datAllStartFramePairs[datAllStartFramePairs$frameDistance,]

## Find Closesto automatically detected frame to the Manually labelled one
vFrameDistToAutoDetectedHuntEvent <- tapply(datAllStartFramePairs_top$frameDistance,datAllStartFramePairs_top$manual,min)
## Get Number of Detected events - use thres between auto detected and manual event 
vTruePositiveDetected <- vFrameDistToAutoDetectedHuntEvent[vFrameDistToAutoDetectedHuntEvent < HUNTEVENT_MATCHING_OFFSET] 
nTruePositiveDetected <- NROW(vTruePositiveDetected)
## The remaining manually labelled events that were not matched are counted as falselly classified as negative
nFalseNegativeDetected <- NROW(vFrameDistToAutoDetectedHuntEvent) - NROW(vTruePositiveDetected)
vValidatedAutoDetectedEvents <- as.numeric(names(vTruePositiveDetected))
## How likely is it that tracker detects a hunt event 
sensitivity <- nTruePositiveDetected/(nTruePositiveDetected + nFalseNegativeDetected)
## How likely is it that it responds specific to genuine hunt events
nFalsePositives <- NROW(datHuntEvents) - NROW(datHuntEventsM)
## Since we are classifying each frame, then here All non-Hunt Frames classified as such are True negatives - problem is the majority of frames are  true negatives are hunt events are generally rare
# Sum Total Automatic Detected Hunt Frames
nDetectedHuntFrames <- sum(datHuntEvents$endFrame - datHuntEvents$startFrame)  
nTrueNegative <- 
## Plot Manual and Automatic

datExpFrames = datAllFrames[datAllFrames$expID == vExpID[1] ,]

plot(datExpFrames$frameN,(datExpFrames$LEyeAngle-datExpFrames$REyeAngle),type="l",ylim=c(0,70),ylab="Eye vergence",xlab="frame N",xlim=c(0,10000))
abline(h=G_THRESHUNTVERGENCEANGLE,lwd=2,lty=2)
points(datAllStartFramePairs_top$automatic,rep(60,NROW(datAllStartFramePairs_top)),pch=2,col="red")
#points(datAllStartFramePairs_top$manual,rep(64,NROW(datAllStartFramePairs_top)),pch=6,col="blue")
points(datHuntEventsM$startFrame,rep(63,NROW(datHuntEventsM)),pch=25,col="blue")
points(vValidatedAutoDetectedEvents,rep(64,NROW(vValidatedAutoDetectedEvents)),pch=25,col="purple")
legend("bottomright",
       legend = c(paste("auto n",NROW(datHuntEvents)),paste("manual n",NROW(datHuntEventsM)),paste("matched n",nTruePositiveDetected) ) ,col=c("red","blue","purple"),pch=c(2,25,25) )
title(paste("F", unique(datHuntEvents$expID) ,"Hunt event detection sensitivity:",prettyNum(sensitivity*100,digits=4)) )

