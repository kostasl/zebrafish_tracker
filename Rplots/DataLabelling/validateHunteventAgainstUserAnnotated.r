### Match Hunt Event Annotation between manual and automated methods
## Kostas 2022
source("config_lib.R")
source("HuntingEventAnalysis_lib.r")
setEnvFileLocations("LAPTOP") #HOME,OFFICE,#LAPTOP


load("/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/Analysis/dat/datAllFrames_Ds-1-4.RData")

vExpID <- c(104)
## Load Manually Labelled Data for Exp
datHuntEventsM <- read.csv(file="/mnt/data/Dropbox/Calculations/zebrafishtrackerData/OliviaHuntEventValidate/ManuallyLabelled/fish104_video_mpeg_fixed_huntEvents.csv", header = T)

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
vTruePositiveDetected <- vFrameDistToAutoDetectedHuntEvent[vFrameDistToAutoDetectedHuntEvent < G_MINGAPBETWEENEPISODES] 
nTruePositiveDetected <- NROW(vTruePositiveDetected)
vValidatedAutoDetectedEvents <- as.numeric(names(vTruePositiveDetected))
sensitivity <- nTruePositiveDetected/(NROW(vFrameDistToAutoDetectedHuntEvent))
## Plot Manual and Automatic

datExpFrames = datAllFrames[datAllFrames$expID == vExpID[1] ,]

plot(datExpFrames$frameN,(datExpFrames$LEyeAngle-datExpFrames$REyeAngle),type="l",ylim=c(0,70),ylab="Eye vergence")
abline(h=G_THRESHUNTVERGENCEANGLE,lwd=2,lty=2)
points(datAllStartFramePairs_top$automatic,rep(60,NROW(datAllStartFramePairs_top)),pch=2,col="red")
#points(datAllStartFramePairs_top$manual,rep(64,NROW(datAllStartFramePairs_top)),pch=6,col="blue")
points(datHuntEventsM$startFrame,rep(63,NROW(datHuntEventsM)),pch=25,col="blue")
points(vValidatedAutoDetectedEvents,rep(64,NROW(vValidatedAutoDetectedEvents)),pch=25,col="purple")
legend("bottomright",legend = c("auto","manual","matched"),col=c("red","blue","purple"),pch=c(2,25,25) )

