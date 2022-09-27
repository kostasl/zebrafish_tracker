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

datAllStartFramePairs$frameDistance <- abs(datAllStartFramePairs$Var1-datAllStartFramePairs$Var2)
idxSort <- order(datAllStartFramePairs$frameDistance,decreasing = FALSE)
datAllStartFramePairs_top <- datAllStartFramePairs[head(idxSort,NROW(datHuntEventsM)),]

## Plot Manual and Automatic

datExpFrames = datAllFrames[datAllFrames$expID == vExpID[1] ,]

plot(datExpFrames$frameN,(datExpFrames$LEyeAngle-datExpFrames$REyeAngle),type="l",ylim=c(0,70),ylab="Eye vergence")
abline(h=G_THRESHUNTVERGENCEANGLE,lwd=2,lty=2)
points(datAllStartFramePairs_top$manual,rep(60,NROW(datAllStartFramePairs_top)),pch=8,col="red")
points(datAllStartFramePairs_top$automatic,rep(64,NROW(datAllStartFramePairs_top)),pch=6,col="blue")
legend("bottomright",legend = c("auto","manual"),col=c("red","blue"),pch=c(8,6) )
