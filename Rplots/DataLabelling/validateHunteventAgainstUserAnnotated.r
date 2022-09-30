### Match Hunt Event Annotation between manual and automated methods
## Kostas 2022
source("config_lib.R")
source("HuntingEventAnalysis_lib.r")
setEnvFileLocations("LAPTOP") #HOME,OFFICE,#LAPTOP


load("/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/Analysis/dat//setn1_Dataset_VAL.RData")

vExpID <- unique(datAllFrames$expID)

for (expID in vExpID)
{
  ## Load Manually Labelled Data for Exp
  strFileUserHuntEvents <- paste0(strDataExportDir,"ManuallyLabelled/fish",expID,"_video_mpeg_fixed_huntEvents.csv") 
  datHuntEventsM <- read.csv(
    file=strFileUserHuntEvents, header = T)
  
  ## Load Automated detection
  datHuntEvents <- detectHuntEvents(datAllFrames,expID,"LR",1)
  
  datExpFrames = datAllFrames[datAllFrames$expID == expID ,]
  datExpEyeV = (datExpFrames$LEyeAngle-datExpFrames$REyeAngle)
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
  ## The remaining manually labelled events that were not matched are counted as falsely classified as negative
  nFalseNegativeDetected <- NROW(datHuntEventsM) - NROW(vTruePositiveDetected)
  vValidatedAutoDetectedEvents <- as.numeric(names(vTruePositiveDetected))
  ## How likely is it that tracker detects a hunt event 
  sensitivity <- nTruePositiveDetected/(nTruePositiveDetected + nFalseNegativeDetected)
  ## How likely is it that it responds specific to genuine hunt events
  nFalsePositives <- NROW(datHuntEvents) - NROW(datHuntEventsM)
  ## Since we are classifying each frame, then here All non-Hunt Frames classified as such are True negatives - problem is the majority of frames are  true negatives are hunt events are generally rare
  # Sum Total Automatic Detected Hunt Frames
  vValidatedAutoHuntEvents <- datAllStartFramePairs_top[datAllStartFramePairs_top$manual %in% as.numeric(names(vTruePositiveDetected)),"automatic"]
  # Count Number of Automatically Correctly Classified Frames
  nTruePositiveHuntFrames <- sum(datHuntEvents[datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$endFrame - 
                               datHuntEvents[datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$startFrame )
  nFalsePositiveFrames <- sum(datHuntEvents[!datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$endFrame - 
                               datHuntEvents[!datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$startFrame )
  # Count Total Frames With EyesVerged>THRESHOLD which have been correctly classified as non-Hunting 
  nTrueNegative <- NROW(datExpEyeV[datExpEyeV > G_THRESHUNTVERGENCEANGLE])-nTruePositiveHuntFrames
  ##  Specificity 
  specificity <- nTrueNegative/(nTrueNegative+nFalsePositiveFrames)
  ## Plot Manual and Automatic
  
  plot(datExpFrames$frameN,datExpEyeV,type="l",ylim=c(0,70),ylab="Eye vergence",xlab="frame N")
  abline(h=G_THRESHUNTVERGENCEANGLE,lwd=2,lty=2)
  points(datAllStartFramePairs_top$automatic,rep(60,NROW(datAllStartFramePairs_top)),pch=2,col="red")
  #points(datAllStartFramePairs_top$manual,rep(64,NROW(datAllStartFramePairs_top)),pch=6,col="blue")
  points(datHuntEventsM$startFrame,rep(63,NROW(datHuntEventsM)),pch=25,col="blue")
  points(vValidatedAutoDetectedEvents,rep(64,NROW(vValidatedAutoDetectedEvents)),pch=25,col="purple")
  legend("bottomright",
         legend = c(paste("auto n",NROW(datHuntEvents)),paste("manual n",NROW(datHuntEventsM)),paste("matched n",nTruePositiveDetected) ) ,col=c("red","blue","purple"),pch=c(2,25,25) )
  title(paste("F", expID ,"Hunt event sensitivity:",prettyNum(sensitivity*100,digits=4)," specificity:",prettyNum(specificity*100,digits=4) ) )
} ## each experiment
