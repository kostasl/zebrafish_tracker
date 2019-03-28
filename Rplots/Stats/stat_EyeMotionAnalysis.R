### Stat Analysis of Eye Movement Data ##


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


#### Load the Tracked Hunts Register ###
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #Processed Registry on which we add )


dev.off()

idx <- 10
plot(lEyeMotionDat[[idx]][,"t"],lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"],type="l")
plot(lEyeMotionDat[[idx]][,"DistToPrey"])
plot(lEyeMotionDat[[idx]][,"t"],lEyeMotionDat[[idx]][,"DistToPrey"])
plot(lEyeMotionDat[[idx]][,"DistToPrey"],lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"] )
## Does max EyeV Correlate with Hunt Duration Or With Distance To Prey Travelled ?
### Create A DAta Struct with EyeV before Strike/Capture - Duration of hunt Event-Eye Vergence - bUt also Distance from Prey travelled

##Get the Eye Vergence before the strike
lMaxEyeV <- list()
for (idx in 1:NROW(lEyeMotionDat) )
{
  vEye <- lEyeMotionDat[[idx]][,"LEyeAngle"]-lEyeMotionDat[[idx]][,"REyeAngle"]
  idx_Start <- min(which ( vEye > G_THRESHUNTVERGENCEANGLE))
  ##Closest To Prey in far time 
  idx_End <- max(which(lEyeMotionDat[[idx]][,"DistToPrey"] < median(lEyeMotionDat[[idx]][,"DistToPrey"]))) ##max(which ( vEye > G_THRESHUNTVERGENCEANGLE))
  maxV <- max(vEye[idx_End-100:idx_End ])
  
  huntDuration_ms <- lEyeMotionDat[[idx]][,"t"][idx_End]-lEyeMotionDat[[idx]][,"t"][idx_Start]
  huntDistance_mm <- lEyeMotionDat[[idx]][,"DistToPreyInit"][idx_End]-lEyeMotionDat[[idx]][,"DistToPrey"][idx_End]
  
  lMaxEyeV[[idx]] <- list(maxEyeV = maxV,duration_ms= huntDuration_ms,distance_mm=huntDistance_mm,RegIdx=lEyeMotionDat[[idx]][,"RegistarIdx"][idx_End] )
}


datMaxEyeV <- data.frame(do.call(rbind,lMaxEyeV))

cor(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$duration_ms) )
cor(unlist(datMaxEyeV$maxEyeV), unlist(datMaxEyeV$distance_mm) ) 


##Does Duration Explain Eye Vergence ? ##
plot(datMaxEyeV$duration_ms,datMaxEyeV$maxEyeV,ylim=c(0,100))
abline(lm(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$duration_ms)), col="red") # regression line (y~x) 
lines(lowess(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$duration_ms) ))

##Distance??
plot(datMaxEyeV$distance_mm,datMaxEyeV$maxEyeV,ylim=c(0,100))
abline(lm(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$distance_mm)), col="red") # regression line (y~x) 
lines(lowess(unlist(datMaxEyeV$maxEyeV)~unlist(datMaxEyeV$distance_mm) ))
