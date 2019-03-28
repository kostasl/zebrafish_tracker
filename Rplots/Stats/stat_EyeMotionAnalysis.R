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
idx <- 20
plot(lEyeMotionDat[[idx]][,"t"],lEyeMotionDat[[idx]][,"REyeAngle"],type="l")