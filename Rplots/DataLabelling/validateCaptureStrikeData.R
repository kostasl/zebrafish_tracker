###  Kostasl 2018 
##  Validate the capture strike data that was automatically detected from the retracked huntevents :
##  * after running runHuntEpisodeAnalysis.r, a bout detection list for each of the hunt episodes is used to obtain what the undershoot 
##  turn was on firstbout, and link it to the final, capture bout intensity.


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") ) ## THis is the Processed Register File On 
## The Original list if the lFirstBout data from runHuntepisode analysis
datCaptureBoutsToValidate <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_CaptureBoutData_ToValidate.rds",sep="")) 

#datHuntEventAllGroupToValidate <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
##Save a backup before changing anything
saveRDS(datCaptureBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_CaptureBoutData_ToValidate_backup.rds",sep="")) ##Save With Dataset Idx Identifier

idxToValidate <- datCaptureBoutsToValidate$RegistarIdx

#idxToValidate <- c(idxReg_DL,idxReg_NL,idxReg_LL)

for (idx in idxToValidate)
{
  rec <- datTrackedEventsRegister[idx,]
  strFormatedEventID <- gsub(" ","0",formatC(format="fg",width=3,as.character(rec$eventID)))
  
  ##Get Respective Video  Filename / Path Uniquely identified by the  expID and event ID
  strVideoFile <- list.files(path =strVideoFilePath, pattern = paste0("_",rec$expID,"_",strFormatedEventID,".mp4") , all.files = FALSE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)  
  stopifnot(NROW(strVideoFile) == 1)
  print(strVideoFile)
}

saveRDS(datCaptureBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_CaptureBoutData_ToValidate.rds",sep="")) ##Save With Dataset Idx Identifier

