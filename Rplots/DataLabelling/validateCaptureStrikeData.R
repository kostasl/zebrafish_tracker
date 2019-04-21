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
datMotionBoutsToValidate <-readRDS(file=paste0(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds") ) 

#datHuntEventAllGroupToValidate <-readRDS(file=paste(strDatDir,"/LabelledSet/",strDataFileName,".rds",sep="" )) ##Save With Dataset Idx Identifier
##Save a backup before changing anything
saveRDS(datMotionBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate_backup.rds",sep="")) ##Save With Dataset Idx Identifier

##Make desired columns if missing 
if (!any(names(datMotionBoutsToValidate) == "MarkValidated"))
  datMotionBoutsToValidate$MarkValidated <- NA


datCaptureBoutsToValidate <- datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ] 
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
  if (NROW(strVideoFile) == 0)
    stop(paste("Could not find video file matching,",rec$filenames, " ,in : ",strVideoFilePath ) )
  if (!file.exists(strVideoFile) )
    stop(paste("Video File ",strVideoFile, "does not exist" ) )
  
  recBout <- datCaptureBoutsToValidate[datCaptureBoutsToValidate$RegistarIdx == idx, ]
  
  message(paste("\n", row.names(rec) ,". Examining Hunt Event -start:",max(0,rec$startFrame + recBout$vMotionBout_On-1)," -End:",rec$endFrame + recBout$vMotionBout_Off, "ExpID:",rec$expID ) )
  strArgs = paste(" --MeasureMode=1 --HideDataSource=0 --ModelBG=0 --SkipTracked=0 --PolygonROI=1 --invideofile=",strVideoFile," --outputdir=",strTrackeroutPath,
                  " --startframe=",max(0,rec$startFrame + recBout$vMotionBout_On-1)," --stopframe=",max(0,rec$endFrame + recBout$vMotionBout_Off-1)," --startpaused=1",sep="")
  message(paste(strTrackerPath,"/zebraprey_track",strArgs,sep=""))
  if (!file.exists(paste(strTrackerPath,"/zebraprey_track",sep="")) )
    stop(paste("Tracker software not found in :",strTrackerPath ))
  
  execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout=NULL,stderr =NULL) ## stdout=FALSE stderr = FALSE
  
  
}

saveRDS(datCaptureBoutsToValidate,file=paste(strDataExportDir,"/huntEpisodeAnalysis_CaptureBoutData_ToValidate.rds",sep="")) ##Save With Dataset Idx Identifier

