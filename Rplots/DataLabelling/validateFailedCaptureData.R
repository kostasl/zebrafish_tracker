##### Extract Capture Swim/Strike Data From 1st Tracking Attempt -
## I used this To complement my capture strike analysis of the successful events with those from failed hunt events 

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel

## Load Tracking Data datAllFrames
load(paste(strDatDir,"datAllFramesFix1_Ds-5-19.RData",sep="/"))

## Load Labelled Data set ###
datHuntEventAllGroupToLabel  <- getLabelledHuntEventsSet()
datHuntEventAllGroupToLabel$huntScore <- convertToScoreLabel(datHuntEventAllGroupToLabel$huntScore)


#vHuntEventLabels <- Lists all scores as factor
datFailedHuntEvents <- datHuntEventAllGroupToLabel[grepl("Fail-With Strike",as.character(datHuntEventAllGroupToLabel$huntScore) ) &
                                                  datHuntEventAllGroupToLabel$groupID %in% c('NL','DL')
                                                  ,]
##Save Empty Prototype Of where we will save capture bout data
#saveRDS(datCaptureBouts,file=paste(strDataExportDir,"/FailedHuntEpisodeAnalysis_CaptureBoutData.rds",sep="")) ##Save With Dataset Idx Identifier
datCaptureBouts <- readRDS(file=paste(strDataExportDir,"/FailedHuntEpisodeAnalysis_CaptureBoutData.rds",sep=""))

#Select The Hunt Events whose capture bouts we have not recorded and validate yet
idxToValidate <- rownames( datFailedHuntEvents)[!(rownames( datFailedHuntEvents) %in% datCaptureBouts[is.na(datCaptureBouts$MarkValidated)  ]$RegistarIdx)]


for (idx in idxToValidate)
{
  rec <- datFailedHuntEvents[idx,]
  strFormatedEventID <- gsub(" ","0",formatC(format="fg",width=3,as.character(rec$eventID)))
  
  ##Get Respective Video  Filename / Path Uniquely identified by the  expID and event ID
  strVideoFile <- list.files(path =strVideoFilePath, pattern = paste0("_",rec$expID,"_",strFormatedEventID,".mp4") , all.files = FALSE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)  
  if (NROW(strVideoFile) == 0)
    stop(paste("Could not find video file matching,",rec$filenames, " ,in : ",strVideoFilePath ) )
  if (!file.exists(strVideoFile) )
    stop(paste("Video File ",strVideoFile, "does not exist" ) )
  
  ##Retrieve Bout - If it exist - 
  recBout <- datCaptureBouts[datCaptureBoutsToValidate$RegistarIdx == idx, ]
  if (nrow(recBout) == 0)
  { newRow <- recBout[NA,] ##Make Empty Row
    newRow$boutRank         = 1 ###we are adding capture bout info
    newRow$vMotionBout_On  = 0 ##We do not know where the capture frame begins relevant to the startFrame of the event
    newRow$vMotionBout_Off = rec$endFrame-rec$startFrame
    newRow$expID           = rec$expID
    newRow$eventID         = rec$eventID
    newRow$groupID         = rec$groupID
    newRow$RegistarIdx     = idx  ##This is now different To the success as this ID refers to the Labelled hunt events reg, and not the retracked events reg
                         
    datCaptureBouts<-rbind(newRow,datCaptureBouts) ##add new row
    recBout <- newRow
  }
  
  message(paste("\n", row.names(rec) ,". Examining Hunt Event -start:",max(0,rec$startFrame + recBout$vMotionBout_On-1)," -End:",rec$startFrame + recBout$vMotionBout_Off, "ExpID:",rec$expID ) )
  strArgs = paste(" --MeasureMode=1 --HideDataSource=0 --ModelBG=0 --SkipTracked=0 --PolygonROI=1 --invideofile=",strVideoFile," --outputdir=",strTrackeroutPath,"/BoutValidate/",
                  " --startframe=",max(0,rec$startFrame + recBout$vMotionBout_On-1)," --stopframe=",max(0,rec$startFrame + recBout$vMotionBout_Off-1)," --startpaused=1",sep="")
  message(paste(strTrackerPath,"/zebraprey_track",strArgs,sep=""))
  if (!file.exists(paste(strTrackerPath,"/zebraprey_track",sep="")) )
    stop(paste("Tracker software not found in :",strTrackerPath ))
  
  execres <- base::system2(command=paste(strTrackerPath,"/zebraprey_track",sep=""),args =  strArgs,stdout=NULL,stderr =NULL) ## stdout=FALSE stderr = FALSE
  
  if (execres != 0)
    stop(execres) ##Stop If Application Exit Status is not success
  
  flush(con=stdin())
  flush(con=stdout())
  ##Show user Prompt for bout validation input
  message(paste("\n\n ### Event's ", row.names(rec) , "  ####" ) )
  
  strDat <-""
  arrDat <- ""
  while ( NROW(arrDat) < 5 && substr(strDat,1,1) != 'c' )
  {
    strDat <- readline(prompt=paste(idx, "# Validation Data [Prey X, Y, Mouth X,Y,Frame Bout On, Bout Off ] :") )
    strDat <- gsub(strDat,pattern = "[][]",replacement = "")
    arrDat <- strsplit(strDat,"," )[[1]]
  }
  
  
  ## extract validation data on prey location and fish mouth position ##
  strDat <- gsub(strDat,pattern = "[][]",replacement = "")
  arrDat <- strsplit(strDat,"," )[[1]]
  
  if (NROW(arrDat) < 5)
    break;
  
  stopifnot(NROW(arrDat) > 5 ) ##Stop if array data not found
  preyX <- as.numeric( arrDat[1]);  preyY  <- as.numeric(arrDat[2]);
  mouthX <- as.numeric(arrDat[3]); mouthY <- as.numeric(arrDat[4]);
  boutOn <- as.numeric(arrDat[5]) - rec$startFrame; boutOff <- as.numeric(arrDat[6]) - rec$startFrame
  fOnSetEyeVergence <- as.numeric(arrDat[7])
  ## \note Angle to prey Not Recalculated , changes to bout frames  invalidates vMotionBoutDistanceTravelled_mm
  recUpd<- within( datCaptureBouts[datCaptureBouts$RegistarIdx == idx, ],{ ##Multiple vars update
    vd_PreyY <- preyY 
    vd_PreyX <- preyX
    vd_MouthX <- mouthX
    vd_MouthY <- mouthY
    vMotionBout_On <- boutOn
    vMotionBout_Off <- boutOff 
    vMotionBoutDistanceToPrey_mm <- DIM_MMPERPX*(sqrt((preyX-mouthX)^2+(preyY-mouthY)^2 ))
    OnSetEyeVergence <- fOnSetEyeVergence
  })
  message(paste(recUpd$RegistarIdx,"Distance to prey:", prettyNum(digits=3,recUpd$vMotionBoutDistanceToPrey_mm)," eyeV:",prettyNum(digits=3,recUpd$OnSetEyeVergence) ))
  ## Pass data to record   // SAVE
  strKeyC <- readline(prompt="### Mark Validated ? (y/n):")
  if (strKeyC == 'y')
  {
    recUpd$MarkValidated <- 2 ## 2 Means Manually set Capture Bouts - Validated against initial tracking data/Not the retracked
    #datFailedHuntEvents[idx,]$markTracked <- NA
    ##update data frame
    datCaptureBouts[datCaptureBouts$RegistarIdx == idx, ] <- recUpd
    
    
    print(paste( recUpd$RegistarIdx, "# Validated") )
    stopifnot(!is.na(datCaptureBouts[datCaptureBouts$RegistarIdx == idx, ]$MarkValidated))
    
    ## print event Hunt Score label  / Allow to change 
    message("Outcome logged on Registry (matched) :", rec$huntScore )
    print(levels(convertToScoreLabel(rec$LabelledScore )))
    strKeyC <- readline(prompt=paste("### Change label to [",rec$huntScore ,"]:") )
    strKeyC <- as.numeric(strKeyC)
    if (is.numeric(strKeyC))
    {
      newLab <- convertToScoreLabel(as.numeric(strKeyC)-1 ) 
      message(paste("new label:",newLab ) )
      datHuntEventAllGroupToLabel[idx,] <- newLab        ####Save onto Origianl Full Record Using IDx
    }
    
    ##Save 
    saveRDS(datCaptureBouts,file=paste(strDataExportDir,"/FailedHuntEpisodeAnalysis_CaptureBoutData.rds",sep="")) ##Save With Dataset Idx Identifier
    saveRDS(datHuntEventAllGroupToLabel, file_LabelledHuntEventsSet ) ## THis is the Processed Register File On 
    message("HuntEvent and CaptureBout Files Updated.")
    
  } ##Mark validated
  
  
  strKeyC <- readline(prompt="### Move to next or quit ? (n/q):")
  if (strKeyC == 'q')
    break;
  
  
} ## end of loop 
# datMotionBoutsToValidate[datMotionBoutsToValidate$boutRank==1, ]  <- datCaptureBoutsToValidate
