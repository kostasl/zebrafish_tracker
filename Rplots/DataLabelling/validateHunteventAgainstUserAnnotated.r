## Match Hunt Event Annotation between manual and automated methods
# \brief:  Script  To Validate effectiveness of Tracker and data processing scripts to detect hunt events - compared to user defined
# I upgraded the tracker to allow user to denote hunt events startframes. Me and Olivia went through 12 videos (Validation Set) .
# The effectiveness of particular hunt event detection parameters in Config.R can be tested by running the script here each time params are modified.
# If zebrafishtracker is modified, then videos need to be reprocessed and csv files re-imported using main_TrackingAnalysis.r 
#
# / Kostas Lagogiannis 2022 /
source("config_lib.R")
source("HuntingEventAnalysis_lib.r")
setEnvFileLocations("HOME") #HOME,OFFICE,#LAPTOP
load(paste0(strDataStore,"/setn1_Dataset_TrackerValidation_081022.RData"))


G_THRESHUNTANGLE          <- 14 #Define Min Angle Both Eyes need to exceed for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE  <- 47 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_HUNTSCORETHRES          <- 0.96 ## Detection Thresh for DNN HUNT event detections (Based on fish image)
G_THRESHCLIPEYEDATA       <- 50 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES   <- G_APPROXFPS/3
G_MINEPISODEDURATION      <- G_APPROXFPS/3
HUNTEVENT_MATCHING_OFFSET <- 3*G_APPROXFPS # Max frames to accept as mismatch when matching manual to auto detected huntevents during tracker validation - frames to Used in validateHuntEventsAgainstUserAnnotated



vExpID <- unique(datAllFrames$expID)

lCompHuntEvents <- list()
vHuntScores <- c(0.1, round(100*seq(0.5,0.90,length=5))/100,0.99)
vEyeThres <- round(100*seq(40,50,length=16))/100

  for (expID in vExpID)
  {
    n=0;
    for (G_THRESHUNTVERGENCEANGLE in vEyeThres)
    {
      for (G_HUNTSCORETHRES in vHuntScores)
      {    
        message(expID)
        ## Load Manually Labelled Data for Exp
        strFileUserHuntEvents <- paste0(strDataExportDir,"/ManuallyLabelled/fish",expID,"_video_mpeg_fixed_huntEvents.csv") 
        if (!file.exists(strFileUserHuntEvents))
        {
          warning("MISSING hunt event file for expID:",expID,"-",strFileUserHuntEvents ,"*Skiped. ")
          next
        }
        datHuntEventsM <- read.csv(
          file=strFileUserHuntEvents, header = T)
        
        ## Load Automated detection
        datHuntEvents <- detectHuntEvents(datAllFrames,expID,"LR",1,
                                          THRESHUNTVERGENCEANGLE=G_THRESHUNTVERGENCEANGLE,
                                          HUNTSCORETHRES=G_HUNTSCORETHRES,
                                          THRESHCLIPEYEDATA = G_THRESHCLIPEYEDATA
                                          )
        datHuntEvents <- datHuntEvents[datHuntEvents$eventID > 0,]
        
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
        
        ## Find Closest to automatically detected frame to the Manually labelled one
        vFrameDistToAutoDetectedHuntEvent <- tapply(datAllStartFramePairs_top$frameDistance,datAllStartFramePairs_top$automatic,min)
        ## Get Number of Detected events - use thres between auto detected and manual event 
        vTruePositiveDetected <- vFrameDistToAutoDetectedHuntEvent[vFrameDistToAutoDetectedHuntEvent < HUNTEVENT_MATCHING_OFFSET] 
        nTruePositiveDetected <- NROW(vTruePositiveDetected)
        ## The remaining manually labelled events that were not matched are counted as falsely classified as negative
        ## COUNT Only uniquely matched manual to auto event - ie unique pairs
        nUniquelyMatchedManualEvents <- min( NROW(unique(datAllStartFramePairs_top[datAllStartFramePairs_top$frameDistance< HUNTEVENT_MATCHING_OFFSET,]$automatic)),
                                             NROW(unique(datAllStartFramePairs_top[datAllStartFramePairs_top$frameDistance< HUNTEVENT_MATCHING_OFFSET,]$manual)) )
        
        nFalseNegativeDetected <- NROW(datHuntEventsM) - nUniquelyMatchedManualEvents
        vValidatedAutoDetectedEvents <- as.numeric(names(vTruePositiveDetected))
        ## How likely is it that tracker detects a hunt event 
        #sensitivity <- nTruePositiveDetected/(nTruePositiveDetected + nFalseNegativeDetected)
        sensitivity<-nUniquelyMatchedManualEvents/(nUniquelyMatchedManualEvents+nFalseNegativeDetected)
        # FalsePositives : Remove number of validated events from total automatically detected events
        # Note: Each automatic event may be matched to multiple Manual Ones and so 
        nFalsePositives <- NROW(datHuntEvents) - nUniquelyMatchedManualEvents #NROW(vTruePositiveDetected)
        stopifnot(nFalsePositives >= 0)
        ## Since we are classifying each frame, then here All non-Hunt Frames classified as such are True negatives - problem is the majority of frames are  true negatives are hunt events are generally rare
        # Sum Total Automatic Detected Hunt Frames
        vValidatedAutoHuntEvents <- datAllStartFramePairs_top[datAllStartFramePairs_top$automatic %in% as.numeric(names(vTruePositiveDetected)),"automatic"]
        # Count Number of Automatically Correctly Classified Frames
        # by summing the duration of Validated automatically selected Hunt Events 
        nTruePositiveHuntFrames <- sum(datHuntEvents[datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$endFrame - 
                                     datHuntEvents[datHuntEvents$startFrame %in% vValidatedAutoHuntEvents,]$startFrame) 
        nFalsePositiveFrames <- sum(datHuntEvents[!(datHuntEvents$startFrame %in% vValidatedAutoHuntEvents),]$endFrame - 
                                     datHuntEvents[!(datHuntEvents$startFrame %in% vValidatedAutoHuntEvents),]$startFrame) 
        # Count Total Experiment Frames which have been correctly classified as non-Hunting 
        # Note: Given Imbalance in number of hunt frames to total frames specificity (fraction of -ve classified that are truly negative) will score very high
        # a subset of likely hunt frames need to be selected
        nTrueNegative <- NROW(datExpEyeV[datExpEyeV >= G_THRESHUNTVERGENCEANGLE*0.98])-nTruePositiveHuntFrames

        ##  Specificity 
        ## How likely is it that it responds specific to genuine hunt events
        stopifnot(nTrueNegative >= 0 & nFalsePositiveFrames >= 0 )
        specificity <- nTrueNegative/(nTrueNegative+nFalsePositiveFrames)
        
        # Positive Predicted Value
        PositivepredictiveValue <- nUniquelyMatchedManualEvents/(nUniquelyMatchedManualEvents + nFalsePositives)
        
        ## Plot Manual and Automatic
        strPlotName = paste(strPlotExportPath,"/fish",expID,"_",n,"-HEventMatching_HC",G_HUNTSCORETHRES,"_EyeV",G_THRESHUNTVERGENCEANGLE,".pdf",sep="")
        pdf(strPlotName)
          plot(datExpFrames$frameN,datExpEyeV,type="l",ylim=c(0,70),ylab="Eye vergence",xlab="frame N")
          abline(h=G_THRESHUNTVERGENCEANGLE,lwd=2,lty=2)
          points(datAllStartFramePairs$automatic,rep(60,NROW(datAllStartFramePairs)),pch='|',cex=1,col="red")
          points(datAllStartFramePairs$automatic,rep(60,NROW(datAllStartFramePairs))-0.5,pch='|',cex=0.3,col="black")
          #points(datAllStartFramePairs_top$manual,rep(64,NROW(datAllStartFramePairs_top)),pch=6,col="blue")
          points(datHuntEventsM$startFrame,rep(63,NROW(datHuntEventsM)),pch=25,col="blue")
          points(vValidatedAutoDetectedEvents,rep(64,NROW(vValidatedAutoDetectedEvents)),pch='|',col="purple")
          legend("bottomright",box.col = "white",bg="white",
                 legend = c(paste("auto n",NROW(datHuntEvents)),paste("manual n",NROW(datHuntEventsM)),
                            paste("matched n",nUniquelyMatchedManualEvents) ) ,
                           col=c("red","blue","purple"),pch=c(2,25,25) )
          title(paste("F", expID ,"H event sensitivity:",prettyNum(sensitivity*100,digits=4),
                      " specificity:",prettyNum(specificity*100,digits=4),
                      "PPV:",prettyNum(PositivepredictiveValue*100,digits=4)) )
  
          lCompHuntEvents[[paste0(expID,"_",n)]] <- data.frame(expID=expID,
                                                               ClassifierThres = G_HUNTSCORETHRES,
                                                               EyeVThres=G_THRESHUNTVERGENCEANGLE,
                                                               TruePositiveHuntFrames=nTruePositiveHuntFrames,
                                                               FalsePositiveFrames=nFalsePositiveFrames,
                                                               TrueNegativeFrames = nTrueNegative,
                                                               FalseNegativeFrames = nFalseNegativeDetected,
                                                               ManualCount=NROW(datHuntEventsM),
                                                               AutomaticCount=NROW(datHuntEvents),
                                                               Matched=nTruePositiveDetected,
                                                               Sensitivity=sensitivity,
                                                               Specificity=specificity,
                                                               PPV=PositivepredictiveValue,
                                                               rSq=NA)
      dev.off()        
        n = n + 1
      } ## Eye V
     } # VHuntScore  
  } ## each experiment
    
  ## Do Linear Fit Through 
  datCompEvents <- do.call(rbind,lCompHuntEvents)

  
  for (testedEyeV in vEyeThres)
  {
    for (testedHuntScore in vHuntScores)
    {    
      
      lmmodel <- lm(AutomaticCount~ManualCount,data=datParamComp)
      ## Extract R square from Linear fit
      rSq <- summary(lmmodel)$r.squared 
      datCompEvents[datCompEvents$EyeVThres == testedEyeV &
                      datCompEvents$ClassifierThres == testedHuntScore,"rSq"] = rSq
      datCompEvents[datCompEvents$EyeVThres == testedEyeV &
                      datCompEvents$ClassifierThres == testedHuntScore,"lm_slope"] = lmmodel$coefficients[2]
      datCompEvents[datCompEvents$EyeVThres == testedEyeV &
                      datCompEvents$ClassifierThres == testedHuntScore,"lm_intercept"] = lmmodel$coefficients[1]
      
      
      datParamComp <-  datCompEvents[datCompEvents$EyeVThres == testedEyeV &
                                       datCompEvents$ClassifierThres == testedHuntScore,  ]
      
      ## Success / Strike Non Strike Percentage ##
      strPlotName = paste(strPlotExportPath,"/trackerHEventVal_HC",testedHuntScore,"_EyeV",testedEyeV,".pdf",sep="")
      pdf(strPlotName,bg="white",
          compress=FALSE,onefile = FALSE, 
          title="Compare User vs Automated Hunt Event Detection  ") #col=(as.integer(filtereddatAllFrames$expID))
      
        mxAxis <- max(c(datCompEvents$AutomaticCount,datCompEvents$ManualCount))
        plot(datParamComp$ManualCount,datParamComp$AutomaticCount,xlim=c(0,mxAxis),ylim=c(0,mxAxis),asp=1,
             xlab="Manual Count",ylab="Automatic",
             main=paste(" HuntEvent Detection H:",testedHuntScore,"V:",testedEyeV))
        abline(lmmodel,col="red",lwd=3,lty=2)
        legend("bottomright",bg="white",legend=c(paste("LM c=",prettyNum(lmmodel$coefficients[1],digits=3),
                                            "b=",prettyNum(lmmodel$coefficients[2],digits=3),
                                            "r2=",prettyNum(rSq,digits=3) ) ) 
               ,lty=2,col="red",lwd=3)
        text(datParamComp$ManualCount,datParamComp$AutomaticCount+4,datParamComp$expID,cex=0.6)
        
      dev.off()
      
          
    }
  }
  ## Export Results from Parameter Scanning
  write.csv(datCompEvents,file=paste0(strDataExportDir,"/datHEventsDetectionAbility.csv") )
  
  ## Marginalize/Integrate Free param and Find mean sensitivity Specificity
  vmuSensitivity <- tapply(datCompEvents$Sensitivity,datCompEvents$ClassifierThres,mean)
  vmuSpecificity <- tapply(datCompEvents$Specificity,datCompEvents$ClassifierThres,mean)
  
  plot(1-vmuSpecificity,vmuSensitivity,type="l",
       main=paste("ROC HUNTSCORETHRES ",min(datCompEvents$ClassifierThres),"-",max(datCompEvents$ClassifierThres) ),
       xlim=c(0,1),ylim=c(0,1))
       
  ## Marginalize/Integrate Eye Threshold param and Find mean sensitivity Specificity
  vmuSensitivity <- tapply(datCompEvents$Sensitivity,datCompEvents$EyeVThres,mean)
  vmuSpecificity <- tapply(datCompEvents$Specificity,datCompEvents$EyeVThres,mean)
  
  plot(1-vmuSpecificity,vmuSensitivity,type="l",
       main=paste("ROC THRESHUNTVERGENCEANGLE ",min(datCompEvents$EyeVThres),"-",max(datCompEvents$EyeVThres) ),xlim=c(0,1),ylim=c(0,1))
  
    
  ## 

  plot(1-datCompEvents$Specificity,datCompEvents$Sensitivity,main="All ROC points")
  