##

### \notes main Tracking handles the initial import of tracked videos into data frames, and the detection of hunting episodes - 
## Hunting episdes are placed in as "setn15-HuntEvents-SB-Updated-Merged2 - and labelling updates this file.
##
### Hunt event Motion Analysis 
##  For close analysis of hunting the DataLabelling directory has a script to draw hunt events of as pecific assigned label or unlabelled ones
##  The output files of fish and prey (food)  are placed in HuntEventS_Retracked and from there they need to be moved to the subdirectory where LF, DF, NF 
##  retracked (Success / Fail separated) events are placed/ 
##  To import them run the HuntEpisodeAnalysis/runimportHuntEpisodeTrackFiles.r - This will merge food and larva motion records and produce
##  an imported hunt event register. From there analysis and a plot for each hunt episode is produced by running runHuntEpisodeAnalysis.r
#
## A script that Randomly and blindly allows to manually label hunt events exists in the DataLablling dir, main_LabellingBlind - this updates the 
## detected huntevent register ("setn15-HuntEvents-SB-Updated-Merged2)
#####################
#library(MASS)
library(RColorBrewer);

####################
##Convert RGB Values To HEX Colour Val.
col2hex <- function(rgbcol) {
  colPal <- vector()
  
  for (j in 1:NCOL(rgbcol))
  {
    
    colPal[j] <- rgb(rgbcol["red",j],rgbcol["green",j],rgbcol["blue",j],rgbcol["alpha",j],maxColorValue = 255)
  }
  return(colPal)
}

## TRACKER ROI ENDS Approximatelly 7.5 mm from dish egde
DIM_PXRADIUS       <- 790 #Is the DIAMETER Of the dish In the Video
DIM_PXDIAMETER       <- 790 #Is the DIAMETER Of the dish In the Video
DIM_MMPERPX        <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
DIM_DISHVOLUME     <- 1000*(0.035/2)^2*pi*0.01 ## in  Liters 
DIM_DISTTOMOUTH_PX <- 14 ## Estimated Distance from Centroid To Mouth based on head template size used in tracker
DIM_DISTTOMOUTH_MM <- DIM_DISTTOMOUTH_PX*DIM_MMPERPX ## Estimated Distance from CEntroid To Mouth based on head template size used in tracker
DIM_ROI_DIAMETER_MM  <- 515*DIM_MMPERPX
G_APPROXFPS              <- 60
G_THRESHUNTANGLE         <- 14 #Define Min Angle Both Eyes need to exceed for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 40 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_HUNTSCORETHRES         <- 0.65 ## Detection Thresh for DNN HUNT event detections (Based on fish image)
G_THRESHCLIPEYEDATA      <- 50 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- G_APPROXFPS/2
G_MINEPISODEDURATION     <- G_APPROXFPS/3
HUNTEVENT_MATCHING_OFFSET <- 2*G_APPROXFPS # Max frames to accept as mismatch when matching manual to auto detected huntevents during tracker validation - frames to Used in validateHuntEventsAgainstUserAnnotated
G_MIN_BOUTSPEED          <- 0.2 ##mm/frame - Need to be above to be considered A Motion Bout
G_THRES_CAPTURE_SPEED    <-  16 ###15##mm/sec ##Theshold on Body Speed above which a hunt event is marked to have a capture strike
G_THRES_MOTION_BOUT_SPEED <- 2.9 ##Got from Clustering #4 ##mm/sec
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments
G_MIN_TURNBOUT_ANGLE     <- 10 ##
G_THRES_TAILFQ_BOUT      <- 9.5 ##Hz
nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags
nEyeFilterWidth          <- nFrWidth

MIN_BOUT_DURATION        <- 10 ##Used in HuntEpisodeAnalysis_lib
MIN_BOUT_PAUSE           <- 25
G_MIN_BOUTSCORE          <- 1 ##Number of coincident event needed to detect bout in frame - (a combo of Tail Fq, Centroid Speed,Turn speed)
MIN_CAPTURE_EVENTS_PCA   <- 5 ## Make HuntBehaviour PCA Using Larvae that have at least 5 Hunt Events
MIN_BOUT_TURN_SPEED = 0.2
G_THRES_MOTION_BOUT_SPEED = 0.2
MIN_PATH_LENGTH     = 0.5 ##Also Set  in Previous Chunk
MIN_TURN_SIZE_THRES = 6 #Degrees - Reject As Non Turns If Overall Turn is less than 6 degrees
MIN_GAP_FRAMES_BETWEEN_TURNS = 50


## Plot Options ##
FONTSZ_AXISLAB <- 1.2
FONTSZ_AXIS    <- 1.2 ##Axis Tik Number Size

line      = 2.8 ## SubFig Label Params
lineAxis  = 2.7
lineTitle = 2.7
lineXAxis = 3.0
cex   = 1.4
adj   = 1.0
padj  <- -8.0
las   <- 1

rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
r <- c(rfc(7),"#FF0000");
pairedPalette <- col2rgb(brewer.pal(8,"Paired"),alpha = 1)
##For the 3 Groups 
###  NF, LF, DF , Black Colouring 
#pairedPalette["alpha",1:8] <- 210 ##Opacity
colourLegE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
colourLegL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
pairedPalette["alpha",1:8] <- 30 ##Opacity
colourHLine <- col2hex(pairedPalette[,c(6,4,2,5,3,1,7)]) ##Transparent For MCMC Samples (Live) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
pairedPalette["alpha",1:8] <-55 ##Opacity
colourHPoint <- col2hex(pairedPalette[,c(6,4,2,5,3,1,7)]) ##Transparent For MCMC Samples (Live) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)

colourDataScheme <- list()
colourDataScheme[["LF"]] <- list(Spont = colourLegE[2],Evoked = colourLegL[2])
colourDataScheme[["NF"]] <- list(Spont = colourLegE[1],Evoked = colourLegL[1])
colourDataScheme[["DF"]] <- list(Spont = colourLegE[3],Evoked = colourLegL[3])
colourDataScheme[["LL"]] <- colourLegL[2]
colourDataScheme[["LE"]] <- colourLegE[2]
colourDataScheme[["NL"]] <- colourLegL[1]
colourDataScheme[["NE"]] <- colourLegE[1]
colourDataScheme[["DL"]] <- colourLegL[3]
colourDataScheme[["DE"]] <- colourLegE[3]



pairedPalette["alpha",1:8] <- 35 ##Opacity
colourHE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
pairedPalette["alpha",1:8] <- 35 ##Opacity
colourHL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
colourH <- colourHL
colourP <- c(rgb(0.8,0.01,0.01,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.01,0.01,0.8,0.5),rgb(0.00,0.00,0.0,0.8)) ##points]
colourR <- c(rgb(0.9,0.01,0.01,0.6),rgb(0.01,0.7,0.01,0.6),rgb(0.01,0.01,0.9,0.6),rgb(0.1,0.1,0.1,0.6)) ##Region (Transparency)

pairedPalette["alpha",1:6] <- 200 ##Opacity
colourD <- col2hex(pairedPalette[,c(5,6,3,4,1,2)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
##<- c("#E60303AA","#03B303FF","#0303E6AA")
colourL        <-c("#03B303AF","#E60303AF","#0303E6AF")
colourClusters <- c("#00AFBB", "#E7B800", "#FC4E07")

pchL <- c(1,2,0,17,15,16,4) ## The style of bullet used for each group DL, LL, NL
pointTypeScheme <- list(DL=pchL[5], LL=pchL[4], NL=pchL[6],DE=pchL[3], LE=pchL[2], NE=pchL[1],LF=pchL[4],NF=pchL[6],DF=pchL[5])
lineTypeL <- c(2,1,3,4) ## The style of bullet used for each group DL, LL, NL
lineTypeScheme <- list(DL=lineTypeL[2], LL=lineTypeL[2], NL=lineTypeL[2], DE=lineTypeL[1], LE=lineTypeL[1], NE=lineTypeL[1], DF=lineTypeL[2], LF=lineTypeL[2], NF=lineTypeL[2])
lineTypeL.DF <- lineTypeL[1]
lineTypeL.LF <- lineTypeL[2]
lineTypeL.NF <- lineTypeL[3]

## Condition Labels
strDataLabels <- expression("NF-s","LF-s","DF-s","NF-e","LF-e","DF-e" )

## Uses Global assignment operateor <<- to set file locations depending on which system I am running the code on
setEnvFileLocations <- function(strSetName)
{
  
  if (strSetName == "HOME")
  {
    message("Set Global path parameters to ", strSetName, " environment")
    ## Required Variables - Locations -- Choose According To 
    # Home Desktop
    setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
    #strVideoFilePath  <<- "/media/kostasl/ARXEIO1TB/Behaviour/" 
    strVideoFilePath  <<- "/media/kostasl/zFish-Heta-T7/HungerExp"
    strTrackerPath    <<- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_15_0_GCC_64bit-Release" 
    strTrackeroutPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/Huntevents_retracked"#/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_Retracked/"
    #strTrackInputPath <<- "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/"
    strTrackInputPath <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/" 
    #strTrackInputPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/tracked_org/" 

    strDatDir         <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat" ##Where Are the Imported RData Stored
    #strDatDir         <<-  "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedOlivia/" ##Where Are the Imported RData Stored
    strDataExportDir  <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat"
    strDataStore      <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat" ##Where Large Data Is stored because Dropbox-Overflows
    strPlotExportPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat/plots" ##Where to source the Tracker csv files from 
  }
  
  
  if (strSetName == "OFFICE")
  {
      
    ## Office PC ##
    setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
    #strVideoFilePath  <<- "/media/kostasl/ARXEIO1TB/Behaviour/" 
    strVideoFilePath  <<- "/media/kostasl/zFish-Heta-T7/HungerExp"
    strTrackerPath    <<- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_15_1_GCC_64bit-Release" ## "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_15_0_GCC_64bit-Release" 
    strTrackeroutPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/Huntevents_retracked"#/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_Retracked/"
    #strTrackInputPath <<- "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/"
    #strTrackInputPath <<- "/media/kostasl/D445GB_ext4/expData/Olivia_assay/" 
    strTrackInputPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/" 
    
    strDatDir         <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat" ##Where Are the Imported RData Stored
    #strDatDir         <<-  "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedOlivia/" ##Where Are the Imported RData Stored
    strDataExportDir  <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat"
    strDataStore      <<-  "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat" ##Where Large Data Is stored because Dropbox-Overflows
    strPlotExportPath <<- "/media/kostasl/zFish-Heta-T7/HungerExp/tracked/Analysis/dat/plots" ##Where to source the Tracker csv files from #Where to source the Tracker csv files from 
  }
  
  if (strSetName == "LAPTOP")
  {
    ## Laptop ##
    setwd("~/workspace/zebrafishtrack/Rplots")
    strVideoFilePath  <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples"
    strTrackerPath    <<-  "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_11_2_GCC_64bit-Release"
    strTrackeroutPath <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/"
    strTrackInputPath <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/"##Where to source the Tracker csv files from 
    strDatDir         <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/" ##Where Are the Imported RData Stored
    strDataExportDir  <<- "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/Analysis/dat/"
    strDataStore      <<-  "/media/kostasl/zFish-Heta-T7/OliviaExp/Appetitesamples/tracked_org/Analysis/dat" ##Where Large Data Is stored because Dropbox-Overflows
    strPlotExportPath <<- "/mnt/data/Dropbox/Calculations/zebrafishtrackerData/plots"
    
  }  
  
  if (strSetName == "LAB")
  {
    ## Laptop ##
    setwd("~/workspace/zebrafishtrack/Rplots")
    strVideoFilePath  <<- "/mnt/Datastore/Olivia/Appetitesamples"
    strTrackerPath    <<-  "/home/meyerlab/workspace/build-zebraprey_track-Desktop_Qt_5_14_2_GCC_64bit-Release"
    strTrackeroutPath <<- "/mnt/Datastore/Olivia/Tracked"
    strTrackInputPath <<- "/mnt/Datastore/Olivia/Appetitesamples/tracked_org"##Where to source the Tracker csv files from 
    strDatDir         <<- "/mnt/Datastore/Olivia/Appetitesamples/tracked_org/Analysis/dat" ##Where Are the Imported RData Stored
    strDataExportDir  <<- "/mnt/Datastore/Olivia/Appetitesamples/tracked_org/Analysis/dat"
    strDataStore      <<-  "/mnt/Datastore/Olivia/Appetitesamples/tracked_org/Analysis/dat" ##Where Large Data Is stored because Dropbox-Overflows
    strPlotExportPath <<- "/mnt/Datastore/Olivia/Appetitesamples/tracked_org/Analysis/plots"
  }  
}
  
## GLOBAL VARS ###

## Util FUNCT ##
ProbValLessThan <- function(pPDF,X){ idxVal <- which(pPDF$x < X); return(sum(pPDF$y[idxVal][1:NROW(diff(pPDF$x[idxVal]))]*diff(pPDF$x[idxVal]) ));}
ProbValGreaterThan <- function(pPDF,X){ idxVal <- which(pPDF$x > X); return(sum(pPDF$y[idxVal][1:NROW(diff(pPDF$x[idxVal]))]*diff(pPDF$x[idxVal]) ));}
ProbValEqualZero <- function(pPDF){ idxVal <- which(pPDF$x >= -min(abs(pPDF$x)) & pPDF$x <= min(abs(pPDF$x))); return(sum(pPDF$y[idxVal][1:NROW(diff(pPDF$x[idxVal]))]*diff(pPDF$x[idxVal]) ));}


