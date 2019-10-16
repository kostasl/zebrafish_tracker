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
library(MASS)


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


DIM_PXRADIUS <- 790 #Is the Radius Of the dish In the Video
DIM_MMPERPX <- 35/DIM_PXRADIUS ##35mm Opening of The viewport Assumed
DIM_DISTTOMOUTH_PX <- 14 ## Estimated Distance from CEntroid To Mouth based on head template size used in tracker
DIM_DISTTOMOUTH_MM <- DIM_DISTTOMOUTH_PX*DIM_MMPERPX ## Estimated Distance from CEntroid To Mouth based on head template size used in tracker
G_APPROXFPS              <- 410
G_THRESHUNTANGLE         <- 20 #Define Min Angle Both Eyes need for a hunting event to be assumed
G_THRESHUNTVERGENCEANGLE <- 45 ## When Eyes pointing Inwards Their Vergence (L-R)needs to exceed this value for Hunting To be considered
G_THRESHCLIPEYEDATA      <- 40 ##Limit To Which Eye Angle Data is filtered to lie within
G_MINGAPBETWEENEPISODES  <- 300
G_MINEPISODEDURATION     <- 100
G_MIN_BOUTSPEED          <- 0.2 ##mm/frame - Need to be above to be considered A Motion Bout
G_THRES_CAPTURE_SPEED    <-  16 ###15##mm/sec ##Theshold on Body Speed above which a hunt event is marked to have a capture strike
G_THRES_MOTION_BOUT_SPEED <- 2.9 ##Got from Clustering #4 ##mm/sec
PREY_COUNT_FRAMEWINDOW   <- 1600 ##Number oF Frames Over which to count Prey Stats at Beginning And End Of Experiments
G_MIN_TURNBOUT_ANGLE     <- 10 ##
G_THRES_TAILFQ_BOUT      <- 9.5 ##Hz
nFrWidth                 <- 20 ## Sliding Window Filter Width - Reduced From 50 to 20 to improve Meanf sliding window speed estimation lags
nEyeFilterWidth          <- nFrWidth*2
MIN_BOUT_DURATION        <- 10 ##Used in HuntEpisodeAnalysis_lib
MIN_BOUT_PAUSE           <- 25
G_MIN_BOUTSCORE          <- 1 ##Number of coincident event needed to detect bout in frame - (a combo of Tail Fq, Centroid Speed,Turn speed)

## Plot Options ##
FONTSZ_AXISLAB <- 1.2
FONTSZ_AXIS    <- 1.2 ##Axis Tik Number Size

line = 2.8 ## SubFig Label Params
lineAxis = 2.7
lineTitle = 2.7
lineXAxis = 3.0
cex = 1.4
adj  = 1.0
padj <- -8.0
las <- 1

rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
r <- c(rfc(7),"#FF0000");
pairedPalette <- col2rgb(brewer.pal(8,"Paired"),alpha = 1)
##For the 3 Groups 
###  NF, LF, DF , Black Colouring 
pairedPalette["alpha",1:8] <- 210 ##Opacity
colourLegE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
colourLegL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
pairedPalette["alpha",1:8] <-100 ##Opacity
colourHLine <- col2hex(pairedPalette[,c(6,4,2,5,3,1,7)]) ##Transparent For MCMC Samples (Live) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
pairedPalette["alpha",1:8] <-55 ##Opacity
colourHPoint <- col2hex(pairedPalette[,c(6,4,2,5,3,1,7)]) ##Transparent For MCMC Samples (Live) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)


pairedPalette["alpha",1:8] <- 5 ##Opacity
colourHE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
pairedPalette["alpha",1:8] <- 5 ##Opacity
colourHL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
colourH <- colourHL
colourP <- c(rgb(0.8,0.01,0.01,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.01,0.01,0.8,0.5),rgb(0.00,0.00,0.0,0.8)) ##points]
colourR <- c(rgb(0.9,0.01,0.01,0.6),rgb(0.01,0.7,0.01,0.6),rgb(0.01,0.01,0.9,0.6),rgb(0.1,0.1,0.1,0.6)) ##Region (Transparency)

pairedPalette["alpha",1:6] <- 200 ##Opacity
colourD <- col2hex(pairedPalette[,c(5,6,3,4,1,2)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
##<- c("#E60303AA","#03B303FF","#0303E6AA")
colourL <-c("#03B303AF","#E60303AF","#0303E6AF")

pchL <- c(1,2,0,17,15,16,4) ## The style of bullet used for each group DL, LL, NL
lineTypeL <- c(2,1,3,4) ## The style of bullet used for each group DL, LL, NL

## Condition Labels
strDataLabels <- expression("NF"["s"],"LF"["s"],"DF"["s"],"NF"["e"],"LF"["e"],"DF"["e"] )

## Uses Global assignment operateor <<- to set file locations depending on which system I am running the code on
setEnvFileLocations <- function(strSetName)
{
  if (strSetName == "HOME")
  {
    ## Required Variables - Locations -- Choose According To 
    # Home Desktop
    setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
    strVideoFilePath  <<- "//mnt/ARXEIO1TB/ExpData/zebrapreyCap/AnalysisSet/"
    strTrackerPath    <<- "/home/kostasl/workspace/build-zebraprey_track-Desktop-Release//" 
    strTrackeroutPath <<- "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_Retracked/"
    strTrackInputPath <<- "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/" ##Same As Working Dir
    strDatDir         <<-  "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
    strDataExportDir  <<-  "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
    strPlotExportPath <<- "/media/kostasl/D445GB_ext4/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots" ##Where to source the Tracker csv files from 
  }
  
  if (strSetName == "OFFICE")
  {
      
      ## Office PC ##
    setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
    strVideoFilePath  <<- "/media/LinuxDat/expDataKostas/AnalysisSetAlpha/" 
    strTrackerPath    <<- "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_11_1_GCC_64bit-Release/"
    strTrackeroutPath <<- "/media/LinuxDat/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_Retracked/"
    #strTrackInputPath <- "/mnt/570dce97-0c63-42db-8655-fbd28d22751d/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from
    strTrackInputPath <<- "/media/LinuxDat/TrackerOut/TrackASetRepeat/" ##Where to source the Tracker csv files from
    strDatDir         <<- "/media/LinuxDat/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
    strDataExportDir  <<- "/media/LinuxDat/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
    strPlotExportPath <<- "/media/LinuxDat/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots" ##Where to source the Tracker csv files from 
  }
  
  if (strSetName == "LAPTOP")
  {
    ## Laptop ##
    setwd("~/workspace/zebrafishtrack/Rplots")
    strVideoFilePath  <<- "/media/kostasl/zfishDataAlpha"
    strTrackerPath    <<-  "/home/kostasl/workspace/build-zebraprey_track-Desktop_Qt_5_11_2_GCC_64bit-Release"
    strTrackeroutPath <<- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/HuntEvents_Retracked/"
    strTrackInputPath <<- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData"##Where to source the Tracker csv files from 
    strDatDir         <<- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/dat/TrackedSessionA" ##Where Are the Imported RData Stored
    strDataExportDir  <<- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/out/"
    strPlotExportPath <<- "/home/kostasl/Dropbox/Calculations/zebrafishtrackerData/plots"
  }  
}
  
## GLOBAL VARS ###
