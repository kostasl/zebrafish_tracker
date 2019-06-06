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
G_THRES_MOTION_BOUT_SPEED <- 3.9 ##Got from Clustering #4 ##mm/sec
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

rfc <- colorRampPalette(rev(brewer.pal(8,'Spectral')));
r <- c(rfc(7),"#FF0000");
pairedPalette <- col2rgb(brewer.pal(8,"Paired"),alpha = 1)
##For the 3 Groups 
###  NF, LF, DF , Black Colouring 
pairedPalette["alpha",1:8] <- 210 ##Opacity
colourLegE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
colourLegL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
pairedPalette["alpha",1:8] <- 200 ##Opacity
colourHE <- col2hex(pairedPalette[,c(5,3,1,7)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
pairedPalette["alpha",1:8] <- 120 ##Opacity
colourHL <- col2hex(pairedPalette[,c(6,4,2,8)]) ##Transparent For MCMC Samples (Live)
colourH <- colourHL
colourP <- c(rgb(0.8,0.01,0.01,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.01,0.01,0.8,0.5),rgb(0.00,0.00,0.0,0.5)) ##points]
colourR <- c(rgb(0.9,0.01,0.01,0.6),rgb(0.01,0.7,0.01,0.6),rgb(0.01,0.01,0.9,0.6),rgb(0.1,0.1,0.1,0.6)) ##Region (Transparency)

pairedPalette["alpha",1:6] <- 200 ##Opacity
colourD <- col2hex(pairedPalette[,c(5,6,3,4,1,2)]) #c(rgb(0.95,0.01,0.01,0.1),rgb(0.01,0.7,0.01,0.1),rgb(0.01,0.01,0.9,0.1),rgb(0.00,0.00,0.0,1.0)) ####Transparent For MCMC Samples (Empty)
##<- c("#E60303AA","#03B303FF","#0303E6AA")
colourL <-c("#E60303AF","#03B303AF","#0303E6AF")

pchL <- c(1,2,0,5,16,17,4) ## The style of bullet used for each group DL, LL, NL
lineTypeL <- c(2,1,3,4) ## The style of bullet used for each group DL, LL, NL

## Condition Labels
strDataLabels <- expression("NF"["s"],"LF"["s"],"DF"["s"],"NF"["e"],"LF"["e"],"DF"["e"] )

## GLOBAL VARS ###
