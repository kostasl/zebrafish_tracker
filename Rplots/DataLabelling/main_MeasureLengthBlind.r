## Script Assisting to Measure Fish Length from each larval experiment - (in order to check correlation with success / Hunt Power)
# Kostasl 2 May 2020

library(tools)
library("MASS");


source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r")
source("Stats/stat_InformationTheoryAndCorrelations_bootstrap_lib.r")

#setEnvFileLocations("HOME") #HOME,OFFICE,#LAPTOP

resample <- function(x, ...) x[sample.int(length(x), ...)]

strProcDataFileName <- paste("setn15-HuntEvents-SB-Updated-Merged3") ##To Which To Save After Loading
message(paste(" Using Hunt Event List to Process... ",strProcDataFileName))

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
strGroupID <- c(levels(datTrackedEventsRegister$groupID),"DE","LE","NE")


datHuntEventAllGroupToLabel  <- getLabelledHuntEventsSet()
datFishSuccessRate <- getHuntSuccessPerFish(datHuntEventAllGroupToLabel)

## The scriptlet to run the labelling process on a set of expID is found in auxFunctions.r
datFlatPxLength <- readRDS(file= paste(strDataExportDir,"/FishLength_Updated3.rds",sep=""))
message(paste(" Loading Measured fish length in pixels data ... "))

 ##<- datHuntEvent
groupsList <- c("DL","NL","LL") ##unique(datHuntEventAllGroupToLabel$groupID)
str_FilterLabel <- vHuntEventLabels##c("UnLabelled","Success","Fail","Fail-No Strike")
#str_FilterLabel <- "NA"
##Select Randomly From THe Already Labelled Set ##
##Main Sample Loop
Keyc <- 'n'
while (Keyc != 'q')
{
  
  Keyc <- readline(prompt="### Press q to exit, 'n' for next, or type event number you wish to label  :")
  
  if (Keyc == 'q')
    break
  
  TargetLabel = which(vHuntEventLabels %in% str_FilterLabel)-1; ##Convert to Number Score
  gc <- resample(groupsList,1)
  idx <- NA
  TargetLabels <- vHuntEventLabels
  ##Select All non Measured ExpID - Within exp. cond tested with  Prey 
  vExpID_ToMeasure <- datFishSuccessRate[!(datFishSuccessRate$expID %in% datFlatPxLength[!is.na(datFlatPxLength$LengthPx),"expID"])
                                         & (datFishSuccessRate$groupID %in% groupsList) &
                                           datFishSuccessRate$HuntEvents > 0 ,"expID"]
  
  message("Larvae with observed hunt events")
  table(datFishSuccessRate[datFishSuccessRate$HuntEvents > 0,]$groupID)
  message(NROW(vExpID_ToMeasure), " Left to Measure")
  if ( NROW(vExpID_ToMeasure) == 0) 
  {
    message("Done Measuring all fish with hunt events in groups ",paste(groupsList,sep=", ") )
    break;
  }
    
  
  if (Keyc == 'n')
  {
    ##Choose From THe Set Of Videos Already Labelled From Another User (Kostasl) So as to Verify The Label # Sample Only From THose ExpID that have not been already verified
    #datHuntEventPool <- datHuntEventAllGroupToValidate[datHuntEventAllGroupToValidate$huntScore != "UnLabelled" & datHuntEventAllGroupToValidate$eventID != 0
    #                                           & (datHuntEventAllGroupToValidate$expID %in% datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$huntScore == TargetLabel,]$expID ),]
    
    datHuntEventPool <- datHuntEventAllGroupToLabel[datHuntEventAllGroupToLabel$eventID != 0 &
                                                      datHuntEventAllGroupToLabel$expID %in% vExpID_ToMeasure ,]
    datHuntEventPool <- datHuntEventPool[ datHuntEventPool$huntScore %in% TargetLabel ,] #& is.na(datHuntEventPool$markTracked)
    if (NROW(datHuntEventPool)  == 0)
    {
      message( paste("Finished with Hunt Events for labels ",TargetLabels[TargetLabel+1], ". Try Again") )
      groupsList <- groupsList[which(groupsList != gc)]
      next
    }
    ##Choose Random Exp/Event from Subset
    expID <- resample(datHuntEventPool$expID,1)
    datHuntEventPool <- datHuntEventPool[datHuntEventPool$expID == expID ,]
    eventID <- resample(datHuntEventPool$eventID,1)
    ###
    TargetLabels <- vHuntEventLabels[vHuntEventLabels %in% str_FilterLabel] ##Convert to Text Label Score to Use for Filtering OUt
  }
  ##Extract If Any Numbers In Input/ Then User Picked a specific Row
  if (!is.na(as.numeric(gsub("[^0-9]","",Keyc)) ) )
  {
    message(paste("Goto Event:",Keyc ) )
    idx <- as.character(Keyc) ##Note It acts as key only as string, numeric would just bring out the respective order idx record
    datHuntEventPool <- datHuntEventAllGroupToLabel[idx,]
    expID <- datHuntEventPool$expID
    eventID <- datHuntEventPool$eventID
    TargetLabels <- vHuntEventLabels
    
    if (is.na(datHuntEventAllGroupToLabel[idx,]$expID))
    {
      message("Event Not Found")
      next
    }
  }
  ##ExPORT 

  
  datHuntEventAllGroupToLabel <- labelHuntEvents(datHuntEventAllGroupToLabel,
                                                 strProcDataFileName,strVideoFilePath,
                                                 strTrackerPath,strTrackeroutPath,
                                                 TargetLabels,expID,eventID,idx,FALSE)
  
  ## Get User Measurement For Length
  KeyDat <- readline(prompt="### Type Fish Length In Pixels  (c exit):")
  
  if (KeyDat == 'c')
  {
    message("Stopping ", KeyDat)
    break;
  }
  
  groupID <- which(strGroupID == unique(datHuntEventPool[datHuntEventPool$expID == expID,]$groupID)  )
  datRow <- data.frame(KeyDat,groupID,expID)
  names(datRow) <- names(datFlatPxLength)
  datFlatPxLength <- rbind(datFlatPxLength,datRow)
  datFlatPxLength$LengthPx <- as.numeric(datFlatPxLength$LengthPx)
  saveRDS(datFlatPxLength,file= paste(strDataExportDir,"/FishLength_Updated3.rds",sep="") )
  message(paste(" Saved updated fish length to FishLength_Updated3.csv. "))

} ## Labelling LOOP

## Merge Length With Success
datFlatPxLength_filterNA <- datFlatPxLength[!is.na(datFlatPxLength$expID) & !is.na(datFlatPxLength$LengthPx), ]
datFlatPxMeanLength <- aggregate(datFlatPxLength_filterNA,
                                 by=list(m.expID=datFlatPxLength_filterNA$expID),"mean")
## Measured 
message("Larvae which we Measured their lengths ")
table(datFlatPxMeanLength$groupID)
message("Larvae with observed hunt events")
table(datFishSuccessRate[datFishSuccessRate$HuntEvents > 0,]$groupID)
## Now Merge Success With Lengths
datSuccessVsSize <- merge(datFlatPxMeanLength,datFishSuccessRate,by.x="m.expID",by.y="expID" )
datSuccessVsSize <- datSuccessVsSize[!is.nan(datSuccessVsSize$Efficiency),]
datSuccessVsSize.LF <- cbind(datSuccessVsSize[datSuccessVsSize$groupID.y == "LL",],Lengthmm=datSuccessVsSize[datSuccessVsSize$groupID.y == "LL",]$LengthPx*DIM_MMPERPX)
datSuccessVsSize.NF <- cbind(datSuccessVsSize[datSuccessVsSize$groupID.y == "NL",],Lengthmm=datSuccessVsSize[datSuccessVsSize$groupID.y == "NL",]$LengthPx*DIM_MMPERPX)
datSuccessVsSize.DF <- cbind(datSuccessVsSize[datSuccessVsSize$groupID.y == "DL",],Lengthmm=datSuccessVsSize[datSuccessVsSize$groupID.y == "DL",]$LengthPx*DIM_MMPERPX)
## Save To summary Stat Output - Used By generate figure 
saveRDS(datSuccessVsSize,file= paste(strDataExportDir,"/FishLengthVsHuntSuccess.rds",sep=""))

message("Correlation Of Efficiency To Size:")
message("(where hunt performance could affect past nutritional values )LF:",
        cor(datSuccessVsSize_LF$LengthPx*DIM_MMPERPX,datSuccessVsSize_LF$Efficiency,method="spearman"))
message("NF:",cor(datSuccessVsSize_NF$LengthPx*DIM_MMPERPX,datSuccessVsSize_NF$Efficiency,method="spearman"))

message("DF:",cor(datSuccessVsSize_DF$LengthPx*DIM_MMPERPX,datSuccessVsSize_DF$Efficiency,method="spearman"))

message("Efficiency Only Correlates with size in the LF group")
##EFFICIENCY 
plot(datSuccessVsSize.LF$LengthPx*DIM_MMPERPX,datSuccessVsSize.LF$Efficiency,col=colourLegL[2],pch=pchL[4],xlim=c(3.9,5),ylab=NA,xlab=NA)



## HUNT POWER LInear Fit 
pdf(paste0(strPlotExportPath,"/stat/fig3S1_stat_LarvalLengthsVsHPI.pdf"),width=7,height=7,title="Correlation Larval size to Hunt Success (HPI) ",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
  par(mar = c(3.9,4.7,1,1))
  lm.LF <- lm(data=datSuccessVsSize.LF,HuntPower~Lengthmm  )
  lm.NF <- lm(data=datSuccessVsSize.NF,HuntPower~Lengthmm  )
  lm.DF <- lm(data=datSuccessVsSize.DF,HuntPower~Lengthmm  )
  plot(datSuccessVsSize.LF$Lengthmm,datSuccessVsSize.LF$HuntPower,col=colourLegL[2],pch=pchL[4],xlim=c(3.9,5),ylab=NA,xlab=NA,ylim=c(0,5))
  abline(lm.LF,col=colourLegL[2],lwd=3)
  summary(lm.LF)
  
  points(datSuccessVsSize.DF$Lengthmm,datSuccessVsSize.DF$HuntPower,col=colourLegL[3],pch=pchL[5],xlim=c(3.9,5),ylab=NA,xlab=NA,ylim=c(0,5))
  abline(lm.DF,col=colourLegL[3],lwd=3)
  
  points(datSuccessVsSize.NF$Lengthmm,datSuccessVsSize.NF$HuntPower,col=colourLegL[1],pch=pchL[6],xlim=c(3.9,5),ylab=NA,xlab=NA)
  abline(lm.NF,col=colourLegL[1],lwd=3)
  mtext(side = 1,cex=cex,cex.main=cex, line = lineXAxis, expression(paste("Larval length (mm)  ") ))
  mtext(side = 2,cex=cex,cex.main=cex, line = lineAxis, expression("Hunt power index "))

dev.off()

hist(datSuccessVsSize_NF$LengthPx*DIM_MMPERPX*datSuccessVsSize_NF$Efficiency)
hist(datSuccessVsSize_LF$LengthPx*DIM_MMPERPX*datSuccessVsSize_LF$Efficiency)
hist(datSuccessVsSize_DF$LengthPx*DIM_MMPERPX*datSuccessVsSize_DF$Efficiency)

#plot(lm.LF)
summary(lm.NF)
summary(lm.DF)

cov(datSuccessVsSize_LF$LengthPx*DIM_MMPERPX,datSuccessVsSize_LF$Efficiency)
cov(datSuccessVsSize_DF$LengthPx*DIM_MMPERPX,datSuccessVsSize_DF$Efficiency)
cov(datSuccessVsSize_NF$LengthPx*DIM_MMPERPX,datSuccessVsSize_NF$Efficiency)

## Bootstrap correlation Analysis - Hunt Power Against Development/Nutrition Measured from Larval Std. Length
XRange <- c(3.9,5)
YRange <- c(0,5)
pBw <- 0.05
stat_SizeVsHuntPower_NF <- bootStrap_stat(datSuccessVsSize.NF$Lengthmm,datSuccessVsSize.NF$HuntPower,1000,XRange,YRange,"spearman")
stat_SizeVsHuntPower_LF <- bootStrap_stat(datSuccessVsSize.LF$Lengthmm,datSuccessVsSize.LF$HuntPower,1000,XRange,YRange,"spearman")
stat_SizeVsHuntPower_DF <- bootStrap_stat(datSuccessVsSize.DF$Lengthmm,datSuccessVsSize.DF$HuntPower,1000,XRange,YRange,"spearman")

pdf(paste0(strPlotExportPath,"/stat/fig3S1_stat_LarvalLengthsToHPI_Correlation.pdf"),width=7,height=7,title="Correlation Larval size to Hunt Success (HPI) ",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
  par(mar = c(3.9,4.7,1,1))

  plot(density(stat_SizeVsHuntPower_NF$corr,kernel="gaussian",bw=pBw),
       col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,7),main=NA, xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
  lines(density(stat_SizeVsHuntPower_LF$corr,kernel="gaussian",bw=pBw),col=colourLegL[2],lwd=3,lty=2)
  lines(density(stat_SizeVsHuntPower_DF$corr,kernel="gaussian",bw=pBw),col=colourLegL[3],lwd=3,lty=3)
  mtext(side = 1,cex=cex,cex.main=cex, line = lineXAxis, expression(paste("Correlation of hunt success (HPI) to larval length  ") ))
  mtext(side = 2,cex=cex,cex.main=cex, line = lineAxis, expression("Density function"))
  
  
  legend("topright",   legend=c( paste0("NF # ",  NROW(datSuccessVsSize_NF$expID) ),
                                  paste0("LF # " , NROW(datSuccessVsSize_LF$expID) ),
                                  paste0("DF # " , NROW(datSuccessVsSize_DF$expID) )
                                  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3),lwd=3,cex=cex)
  
dev.off() 
 
################## Labelling Record Manual Entry ###
## Extra Record ##
#3951 106.508
#3941 93.10
#3421 86.03
#3751 95.52
#3811 94.59
#3731 98.47
#3871 92.97
#284 83.8153
#3541  97.63
#4612 96.0417
#3931 91.54
#4211 95.27
#4001 91.41 91.44
#3413 91.23 91.70 94.37 91.09
#3961 89.2749 89.96 90.60
#3441 92.72 93.13
#4041 106.2 107.4
#3681 97.94 100.62 98.61 100.49
#4071 100.49 100.88
#4051 98.00 100.285 96.93 98
#3511  95.52 97.80 92.913 93.98 95.00
#3991 93.338 95.078 97 95.42
#3591 110.11 111.73
#3591 110.60
#4251 88.45 94.13 92.61 91.52
#4561 95.80 95.210 94.00
#252 88.60 91.06 89 92.1
#3921 95.33 91.08 93.34 95 
#3781 94.429 90.79 94.11 94
#301 100.28 97.61 101.17
#285 93.193 92.633 94.047
#2320 91.394 91.394 97.575 96.104 99.403 95.189 96.56
#240 96.834 96.208 96.840
#3721 97.406 96.881 94.260 97.83
#4551 99.29 96.93 96.76 
#4081 97.739 97.015 97.80
#306 87.664 92.097 92.309
#3971  88.600 92.347 88.45 92.11 89.80 5
#3801 103.368 100.84 103.47 101.6
#2951 91.547 93.557 93.74
#3471 100.623 101.356 101.02
#3981 91.923 94.255 96.426 94.403
#3691 93.8136 94.339 93.230
#282 92.200 92.913 89.5545
#3911 87.458 89.196 90.005
#277 95.524 98.08 93.9415
#4581 89.988 97.123 87.09 92.44
#4241 86.608 89.185 87.641
#4091 89.49 89.375 91.92
#254 97.046 104.48 101.789
#3551 90.4268 90.088 94.794 94.762 92.195 90.6863
#3662 92.444 91.78 92.89
#3581 93.6056 96.648 93.509 92.59
#3561 94.1329 92.13 99.609 92.574 98.4124 97.621
#3641 90.520 88.23 88.45 90.271
#319 99.2975 98.488 101.178
#4141 106.9 106.066 105.948 
#318 87.206 92.114 90.338 92.135 90.249
#3891 106.33 105.721 102.421 103.769
#287 93.4077 93.145 93.26
#4620 91.7061 92.6553 94.201 
#3600 98.234 100.12 98.183 100.24
#4031 99.1262 102.533  98.4073 100.424
#334 97.9439 99.0202 101.07
#4021 92.72 93.7763 92.1954
#4181 94.752 92.763 90.6863 93.086
#225 86.3713 90.4489 89.140 89.1403
#3821 93.230 95.210 95.0158
#228 103.238 107.019 108.167 107.448 106.17
#327 96.046 95.901 96.798
#326 97.6371 100 99.04
#3521 101.597 101.39 97.61 96.176
#3461 93.7443 94.7945 91.4822
#3831 99.247 99.282 99.085
#3451 94.8103 96.3328 94.201 95.462
#253 97.17 98.994 96.00
#224 93.96 91.26 92.005 93.300
#336 97.984 95.60 94.752 97.867 95.51
#288 100.76 96.896 97.989 101.51
#3771 93.94 94.868 93.605 
#3671 100.17 99.80 100.22 88.684 99.624
#3431 91.93 92.097 93.380
#229 89.82 90.75 88.865 90.70
#4361 89.738 91.92 91.04
#272 94.175 100.1 99.724 100
#4521 96.56 96.176 97.00 95.96
#3381 103.44 102.00 106.042 101.9 105.47
#3652 94.25 94.493 96.66 95.671
#298 94.641 93.477 92.195 95.414 
#3601 100 101.257 100.896
#226 92.800 93.407 90.91 90.443
#3701 91.181 89.69 88.61
#4121 96.260 96.046 93.107 93.637 91.54 94. 4
#2950 94.26 93.343 96.208 97.621 94.868
#263 102.12 102.318 104.0 102.20
#4621 85.0235 88.323 87.72  86.02
#3611 100.464 98.045 102.62 97.416 95.71 98.27
#309 92.59 90.05 89.89
#4391 92.54 94.89 94.37 93.477
#3741 102.92 103.23 102.30
#4191 98.249 100.24 99.36 99.92
#274 88.056 93.47 93.744 94.148
#308 97.08 97.693 99.045
#3851 99.02 100.00 103.07
#262  99.005 100.24 99.181 99.29
#3751 97.616 97.94 98.310 
#3761 95.273 94.66 96.40

#296 91.2853 89.185 92.135 90.426 90.603
#273 95.6033 92.195 96.150 91.934 94.339
#3571 95.25 95.189 96.332 
#3621 98.020 102.95 102.72 97.693
#239 91.934 92.95 94.021
#4161 102.46 97.86 98.84
#4451 104.3 103.46 104.86
#4151 91.678 88.19 89.988
#311 96.840 98.35 96.648
#3841 101.553 102.1 100
#283 92.417 93.150 91.96 91.21 90.603
#4461 93.648  92.347 92.655
#317 108.24 107.89 110.82
#236 104.48 104.019 102.201
#248 101.533 102.883 101.769
#--3501 92.6553
#4631 98.2497 93.637 98.0204 94.831
#3861 109.417 107.49 109.0 109.49
#238 90.9725 90.824 90.917 90.24
#231 90.354 91.706 92.957 91.394
#4011 99.93 100.125 99.403
#261 95.0789 94.429 94.847
#3881 88.0057 90.08 88.814 88.955
#297 92.179 91.416 95.692 94.260 93.861
#4061 89.88 92.784 92.309
#3501 93.962 95.21 93.493 90.80 92.617
#4341 100.18 97.514 98.615
#4321 99.724 97.862
#234 93.193 96.896 93.493 92.11
#337 95.341 97.308 97.949
#310 102.181 103.078 102.95
#4591 95.8801 97.128 98.858 97.128 96.0208
#4311 95.273  89.185 91.082 95.817 90.426 93.348 92.114
#260 102.489 102.24 104 104.4 104.73 103.96
#251 94.641 95.341 93.086 95.7706 
#4542 90.824 94.578 93.493 89.693 91.809
#4571 99.201 98.112 97.867 101.83
#4131 95.015 97.59 98.59 97.6217
#3531 101.045 97.416 98.081 98.00 97.185
#335  100.603 102 101.316 91.706 101.24
#3901 96.602 95.078 97.509 93.904 93.40 97.529
#4422 92.5743 94.593 94.25 93.680
#264 96.896 90.824 96.881 85.7963 95.50 94.752 95.566
#4231 101.11  102.303 100.846 99.724
#3481 97.862 95.859 96.602