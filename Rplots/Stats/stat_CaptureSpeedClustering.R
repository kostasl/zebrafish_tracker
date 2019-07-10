### Kostas Lagogiannis 2019-04-17 
## Discovered relationship between the last bout speed - ie Capture speed and the undershoot ratio -
## higher undershoot predicts higher capture speeds, while undershoot also seems to predict higher distance from prey 
##  suggesting that LF stays further away from prey, so it does stronger capture bouts and it is undershoot that allows it to do it.
## Their ability to judge distance is also revealed in the the eye vergence prior to capture, where there is a relationship between EyeV and distance to prey   is shown 
## Stohoi:
## S1 Establish whether undershoot covaries with capture speed
## S2 Compare cap.Speed vs UNdershoot models between groups - Do they also covary in all groups?
## S3 compare accuracy of capture speed vs distance to prey between groups (use covariance distributions)
## Aitiology :
## 
## A: Does undershoot explain capture speed and distance to prey accuracy?

### Stat Model on Capture speed vs undershoot
library(rjags)
library(runjags)

source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

strmodel_capspeedCluster <- "
model {

##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model  as determined by mod flag
  c[i] ~ dnorm(mu[mID[i]+1],1/(sigma[mID[i]+1])^2  ) ## data in column 1 and 2
  mID[i] ~ dbern(0.5) ##Se Gaussian class membership randomly
  
}

## Fit Bernouli distribution on Number of Hunt |Events that have a high-speed strike 
## Probability of Strike Swim 
pS  ~ dnorm(sum(mID)/N,1000)T(0,1)
mStrikeCount ~ dbin(pS,N )

##Prior
  ## Low Speed Captcha cluster
  mu[1] ~ dnorm(5,0.1)T(0,) ##cap speed
  sigma[1] ~ dunif(0,2) ##the low cap speed sigma 

  ## High speed Capture Cluster
  mu[2] ~ dnorm(35,0.1)T(mu[1],) ##cap speed
  sigma[2] ~ dunif(0,10) ##the high cap speed sigma 

} "


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotCaptureSpeedFit <- function(datSpeed,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,70,1)
  XLIM <- c(0,60)
  YLIM <- c(0,0.15)
  pdistBW <- 2 ## mm/sec
  strKern <- "gaussian"
  #ntail <- NROW(drawMCMC$mu[1,2,,nchain])*0.10
  ntail <- min(50,NROW(drawMCMC$mu[1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM,cex=cex,cex.axis=cex 
       ,main=NA,xlab = NA,ylab=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=2 )
  }
  
  dens<- density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern)
  lines(dens,col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",cex=cex,
         legend=c( paste0("",dens$n, "# Data density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourLegL[colourIdx],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
  mtext(side = 1,cex=cex, line = 3.2, expression("Capture speed (mm/sec) " ))
  mtext(side = 2,cex=cex, line = 2.5, expression("Density function " ))
  
}




strMainPDFFilename <- "/stat/UndershootAnalysis/fig4_stat_modelMixCaptureSpeedVsDistToPrey.pdf"; ## Used Fig 4
strModelVarPDFFilename <- "/stat/UndershootAnalysis/stat_modelMixCaptureSpeedVsDistToPrey_Variances.pdf";
strModelCoVarPDFFilename <- "/stat/UndershootAnalysis/fig4S1_stat_modelMixCaptureSpeedVsDistToPrey_COVariances.pdf";
strDataPDFFileName <- "/stat/UndershootAnalysis/fig4S2_PreyDistanceCaptureSpeed_scatterValid.pdf"
strClusterOccupancyPDFFileName <- "/stat/UndershootAnalysis/stat_modelCaptureStrike_ClusterOccupancy.pdf"

strCaptSpeedDensityPDFFileName <- "/stat/UndershootAnalysis/fig4_stat_modelMixCaptureSpeed_Valid.pdf" ## Used in Fig 4

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

### Capture Speed vs Distance to prey ###
datDistanceVsStrikeSpeed_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
datDistanceVsStrikeSpeed_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datDistanceVsStrikeSpeed_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )

###Subset Validated Only

###Validated Only
replace(datDistanceVsStrikeSpeed_NL$Validated, is.na(datDistanceVsStrikeSpeed_NL$Validated), 0)
replace(datDistanceVsStrikeSpeed_LL$Validated, is.na(datDistanceVsStrikeSpeed_LL$Validated), 0)
replace(datDistanceVsStrikeSpeed_DL$Validated, is.na(datDistanceVsStrikeSpeed_DL$Validated), 0) 

datDistanceVsStrikeSpeed_NL <- datDistanceVsStrikeSpeed_NL[datDistanceVsStrikeSpeed_NL$Validated == 1, ]
datDistanceVsStrikeSpeed_LL <- datDistanceVsStrikeSpeed_LL[datDistanceVsStrikeSpeed_LL$Validated == 1, ]
datDistanceVsStrikeSpeed_DL <- datDistanceVsStrikeSpeed_DL[datDistanceVsStrikeSpeed_DL$Validated == 1, ]

datDistanceVsStrikeSpeed_ALL <- rbind(datDistanceVsStrikeSpeed_NL,datDistanceVsStrikeSpeed_LL,datDistanceVsStrikeSpeed_DL)
##

##  Init  datastruct that we pass to model ##

##For Random allocation to model use: rbinom(n=10, size=1, prob=0.5)
steps <- 5500 #105500
str_vars <- c("mu","sigma","mID","mStrikeCount","pS","RegistarIdx")
ldata_LF <- list(c=datDistanceVsStrikeSpeed_LL$CaptureSpeed,N=NROW(datDistanceVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datDistanceVsStrikeSpeed_NL$CaptureSpeed,N=NROW(datDistanceVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datDistanceVsStrikeSpeed_DL$CaptureSpeed,N=NROW(datDistanceVsStrikeSpeed_DL)) ##Dry fed
ldata_ALL <- list(c=datDistanceVsStrikeSpeed_ALL$CaptureSpeed,N=NROW(datDistanceVsStrikeSpeed_ALL)) ##Dry fed

jags_model_LF <- jags.model(textConnection(strmodel_capspeedCluster), data = ldata_LF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_LF, 500)
draw_LF=jags.samples(jags_model_LF,steps,thin=2,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedCluster), data = ldata_NF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_NF)
draw_NF=jags.samples(jags_model_NF,steps,thin=2,variable.names=str_vars)

##  DRY  Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedCluster), data = ldata_DF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_DF, 500)
draw_DF=jags.samples(jags_model_DF,steps,thin=2,variable.names=str_vars)

## All groups combined data points
#jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_ALL, 
#                            n.adapt = 500, n.chains = 3, quiet = F)
#update(jags_model_ALL, 500)
#draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)

save(draw_LF,draw_NF,draw_DF,file =paste(strDataExportDir,"stat_CaptSpeedCluster_RJags.RData",sep=""))




### Load Pre Calc Results
load(file =paste(strDataExportDir,"stat_CaptSpeedCluster_RJags.RData",sep=""))
#### Main Figure 4 - Show Distance Vs Capture speed clusters for all groups - and Prob Of Capture Strike###


#######################################################
### PLOT  
####
########################################################
###        Distance Vs Capture speed               ###



#### FIG 4 / Capture Speed Only Model And Data ##
pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName ,sep=""))

par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),1,3, byrow = FALSE))
npchain<-3
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_NL,draw_NF,1,npchain)
#title(main="Model capture Speed")
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_LL,draw_LF,2,npchain)
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_DL,draw_DF,3,npchain)


dev.off()
#embed_fonts(strCaptSpeedDensityPDFFileName)









#### Plot Prey Location  ###########
## The Original list if the lFirstBout data from runHuntepisode analysis
source("plotTrackScatterAndDensities.r")
#plotCaptureBoutPreyPositions