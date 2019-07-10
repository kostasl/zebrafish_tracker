## Organize Manuscript Figures ### 
#### Kostas Lagogiannis 2019 
## \brief Make a scipt clarifying the script files used to produce each figure Used in the MS 



library(tools)
library(RColorBrewer);
library("MASS");
library(extrafont) ##For F
library(mvtnorm)

source("config_lib.R")
setEnvFileLocations("OFFICE") #OFFICE,#LAPTOP

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

### Capture Speed vs Distance to prey ###
datCapture_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
datCapture_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datCapture_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )


##Make A Contour From Params
xval=seq(1,10,0.1)
yval=seq(1,10,0.1)
grid=expand.grid(xval,yval)
m=matrix(dmvnorm(grid,mean=c(5,5),sigma=draw_NF$sigma[,,1,2]),length(xval),length(yval))
contour(m)

####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####

#### Plot Raw Capture Data Indicating Low/High Speed Clustering for each
### Load Pre Calc RJAgs Model Results
##   stat_CaptSpeedVsDistance_RJags.RData ##stat_CaptSpeedCluster_RJags.RData
load(file =paste(strDataExportDir,"stat_CaptSpeedVsDistance_RJags.RData",sep=""))


outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 3.2
lineXAxis = 3.0
cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1


### PLOT EMPIRICAL 
####
########################################################
###        Distance Vs Capture speed               ###
###
## Denote Fast/Slow CLuster Membership of Data Points - 
##Make List For Mean Number of Times Strike Was Classed as fast (score likelihood this is a fast one), and the RegIDx and Plot Point type,
minClusterLikelyhood <- 0.95 
steps <- NROW(draw_LF$mID[1,,1])
nsamples <- min(steps,1)

lClustScore_LF <- list(fastClustScore=apply(draw_LF$mID[,(steps-nsamples):nsamples,1],1,mean) ,RegistarIdx=datCapture_LL$RegistarIdx,pchL=rep_len(1,NROW(datCapture_LL)))
lClustScore_LF$pchL[lClustScore_LF$fastClustScore > minClusterLikelyhood] <- 16

lClustScore_NF <- list(fastClustScore=apply(draw_NF$mID[,(steps-nsamples):nsamples,1],1,mean) ,RegistarIdx=datCapture_NL$RegistarIdx,pchL=rep_len(1,NROW(datCapture_NL)))
lClustScore_NF$pchL[lClustScore_NF$fastClustScore > minClusterLikelyhood] <- 16

lClustScore_DF <- list(fastClustScore=apply(draw_DF$mID[,(steps-nsamples):nsamples,1],1,mean) ,RegistarIdx=datCapture_DL$RegistarIdx,pchL=rep_len(1,NROW(datCapture_DL)))
lClustScore_DF$pchL[lClustScore_DF$fastClustScore > minClusterLikelyhood] <- 16

##Make SPeed Density Of Each Cluster
dens_dist_NF_all <- density(datCapture_NL$DistanceToPrey)
dens_dist_NF_fast <- density(datCapture_NL$DistanceToPrey[lClustScore_NF$pchL == 16])
dens_dist_NF_slow <- density(datCapture_NL$DistanceToPrey[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_LF_all <- density(datCapture_LL$DistanceToPrey)
dens_dist_LF_fast <- density(datCapture_LL$DistanceToPrey[lClustScore_LF$pchL == 16])
dens_dist_LF_slow <- density(datCapture_LL$DistanceToPrey[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_DF_all <- density(datCapture_DL$DistanceToPrey)
dens_dist_DF_fast <- density(datCapture_DL$DistanceToPrey[lClustScore_DF$pchL == 16])
dens_dist_DF_slow <- density(datCapture_DL$DistanceToPrey[lClustScore_DF$pchL == 1])



plot(dens_dist_NF_all,xlim=c(0.0,0.5),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,5),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_dist_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_dist_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_dist_LF_all,xlim=c(0.0,0.5),col=colourLegL[2],lwd=4,lty=1)
lines(dens_dist_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_dist_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_dist_DF_all,xlim=c(0.0,0.5),col=colourLegL[3],lwd=4,lty=1)
lines(dens_dist_DF_fast,col=colourLegL[3],lwd=2,lty=2)
lines(dens_dist_DF_slow,col=colourLegE[3],lwd=2,lty=2)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(NROW(datCapture_NL))  ),
                  bquote(LF[""] ~ '#' ~ .(NROW(datCapture_LL))  ),
                  bquote(DF[""] ~ '#' ~ .(NROW(datCapture_DL))  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density ") )
mtext(side = 1,cex=cex, line = lineXAxis, expression(paste("Probability of high speed capture  ["~p["s"]~"]" ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)



###### Capture Speed  ###

lineAxis = 2.4
lineXAxis = 2.7
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,2,1))
##Make SPeed Density Of Each Cluster
dens_speed_NF_all <- density(datCapture_NL$CaptureSpeed)
dens_speed_NF_fast <- density(datCapture_NL$CaptureSpeed[lClustScore_NF$pchL == 16])
dens_speed_NF_slow <- density(datCapture_NL$CaptureSpeed[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_speed_LF_all <- density(datCapture_LL$CaptureSpeed)
dens_speed_LF_fast <- density(datCapture_LL$CaptureSpeed[lClustScore_LF$pchL == 16])
dens_speed_LF_slow <- density(datCapture_LL$CaptureSpeed[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_speed_DF_all <- density(datCapture_DL$CaptureSpeed)
dens_speed_DF_fast <- density(datCapture_DL$CaptureSpeed[lClustScore_DF$pchL == 16])
dens_speed_DF_slow <- density(datCapture_DL$CaptureSpeed[lClustScore_DF$pchL == 1])

## Plot Density Speed ##
plot(dens_speed_NF_all,xlim=c(0.0,60),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,0.1),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_speed_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_speed_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_speed_LF_all,xlim=c(0.0,60),col=colourLegL[2],lwd=4,lty=1,ylim=c(0,0.1))
lines(dens_speed_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_speed_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_speed_DF_all,xlim=c(0.0,60),col=colourLegL[3],lwd=4,lty=1,ylim=c(0,0.1))
lines(dens_speed_DF_fast,col=colourLegL[3],lwd=2,lty=2)
lines(dens_speed_DF_slow,col=colourLegE[3],lwd=2,lty=2)
####### END OF Speed ###

#'######### TURN RATIO ##########
##Make SPeed Density Of Each Cluster
dens_turn_NF_all <- density(datCapture_NL$Undershoot)
dens_turn_NF_fast <- density(datCapture_NL$Undershoot[lClustScore_NF$pchL == 16])
dens_turn_NF_slow <- density(datCapture_NL$Undershoot[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_turn_LF_all <- density(datCapture_LL$Undershoot)
dens_turn_LF_fast <- density(datCapture_LL$Undershoot[lClustScore_LF$pchL == 16])
dens_turn_LF_slow <- density(datCapture_LL$Undershoot[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
fracSlow_DF <- table(lClustScore_DF$pchL )[1]/NROW(lClustScore_DF$pchL)
dens_turn_DF_all <- density(datCapture_DL$Undershoot)
dens_turn_DF_fast <- density(datCapture_DL$Undershoot[lClustScore_DF$pchL == 16] )
dens_turn_DF_slow <- density(datCapture_DL$Undershoot[lClustScore_DF$pchL == 1])



## Plot TURN RATIO  ##
plot(dens_turn_NF_all,xlim=c(0.0,2),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,3),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_turn_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_turn_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_turn_LF_all,xlim=c(0.0,2),col=colourLegL[2],lwd=4,lty=1,ylim=c(0,3))
lines(dens_turn_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_turn_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_turn_DF_all,xlim=c(0.0,2),col=colourLegL[3],lwd=4,lty=1,ylim=c(0,3))
lines(dens_turn_DF_fast$x,dens_turn_DF_fast$y*(1-fracSlow_DF),col=colourLegL[3],lwd=2,lty=2)
lines(dens_turn_DF_slow$x,dens_turn_DF_slow$y*fracSlow_DF,col=colourLegE[3],lwd=2,lty=2)
####### END OF Speed ###



pdf(file= paste(strPlotExportPath,strDataPDFFileName,sep=""))

lineAxis = 2.4
lineXAxis = 2.7
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,2,1))
  
##For Colouring Data based on Fish/ExpID to examine Systematic differences
 colIdx_NL<- rainbow(max(colIdx_NL) )[as.numeric(as.factor(as.numeric(datTrackedEventsRegister[datCapture_NL$RegistarIdx,]$expID)))]
 colIdx_DL<- rainbow(max(colIdx_DL) )[as.numeric(as.factor(as.numeric(datTrackedEventsRegister[datCapture_DL$RegistarIdx,]$expID)))]
 colIdx_LL<- rainbow(max(colIdx_LL) )[as.numeric(as.factor(as.numeric(datTrackedEventsRegister[datCapture_LL$RegistarIdx,]$expID)))]
 
plot(datCapture_NL$DistanceToPrey, datCapture_NL$CaptureSpeed,col=colourP[1]  ,pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,1),main=NA,cex=cex)

lFit <- lm(datCapture_NL$CaptureSpeed[lClustScore_NF$pchL == 16] ~ datCapture_NL$DistanceToPrey[lClustScore_NF$pchL == 16])
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex  )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datCapture_LL$DistanceToPrey, datCapture_LL$CaptureSpeed,col=colourP[2],pch=lClustScore_LF$pchL,
     ylim=c(0,60),xlim=c(0,1),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datCapture_LL$CaptureSpeed[lClustScore_LF$pchL == 16] ~ datCapture_LL$DistanceToPrey[lClustScore_LF$pchL == 16])
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex  ) 


plot(datCapture_DL$DistanceToPrey, datCapture_DL$CaptureSpeed,col=colourP[3],pch=lClustScore_DF$pchL,
     ylim=c(0,60),xlim=c(0,1),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datCapture_DL$CaptureSpeed[lClustScore_DF$pchL == 16] ~ datCapture_DL$DistanceToPrey[lClustScore_DF$pchL == 16])
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Distance To Prey ["~d~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


dev.off()


############## Capture Speed Vs Turn Ratio #### 



pdf(file= paste(strPlotExportPath,"distal",strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.5,1,1))

plot(datCapture_NL$Undershoot, datCapture_NL$CaptureSpeed,col=colourLegL[1],pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datCapture_NL$CaptureSpeed ~ datCapture_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex)  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datCapture_LL$Undershoot, datCapture_LL$CaptureSpeed,col=colourLegL[2],pch=lClustScore_LF$pchL,
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datCapture_LL$CaptureSpeed ~ datCapture_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis-0.7, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datCapture_DL$Undershoot, datCapture_DL$CaptureSpeed,col=colourLegL[3],pch=lClustScore_DF$pchL,
     ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datCapture_DL$CaptureSpeed ~ datCapture_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


dev.off()



## EMPIRICAL - UNdeshoot vs Prey Distance 
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig7-UndershootDistanceCV_Distal_scatter.pdf",sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(4.5,4.3,0.5,1))

plot(datCapture_NL$Undershoot, datCapture_NL$DistanceToPrey,col=colourLegL[1],pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,1.0),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datCapture_NL$DistanceToPrey ~ datCapture_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datCapture_LL$Undershoot, datCapture_LL$DistanceToPrey,col=colourLegL[2],pch=lClustScore_LF$pchL,
     ylim=c(0,1),xlim=c(0,2.0),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datCapture_LL$DistanceToPrey ~ datCapture_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
mtext(side = 2,cex=cex, line = 2.2, expression("Distance to prey  (mm) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datCapture_DL$Undershoot, datCapture_DL$DistanceToPrey,col=colourLegL[3],pch=lClustScore_DF$pchL,
     ylim=c(0,1.0),xlim=c(0,2),   xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datCapture_DL$DistanceToPrey ~ datCapture_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex) 

dev.off()




########## END oF CaptureSpeed vs Distance  ###

##Fig 1  Epxperimental TimeLine manually designed  

### Fig 2A ####
## The kinematics was produced by selecting one of the figure produced
## from Hunt event analysis loop in : runHuntEpisodeAnalysis.r
## 

### Fig 2B ####
## 
#source("Stats/stat_HuntRateInPreyRange.R")
#source("Stats/stat_HuntDuration.R")


### Fig 3 ####
#source("Stats/stat_HuntEfficiency.r")

### Fig 4 ####
#source("DataLabelling/plotLabelledDataResults.R")
#source("Stats/stat_CaptureSpeedVsDistanceToPrey.R")

### Fig 5 ####
#source("Stats/stat_LinRegression_TurnVsBearing.R")

## Fig 6  clustering Speed, TurnRatio and Distance to Prey using 3D 2xGaussian mixture method####
#source("Stats/stat_ClusterCaptureSpeedVsUndershootAndDistance.r")

### Fig 7 Show Covariance using Gaussian 3D non clustering model aong wit ####
#source("Stats/stat_CaptureSpeedVsUndershootAndDistance.r")


