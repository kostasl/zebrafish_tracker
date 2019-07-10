## Organize Manuscript Figures ### 
#### Kostas Lagogiannis 2019 
## \brief Make a scipt clarifying the script files used to produce each figure Used in the MS 



library(tools)
library(RColorBrewer);
library("MASS");
library(extrafont) ##For F


source("config_lib.R")
setEnvFileLocations("LAPTOP") #OFFICE,#LAPTOP

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

### Capture Speed vs Distance to prey ###
datDistanceVsStrikeSpeed_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
datDistanceVsStrikeSpeed_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datDistanceVsStrikeSpeed_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"] )

###Subset Validated Only



####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####

#### Plot Raw Capture Data Indicating Low/High Speed Clustering for each
### Load Pre Calc RJAgs Model Results
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
lClustScore_LF <- list(fastClustScore=apply(draw_LF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_LL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_LL)))
lClustScore_LF$pchL[lClustScore_LF$fastClustScore > minClusterLikelyhood] <- 16

lClustScore_NF <- list(fastClustScore=apply(draw_NF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_NL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_NL)))
lClustScore_NF$pchL[lClustScore_NF$fastClustScore > minClusterLikelyhood] <- 16

lClustScore_DF <- list(fastClustScore=apply(draw_DF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_DL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_DL)))
lClustScore_DF$pchL[lClustScore_DF$fastClustScore > minClusterLikelyhood] <- 16

##Make SPeed Density Of Each Cluster
dens_dist_NF_fast <- density(datDistanceVsStrikeSpeed_NL$DistanceToPrey[lClustScore_NF$pchL == 16])
dens_dist_NF_slow <- density(datDistanceVsStrikeSpeed_NL$DistanceToPrey[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_LF_fast <- density(datDistanceVsStrikeSpeed_LL$DistanceToPrey[lClustScore_LF$pchL == 16])
dens_dist_LF_slow <- density(datDistanceVsStrikeSpeed_LL$DistanceToPrey[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_DF_fast <- density(datDistanceVsStrikeSpeed_DL$DistanceToPrey[lClustScore_DF$pchL == 16])
dens_dist_DF_slow <- density(datDistanceVsStrikeSpeed_DL$DistanceToPrey[lClustScore_DF$pchL == 1])



plot(dens_dist_NF_fast,xlim=c(-1.0,1),col=colourLegL[1],lwd=3,lty=1,ylim=c(0,5),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_dist_NF_slow,col=colourLegE[1],lwd=3,lty=1)

lines(dens_dist_LF_fast,col=colourLegL[2],lwd=3,lty=2)
lines(dens_dist_LF_slow,col=colourLegE[2],lwd=3,lty=2)

lines(dens_dist_DF_fast,col=colourLegL[3],lwd=3,lty=3)
lines(dens_dist_DF_slow,col=colourLegE[3],lwd=3,lty=3)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(NROW(datDistanceVsStrikeSpeed_NL))  ),
                  bquote(LF[""] ~ '#' ~ .(NROW(datDistanceVsStrikeSpeed_LL))  ),
                  bquote(DF[""] ~ '#' ~ .(NROW(datDistanceVsStrikeSpeed_DL))  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density ") )
mtext(side = 1,cex=cex, line = lineXAxis, expression(paste("Probability of high speed capture  ["~p["s"]~"]" ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)


pdf(file= paste(strPlotExportPath,strDataPDFFileName,sep=""))

lineAxis = 2.4
lineXAxis = 2.7
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,2,1))

plot(datDistanceVsStrikeSpeed_NL$DistanceToPrey, datDistanceVsStrikeSpeed_NL$CaptureSpeed,col=colourP[1],pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,1),main=NA,cex=cex)

lFit <- lm(datDistanceVsStrikeSpeed_NL$CaptureSpeed ~ datDistanceVsStrikeSpeed_NL$DistanceToPrey)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex  )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datDistanceVsStrikeSpeed_LL$DistanceToPrey, datDistanceVsStrikeSpeed_LL$CaptureSpeed,col=colourP[2],pch=lClustScore_LF$pchL,
     ylim=c(0,60),xlim=c(0,1),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datDistanceVsStrikeSpeed_LL$CaptureSpeed ~ datDistanceVsStrikeSpeed_LL$DistanceToPrey)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex  ) 


plot(datDistanceVsStrikeSpeed_DL$DistanceToPrey, datDistanceVsStrikeSpeed_DL$CaptureSpeed,col=colourP[3],pch=lClustScore_DF$pchL,
     ylim=c(0,60),xlim=c(0,1),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datDistanceVsStrikeSpeed_DL$CaptureSpeed ~ datDistanceVsStrikeSpeed_DL$DistanceToPrey)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Distance To Prey ["~d~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


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


