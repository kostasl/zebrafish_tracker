### Kostas Lagogiannis 2019-04-15 
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


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

strmodel_capspeedVsUndershoot <- "
model {
  ##Draw capt speed from 2d gaussian
  for (i in 1:N)
  {
    c[i,1:2] ~ dmnorm(mu[],prec[ , ])
  }


  ##Covariance matrix and its inverse -> the precision matrix
  prec[1:2,1:2] <- inverse(cov[,])
  cov[1,1] <- sigma[1]*sigma[1]
  cov[1,2] <- sigma[1]*sigma[2]*rho
  cov[2,1] <- sigma[1]*sigma[2]*rho
  cov[2,2] <- sigma[2]*sigma[2]
  
  ## Priors 
  sigma[1] ~ dunif(0,1) ##the undershoot - Keep it broad within the expected limits 
  sigma[2] ~ dunif(0,100) ##the cap speed sigma 
  rho ~ dunif(-1,1) ##The covar coefficient
  mu[1] ~ dnorm(1,0.01) ##undershoot
  mu[2] ~ dnorm(0,0.01) ##cap speed
    
  ## Synthesize data from the distribution
  x_rand ~ dmnorm(mu[],prec[,])

} "

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_SetC",".rds",sep="")) 

datTurnVsStrikeSpeed_NL <- data.frame( cbind(Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]) )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]) )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]) )

##
steps <- 5000
str_vars <- c("mu","rho","sigma","x_rand")
ldata_LF <- list(c=datTurnVsStrikeSpeed_LL,N=NROW(datTurnVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datTurnVsStrikeSpeed_NL,N=NROW(datTurnVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datTurnVsStrikeSpeed_DL,N=NROW(datTurnVsStrikeSpeed_DL)) ##Dry fed

jags_model_LF <- jags.model(textConnection(strmodel_capspeedVsUndershoot), data = ldata_LF, 
                         n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_LF, 500)
draw_LF=jags.samples(jags_model_LF,steps,thin=2,variable.names=str_vars)

##Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedVsUndershoot), data = ldata_NF, 
                         n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_NF)
draw_NF=jags.samples(jags_model_NF,steps,thin=2,variable.names=str_vars)

##Not Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedVsUndershoot), data = ldata_DF, 
                         n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_DF, 500)
draw_DF=jags.samples(jags_model_DF,steps,thin=2,variable.names=str_vars)

### Estimate  densities  ###
nContours <- 5
zLL <- kde2d(c(tail(draw_LF$mu[1,,1],ntail)), c(tail(draw_LF$mu[2,,1],ntail)),n=80)
zNL <- kde2d(c(tail(draw_NF$mu[1,,1],ntail)), c(tail(draw_NF$mu[2,,1],ntail)),n=80)
zDL <- kde2d(c(tail(draw_DF$mu[1,,1],ntail)), c(tail(draw_DF$mu[2,,1],ntail)),n=80)

## Check out the covar coeffient , compare estimated densities
dLLb_rho<-density(tail(draw_LF$rho[,,1],ntail),kernel="gaussian",bw=pBw)
dNLb_rho<-density(tail(draw_NF$rho[,,1],ntail),kernel="gaussian",bw=pBw)
dDLb_rho<-density(tail(draw_DF$rho[,,1],ntail),kernel="gaussian",bw=pBw)


ntail <-2000
pBw   <- 0.1 
##Get the synthesized data:
plot(draw_NF$x_rand[1,(steps-ntail):steps,1],draw_NF$x_rand[2,(steps-ntail):steps,1],col=colourH[1])
points(draw_LF$x_rand[1,(steps-ntail):steps,1],draw_LF$x_rand[2,(steps-ntail):steps,1],col=colourH[2])
points(draw_DF$x_rand[1,(steps-ntail):steps,1],draw_DF$x_rand[2,(steps-ntail):steps,1],col=colourH[3])

####################################
## PLot Model / Means and covariance ##
## Open Output PDF 
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/stat_modelCaptureSpeedVsUndershoot_SetC2.pdf",sep=""),width=14,height=7,title="A statistical model for Capture Strike speed / Undershoot Ratio")

outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 3.5
padj <- -23.0
las <- 1

layout(matrix(c(1,2),1,2, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

## Plot the mean of the 2D Models ##
ntail <- 1000
plot(tail(draw_NF$mu[1,,1],ntail),tail(draw_NF$mu[2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,2),ylim=c(0,60),ylab=NA,xlab=NA )
points(tail(draw_LF$mu[1,,1],ntail),tail(draw_LF$mu[2,,1],ntail),col=colourH[2],pch=pchL[2])
points(tail(draw_DF$mu[1,,1],ntail),tail(draw_DF$mu[2,,1],ntail),col=colourH[3],pch=pchL[1])
mtext(side = 1,cex=0.8, line = 2.2, expression("Undershoot "~(gamma) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec)  " ))

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL, col=colourLegL)
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

## Plot the covariance ##
plot(dNLb_rho,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,3),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_rho,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rho,col=colourLegL[3],lwd=3,lty=3)
legend("topright",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3),lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Capture speed to Undershoot Covariance ",rho) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

dev.off()

#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)

### PLOT EMPIRICAL 
####
########################################################
###        UNdershoot Vs Capture speed               ###



densNL <-  kde2d(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,n=80)

covLL <- cov( 1/datTurnVsStrikeSpeed_LL$Undershoot,datTurnVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( 1/datTurnVsStrikeSpeed_DL$Undershoot,datTurnVsStrikeSpeed_DL$CaptureSpeed)
covNL  <- cov( 1/datTurnVsStrikeSpeed_NL$Undershoot,datTurnVsStrikeSpeed_NL$CaptureSpeed)


pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/UndershootCaptureSpeedCV_scatter.pdf",sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,col=colourP[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA)
lFit <- lm(datTurnVsStrikeSpeed_NL$CaptureSpeed ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourH[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[1],lty=2,lwd=3)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,col=colourP[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA)
lFit <- lm(datTurnVsStrikeSpeed_LL$CaptureSpeed ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourH[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[2],lty=2,lwd=3)
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,col=colourP[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA)
lFit <- lm(datTurnVsStrikeSpeed_DL$CaptureSpeed ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourH[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[3],lty=2,lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression("Undershoot "~(gamma) ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


dev.off()
