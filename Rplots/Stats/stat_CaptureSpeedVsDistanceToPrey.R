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

strmodel_capspeedVsDistance <- "
var x_rand[2,2];

model {


##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model  as determined by mod flag
  c[i,1:2] ~ dmnorm(mu[mID[i]+1,],prec[mID[i]+1, , ]) ## data in column 1 and 2
  mID[i] ~ dbern(0.5) ##Se Gaussian class membership randomly
}


##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture (1 and 2)
for  (g in 1:2)
{
  prec[g,1:2,1:2] <- inverse(cov[g,,])
  
  cov[g,1,1] <- sigma[g,1]*sigma[g,1]
  cov[g,1,2] <- sigma[g,1]*sigma[g,2]*rho[g]
  cov[g,2,1] <- sigma[g,1]*sigma[g,2]*rho[g]
  cov[g,2,2] <- sigma[g,2]*sigma[g,2]
  
  ## Priors 
  sigma[g,1] ~ dunif(0,1) ##dist prey - Keep it broad within the expected limits 
  sigma[g,2] ~ dunif(0,20) ##the cap speed sigma 
  rho[g] ~ dunif(-1,1) ##The covar coefficient
}

  mu[1,1] ~ dnorm(0,0.01) ##Distance prey
  mu[1,2] ~ dnorm(10,0.01) ##cap speed

  mu[2,1] ~ dnorm(0.5,0.01) ##Distance prey
  mu[2,2] ~ dnorm(30,0.01) ##cap speed


## Synthesize data from the distribution
x_rand[1,] ~ dmnorm(mu[1,],prec[1,,])
x_rand[2,] ~ dmnorm(mu[2,],prec[2,,])

} "

strModelPDFFilename <- "/stat/UndershootAnalysis/stat_modelMixCaptureSpeedVsDistToPrey_Valid.pdf";
strDataPDFFileName <- "/stat/UndershootAnalysis/PreyDistanceCaptureSpeed_scatterValid.pdf"
strCaptSpeedDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelMixCaptureSpeed_Valid.pdf"

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

### Capture Speed vs Distance to prey ###
datDistanceVsStrikeSpeed_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],Validated= lFirstBoutPoints$NL[,"Validated"] ) )
datDistanceVsStrikeSpeed_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),Validated= lFirstBoutPoints$LL[,"Validated"] )
datDistanceVsStrikeSpeed_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),Validated= lFirstBoutPoints$DL[,"Validated"] )

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
steps <- 1500
str_vars <- c("mu","rho","sigma","x_rand","mID")
ldata_LF <- list(c=datDistanceVsStrikeSpeed_LL,N=NROW(datDistanceVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datDistanceVsStrikeSpeed_NL,N=NROW(datDistanceVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datDistanceVsStrikeSpeed_DL,N=NROW(datDistanceVsStrikeSpeed_DL)) ##Dry fed
ldata_ALL <- list(c=datDistanceVsStrikeSpeed_ALL,N=NROW(datDistanceVsStrikeSpeed_ALL)) ##Dry fed

jags_model_LF <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_LF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_LF, 500)
draw_LF=jags.samples(jags_model_LF,steps,thin=2,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_NF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_NF)
draw_NF=jags.samples(jags_model_NF,steps,thin=2,variable.names=str_vars)

##  DRY  Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_DF, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_DF, 500)
draw_DF=jags.samples(jags_model_DF,steps,thin=2,variable.names=str_vars)

## All groups combined data points
jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_ALL, 
                            n.adapt = 500, n.chains = 3, quiet = F)
update(jags_model_ALL, 500)
draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)


### Estimate  densities  ###
nContours <- 5
ntail <-2000
pBw   <- 0.02 


zLL <- kde2d(c(tail(draw_LF$mu[,1,,1],ntail)), c(tail(draw_LF$mu[,2,,1],ntail)),n=80)
zNL <- kde2d(c(tail(draw_NF$mu[,1,,1],ntail)), c(tail(draw_NF$mu[,2,,1],ntail)),n=80)
zDL <- kde2d(c(tail(draw_DF$mu[,1,,1],ntail)), c(tail(draw_DF$mu[,2,,1],ntail)),n=80)
zALL <- kde2d(c(tail(draw_ALL$mu[,1,,1],ntail)), c(tail(draw_ALL$mu[,2,,1],ntail)),n=80)


## Check out the covar coeffient , compare estimated densities

dLLb_rho <-density(tail(draw_LF$rho[,,1],ntail),kernel="gaussian",bw=0.05)
dNLb_rho <-density(tail(draw_NF$rho[,,1],ntail),kernel="gaussian",bw=0.05)
dDLb_rho <-density(tail(draw_DF$rho[,,1],ntail),kernel="gaussian",bw=0.05)
dALLb_rho <-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=0.05)
dALLb_rho[[2]] <-density(tail(draw_ALL$rho[2,,1],ntail),kernel="gaussian",bw=0.05)

#dLLb_rho[[2]]<-density(tail(draw_LF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
#dNLb_rho[[2]]<-density(tail(draw_NF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
#dDLb_rho[[2]]<-density(tail(draw_DF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
#dALLb_rho[[2]]<-density(tail(draw_ALL$rho[1,,1],ntail),kernel="gaussian",bw=0.1)


## Check out the dist to prey variance  , compare estimated densities
#dLLb_sigmaD <- list();dNLb_sigmaD<-list();dDLb_sigmaD<-list();dALLb_sigmaD<-list()
dLLb_sigmaD <-density(tail(draw_LF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
dNLb_sigmaD <-density(tail(draw_NF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
dDLb_sigmaD <-density(tail(draw_DF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
dALLb_sigmaD <-density(tail(draw_ALL$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)


## Check out the dist to prey variance  , compare estimated densities
#dLLb_sigmaD[[2]] <-density(tail(draw_LF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dNLb_sigmaD[[2]]<-density(tail(draw_NF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dDLb_sigmaD[[2]]<-density(tail(draw_DF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dALLb_sigmaD[[2]]<-density(tail(draw_ALL$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)


dLLb_sigmaC <- list();dNLb_sigmaC<- list();dDLb_sigmaC<-list();dALLb_sigmaC<-list()
dLLb_sigmaC<-density(tail(draw_LF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
dNLb_sigmaC<-density(tail(draw_NF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
dDLb_sigmaC<-density(tail(draw_DF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
dALLb_sigmaC<-density(tail(draw_ALL$sigma[,2,,1],ntail),kernel="gaussian",bw=1)



##Get the synthesized data:
plot(tail((draw_NF$x_rand[,1,,1]) , ntail),tail((draw_NF$x_rand[,2,,1]) , ntail),col=colourH[1])
points(tail((draw_LF$x_rand[,1,,1]) , ntail),tail((draw_LF$x_rand[,2,,1]) , ntail),col=colourH[2])
points(tail((draw_DF$x_rand[,1,,1]) , ntail),tail((draw_DF$x_rand[,2,,1]) , ntail),col=colourH[3])
points(tail((draw_ALL$x_rand[,2,,1]) , ntail),tail((draw_ALL$x_rand[,1,,1]) , ntail),col=colourH[4],pch=1,cex=1.6)

####################################
## PLot Model / Means and covariance ##
## Open Output PDF 

pdf(file= paste(strPlotExportPath,strModelPDFFilename,sep=""),width=14,height=7,title="A statistical model for Capture Strike speed / Undershoot Ratio")

outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 3.5
padj <- -23.0
las <- 1

layout(matrix(c(1,2,3,4),2,2, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

## Plot the mean of the 2D Models ##
ntail <- 1000
plot(tail(draw_NF$mu[1,1,,1],ntail),tail(draw_NF$mu[1,2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,0.6),ylim=c(10,60),ylab=NA,xlab=NA )
points(tail(draw_NF$mu[2,1,,1],ntail),tail(draw_NF$mu[2,2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,0.6),ylim=c(10,60),ylab=NA,xlab=NA )

points(tail(draw_LF$mu[1,1,,1],ntail),tail(draw_LF$mu[1,2,,1],ntail),col=colourH[2],pch=pchL[2])
points(tail(draw_LF$mu[2,1,,1],ntail),tail(draw_LF$mu[2,2,,1],ntail),col=colourH[2],pch=pchL[2])

points(tail(draw_DF$mu[1,1,,1],ntail),tail(draw_DF$mu[1,2,,1],ntail),col=colourH[3],pch=pchL[3])
points(tail(draw_DF$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[3],pch=pchL[3])

points(tail(draw_ALL$mu[,1,,1],ntail),tail(draw_ALL$mu[,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=0.8, line = 2.2, expression("Distance to Prey (mm) "~(delta) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec)  " ))

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  ),
                  bquote(All ~ '#' ~ .(ldata_ALL$N)  ) ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL, col=colourLegL)
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

## Plot the covariance ##
plot(dNLb_rho,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,3),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_rho,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rho,col=colourLegL[3],lwd=3,lty=3)
lines(dALLb_rho,col=colourLegL[4],lwd=3,lty=4)

legend("topright",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  ),
                  bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3),lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Cov. Capture speed to Prey Distance  ",rho) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

### aDDD DISTANCE TO PREY VARIANCE COMPARISON

plot(dNLb_sigmaD,col=colourLegL[1],xlim=c(0,0.5),lwd=3,lty=1,ylim=c(0,20),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_sigmaD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_sigmaD,col=colourLegL[3],lwd=3,lty=3)
lines(dALLb_sigmaD,col=colourLegL[4],lwd=3,lty=4)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Variance Prey Distance  ",delta) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )

### PloT CAPT SPEED VARIANCE 

plot(dNLb_sigmaC,col=colourLegL[1],xlim=c(0.0,30),lwd=3,lty=1,ylim=c(0,0.3),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_sigmaC,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_sigmaC,col=colourLegL[3],lwd=3,lty=3)
lines(dALLb_sigmaC,col=colourLegL[4],lwd=3,lty=4)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Variance Capture Speed  ") ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )


dev.off()

#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)



### PLOT EMPIRICAL 
####
########################################################
###        Distance Vs Capture speed               ###


densNL <-  kde2d(datDistanceVsStrikeSpeed_NL$DistanceToPrey, datDistanceVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datDistanceVsStrikeSpeed_LL$DistanceToPrey, datDistanceVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datDistanceVsStrikeSpeed_DL$DistanceToPrey, datDistanceVsStrikeSpeed_DL$CaptureSpeed,n=80)

covNL  <- cov(datDistanceVsStrikeSpeed_NL$DistanceToPrey,datDistanceVsStrikeSpeed_NL$CaptureSpeed)
covLL <- cov( datDistanceVsStrikeSpeed_LL$DistanceToPrey,datDistanceVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( datDistanceVsStrikeSpeed_DL$DistanceToPrey,datDistanceVsStrikeSpeed_DL$CaptureSpeed)



pdf(file= paste(strPlotExportPath,strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

plot(datDistanceVsStrikeSpeed_NL$DistanceToPrey, datDistanceVsStrikeSpeed_NL$CaptureSpeed,col=colourP[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA)
lFit <- lm(datDistanceVsStrikeSpeed_NL$CaptureSpeed ~ datDistanceVsStrikeSpeed_NL$DistanceToPrey)
abline(lFit,col=colourH[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[1],lty=2,lwd=3)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datDistanceVsStrikeSpeed_LL$DistanceToPrey, datDistanceVsStrikeSpeed_LL$CaptureSpeed,col=colourP[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA)
lFit <- lm(datDistanceVsStrikeSpeed_LL$CaptureSpeed ~ datDistanceVsStrikeSpeed_LL$DistanceToPrey)
abline(lFit,col=colourH[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[2],lty=2,lwd=3)
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


plot(datDistanceVsStrikeSpeed_DL$DistanceToPrey, datDistanceVsStrikeSpeed_DL$CaptureSpeed,col=colourP[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA)
lFit <- lm(datDistanceVsStrikeSpeed_DL$CaptureSpeed ~ datDistanceVsStrikeSpeed_DL$DistanceToPrey)
abline(lFit,col=colourH[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[3],lty=2,lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression("Distance To Prey "~(d) ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ) 


dev.off()



#### Capture Speed Only Model And Data ##

## pLOT THE Capture Speed

pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName,sep=""))

layout(matrix(c(1,2,3),3,1, byrow = FALSE))
xquant <- seq(0,60,1)
XLIM <- c(0,60)
pdistBW <- 2 ## mm/sec
strKern <- "gaussian"
ntail <- NROW(draw_NF$mu[1,2,,1])*0.10

plot(density(datDistanceVsStrikeSpeed_NL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM ,main="Capture Speed on capture strike")
for (i in 1:(ntail-1) )
{
  lines(xquant,dnorm(xquant,mean=tail(draw_NF$mu[1,2,ntail-i,1],1),sd=tail(draw_NF$sigma[1,2,ntail-i,1],1)),type='l',col=colourH[1] )
  lines(xquant,dnorm(xquant,mean=tail(draw_NF$mu[2,2,ntail-i,1],1),sd=tail(draw_NF$sigma[2,2,ntail-i,1],1)),type='l',col=colourH[1] )
}

lines(density(datTurnVsStrikeSpeed_NL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
legend("topright",title="NF",
       legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                 paste("model " ) ),
       col=c("black",colourH[3]),lwd=c(3,1) ) 


plot(density(datDistanceVsStrikeSpeed_LL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,main=NA)
for (i in 1:(ntail-1) )
{
  lines(xquant,dnorm(xquant,mean=tail(draw_LF$mu[1,2,ntail-i,1],1),sd=tail(draw_LF$sigma[1,2,ntail-i,1],1)),type='l',col=colourH[2] )
  lines(xquant,dnorm(xquant,mean=tail(draw_LF$mu[2,2,ntail-i,1],1),sd=tail(draw_LF$sigma[2,2,ntail-i,1],1)),type='l',col=colourH[2] )
}
lines(density(datTurnVsStrikeSpeed_LL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,main=NA)

legend("topright",title="LF",
       legend=c( paste("Data Density") , #(Bw:",prettyNum(digits=2, pdistBW ),")"
                 paste("model " ) ),
       col=c("black",colourH[2]),lwd=c(3,1) ) 


plot(density(datDistanceVsStrikeSpeed_DL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,main=NA)
for (i in 1:(ntail-1) )
{
  lines(xquant,dnorm(xquant,mean=tail(draw_DF$mu[1,2,ntail-i,1],1),sd=tail(draw_DF$sigma[1,2,ntail-i,1],1)),type='l',col=colourH[3] )
  lines(xquant,dnorm(xquant,mean=tail(draw_DF$mu[2,2,ntail-i,1],1),sd=tail(draw_DF$sigma[2,2,ntail-i,1],1)),type='l',col=colourH[3] )
}
lines(density(datTurnVsStrikeSpeed_DL$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,main=NA)

legend("topright",title="DF",
       legend=c( paste("Data  Density") , #(Bw:",prettyNum(digits=2, pdistBW ),")"
                 paste("model " ) ),
       col=c("black",colourH[3]),lwd=c(3,1) ) 

mtext(side = 1,cex=0.8, line = 2.2, expression("Distance To Prey (mm)" ))

dev.off()
embed_fonts(strDistDensityPDFFileName)
