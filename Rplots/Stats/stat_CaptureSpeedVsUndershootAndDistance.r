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

source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotCaptureSpeedFit <- function(datSpeed,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,70,1)
  XLIM <- c(0,60)
  YLIM <- c(0,0.15)
  pdistBW <- 2 ## mm/sec
  strKern <- "gaussian"
  #ntail <- NROW(drawMCMC$mu[1,2,,nchain])*0.10
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourR[4],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourH[colourIdx],lty=2 )
  }
  
  lines(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotUndeshootClusterFit <- function(datTurn,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,2,0.02)
  XLIM <- c(0,2)
  YLIM <- c(0,3)
  pdistBW <- 0.1 ## mm/sec
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,1,ntail-i,nchain],1)),type='l',col=colourR[4],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,1,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  
  lines(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


plotDistanceClustFit <- function(datDist,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(-0.1,0.8,0.05)
  pdistBW <- DIM_MMPERPX ## Manuall annotation  error is at least 1 px error , so smoothing with this bw is relevant
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  plot(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8),ylim=c(0,5) ,
       main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,3,ntail-i,nchain],1)),type='l',col=colourR[4],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,3,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  lines(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8) )
  legend("topright",title=NA,
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 

}


initfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(mID=c(rbinom(N,1,0.5)), ##Base Line Vergence Prior to HuntOn
#                                     sigma = matrix(c (  c(runif(1,min=0,max=0.1),runif(1,min=0,max=2)),
#s                                                         c(runif(1,min=0,max=0.1),runif(1,min=0,max=15))  ),nrow=2,byrow=T  ),
#                                     mu  = matrix(c (  c( rnorm(1,mean=1,sd=sqrt(1/10) ), rnorm(1,mean=8,sd=sqrt(1/2) ) ),
#                                                        c( rnorm(1,mean=1, sd=sqrt(1/10) ) , rnorm(1,mean=30, sd=sqrt(1/0.1) )    ) )
#                                                     ,nrow=2,byrow = T  ),
                                     ".RNG.name"="base::Super-Duper",
                                     ".RNG.seed"=round(runif(1,0,60000)) ),
                                     simplify=FALSE)
  return(initlist)
}


strmodel_capspeedVsUndershoot_Mixture <- "
var x_rand[2,3];

model {

##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model  as determined by mod flag
  c[i,1:3] ~ dmnorm(mu[mID[i]+1,],prec[mID[i]+1, , ]) ## data in column 1 and 2
  mID[i] ~ dbern(0.5) ##Se Gaussian class membership randomly
  
}
###14
## Fit Bernouli distribution on Number of Hunt |Events that have a high-speed strike 
## Probability of Strike Swim 
pS  ~ dnorm(sum(mID)/N,1000)T(0,1)
mStrikeCount ~ dbin(pS,N )

##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture (1 and 2)
for  (g in 1:2)
{
  prec[g,1:3,1:3] <- inverse(cov[g,1:3,1:3])
  
  cov[g,1,1] <- sigma[g,1]*sigma[g,1]
  cov[g,1,2] <- sigma[g,1]*sigma[g,2]*rho[g,1] ## Undershoot-Speed Covar
  cov[g,1,3] <- sigma[g,1]*sigma[g,3]*rho[g,3] ##Undeshoot-Dist Covar
  
  cov[g,2,1] <- sigma[g,1]*sigma[g,2]*rho[g,1] #UNdershoot-Speed
  cov[g,2,2] <- sigma[g,2]*sigma[g,2]
  cov[g,2,3] <- sigma[g,2]*sigma[g,3]*rho[g,2] #Speed-Dist Covar
  
  cov[g,3,1] <- sigma[g,1]*sigma[g,3]*rho[g,3] ##Undeshoot-Dist Covar
  cov[g,3,2] <- sigma[g,2]*sigma[g,3]*rho[g,2]
  cov[g,3,3] <- sigma[g,3]*sigma[g,3]

  #Sigmainv[g,1:3,1:3] ~ dwish(cov[g,,],3)
###37
  ##the sum of all the entries in a covariance matrix is the variance of the sum of the n random variables
  rho[g,1]  ~ dunif(-0.5,0.5) ##The Undershoot Speed covar coefficient
  rho[g,2] ~ dunif(-0.5,0.5) ##The Speed - Distance covar coefficient
  rho[g,3] ~ dunif(-0.5,0.5) ##The UNdershoot Distance covar coefficient

}
## Low Speed Captcha cluster

mu[1,1] ~ dnorm(1,0.00001)T(0.001,2) ##undershoot 
mu[1,2] ~ dnorm(5,0.1)T(0,) ## High cap speed
mu[1,3] ~ dnorm(0.5,0.01)T(0.0,) ##Distance prey

sigma[1,1] ~ dunif(0.0,0.20) ##Overshoot prey - Keep it broader within the expected limits
sigma[1,2] ~ dunif(0.00,4) ##the low cap speed sigma 
sigma[1,3] ~ dunif(0.0,1) ##dist prey - Keep it broad within the expected limits 
    

## High speed Capture Cluster
mu[2,1] ~ dnorm(1,0.00001)T(0.0,2) ##undershoot
mu[2,2] ~ dnorm(35,0.1)T(mu[1,2],) ##cap speed
mu[2,3] ~ dnorm(0.5,0.01)T(0,) ##Distance prey

sigma[2,1] ~ dunif(0.0,0.20) ##undershoot prey - Keep it narrow within the expected limits
sigma[2,2] ~ dunif(0.0,10) ##the cap speed sigma 
sigma[2,3] ~ dunif(0.0,1) ##dist prey - Keep it broad within the expected limits 

## Synthesize data from the distribution
x_rand[1,] ~ dmnorm(mu[1,],prec[1,,])
x_rand[2,] ~ dmnorm(mu[2,],prec[2,,])

} "


strModelPDFFileName <- "/stat/UndershootAnalysis/fig6-stat_modelCaptureSpeedVsUndershootAndDistance_Valid.pdf"
strDataPDFFileName <- "/stat/UndershootAnalysis/UndershootCaptureSpeedCV_scatter_Valid.pdf"
strCaptSpeedDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelCaptureSpeed_Valid.pdf"
strUndershootDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelUndershoot_Valid.pdf"
strDistanceDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelDistance_Valid.pdf"

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds","",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

datTurnVsStrikeSpeed_NL <- data.frame( cbind(Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],Validated= lFirstBoutPoints$NL[,"Validated"] )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],Validated= lFirstBoutPoints$DL[,"Validated"] )


###Validated Only
replace(datTurnVsStrikeSpeed_NL$Validated, is.na(datTurnVsStrikeSpeed_NL$Validated), 0)
replace(datTurnVsStrikeSpeed_LL$Validated, is.na(datTurnVsStrikeSpeed_LL$Validated), 0)
replace(datTurnVsStrikeSpeed_DL$Validated, is.na(datTurnVsStrikeSpeed_DL$Validated), 0) 

###Validated Only
datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$Validated == 1, ]
datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$Validated == 1, ]
datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$Validated == 1, ]

datTurnVsStrikeSpeed_ALL <- rbind(datTurnVsStrikeSpeed_NL,datTurnVsStrikeSpeed_LL,datTurnVsStrikeSpeed_DL)

##
##
steps <- 100000
nchains <- 5
nthin <- 2
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","rho","cov","sigma","x_rand","mID","mStrikeCount","pS") #Mixture Model
ldata_LF <- list(c=datTurnVsStrikeSpeed_LL,N=NROW(datTurnVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datTurnVsStrikeSpeed_NL,N=NROW(datTurnVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datTurnVsStrikeSpeed_DL,N=NROW(datTurnVsStrikeSpeed_DL)) ##Dry fed
ldata_ALL <- list(c=datTurnVsStrikeSpeed_ALL,N=NROW(datTurnVsStrikeSpeed_ALL)) ##Dry fed


jags_model_LF <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_LF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_LF$N))
update(jags_model_LF, 300)
draw_LF=jags.samples(jags_model_LF,steps,thin=nthin,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_NF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_NF$N)) 
update(jags_model_NF,300)
draw_NF=jags.samples(jags_model_NF,steps,thin=nthin,variable.names=str_vars)

## Dry  Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_DF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_DF$N))
update(jags_model_DF, 300)
draw_DF=jags.samples(jags_model_DF,steps,thin=nthin,variable.names=str_vars)

save.image(file = paste0(strDataExportDir,"stat_CaptSpeedVsUndershootAndDistance_RJags.RData"))
## ALL  groups
#jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_ALL, 
                            #n.adapt = 500, n.chains = 3, quiet = F)
#update(jags_model_ALL, 300)
#draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)

### Estimate  densities  ###

load(paste0(strDataExportDir,"stat_CaptSpeedVsUndershootAndDistance_RJags.RData"))

nContours <- 6
ntail <- 1200 #NROW(draw_NF$mu[1,1,,1])*0.20



zLL <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,2,,],ntail)),n=180)
zNL <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,2,,],ntail)),n=180)
zDL <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,2,,],ntail)),n=180)
#zALL <- kde2d(c(tail(draw_ALL$mu[,1,,1],ntail)), c(tail(draw_ALL$mu[,2,,1],ntail)),n=80)

zLLD <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,3,,],ntail)),n=180)
zNLD <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,3,,],ntail)),n=180)
zDLD <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,3,,],ntail)),n=180)


zLLS <- kde2d(c(tail(draw_LF$mu[,3,,],ntail)), c(tail(draw_LF$mu[,2,,],ntail)),n=180)
zNLS <- kde2d(c(tail(draw_NF$mu[,3,,],ntail)), c(tail(draw_NF$mu[,2,,],ntail)),n=180)
zDLS <- kde2d(c(tail(draw_DF$mu[,3,,],ntail)), c(tail(draw_DF$mu[,2,,],ntail)),n=180)

## Check out the covar coeffient , compare estimated densities
pBw   <- 0.1
## Strike Cluster Only (Fast speed) draw_NF$mu[,2,,]
## The Undershoot To Capt. Speed covar coefficient
dLLb_rhoUS<-density(tail(draw_LF$rho[2,1,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoUS<-density(tail(draw_NF$rho[2,1,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUS<-density(tail(draw_DF$rho[2,1,,],ntail),kernel="gaussian",bw=pBw)
#The Speed - Distance covar coefficient
dLLb_rhoSD<-density(tail(draw_LF$rho[2,2,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoSD<-density(tail(draw_NF$rho[2,2,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoSD<-density(tail(draw_DF$rho[2,2,,],ntail),kernel="gaussian",bw=pBw)
##The UNdershoot Distance covar 
dLLb_rhoUD<-density(tail(draw_LF$rho[2,3,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoUD<-density(tail(draw_NF$rho[2,3,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUD<-density(tail(draw_DF$rho[2,3,,],ntail),kernel="gaussian",bw=pBw)

#dALLb_rho<-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=pBw)

###Check COnv
draw <- draw_NF
plot(draw$mu[1,1,,1],type='l',ylim=c(0,2),col=rfc(nchains)[1] )
lines(draw$mu[1,1,,2],type='l',ylim=c(0,2),col=rfc(nchains)[2] )
lines(draw$mu[1,1,,3],type='l',ylim=c(0,2),col=rfc(nchains)[3] )
lines(draw$mu[1,1,,4],type='l',ylim=c(0,2),col=rfc(nchains)[4] )
lines(draw$mu[1,1,,5],type='l',ylim=c(0,2),col=rfc(nchains)[5] )
##Get the synthesized data:

#plot(tail(draw_NF$x_rand[1,,1],ntail ),tail(draw_NF$x_rand[2,,1],ntail ),col=colourH[1])
#points(tail(draw_LF$x_rand[1,,1],ntail ),tail(draw_LF$x_rand[2,,1],ntail ),col=colourH[2])
#points(tail(draw_DF$x_rand[1,,1],ntail ),tail(draw_DF$x_rand[2,,1],ntail ),col=colourH[3])

#################################### MAIN FIGURE #####
## PLot Model / Means and covariance ##
## Open Output PDF 
pdf(file= paste(strPlotExportPath,strModelPDFFileName,sep=""),width=14,height=7,title="A statistical model for Capture Strike speed / Undershoot Ratio")

### Show Speed Fit ###
outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 3.5
padj <- -10.0
las <- 1


layout(matrix(c(1,2,3,4),2,2, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

## Plot the mean of the 2D Models ##
##Collect Draws from all chains
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1], xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=0.8, line = 2.2, expression("Undershoot "~(gamma) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Capture Speed (mm/sec)  " ))
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex)

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)

legend("topright",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
                  ),
       pch=pchL, col=colourLegL)
###############

## Distance To Prey Vs Speed ##
plot(tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA )
points(tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=0.8, line = 2.2, expression("Distance to prey (mm)  " ))
mtext(side = 2,cex=0.8, line = 2.2, expression(" Capture Speed (mm/sec)" ))

contour(zDLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


#legend("topleft",
#       legend=c(  expression (),
#                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
#                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
#                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
#                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
#      ),
#       pch=pchL, col=colourLegL)
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex)


## Distance To Prey ##
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0.5,1.5),ylim=c(0,0.5),ylab=NA,xlab=NA )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=0.8, line = 2.2, expression("Undershoot "~(gamma) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Distance to prey (mm)  " ))

contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


#legend("topright",
#       legend=c(  expression (),
#                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
#                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
#                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
#                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
#       ),
#       pch=pchL, col=colourLegL)
mtext("C",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex)



### Show Cluster Membership
plot(density(tail(draw_NF$pS,ntail),pBw=0.05),col=colourLegL[1],xlim=c(0,1),ylim=c(0.4,10),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA)
lines(density(tail(draw_LF$pS,ntail) ),col=colourLegL[2],lwd=3,lty=2)
lines(density(tail(draw_DF$pS,ntail) ),col=colourLegL[3],lwd=3,lty=3)
#lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3)

mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Probability of high speed capture  ",(p["s"]) ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex)


dev.off()

##################################################
####################################################


##### Individual Rand Vars Fit ###
## Capture Speeds ##
pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
npchain<-3
plotCaptureSpeedFit(datTurnVsStrikeSpeed_NL,draw_NF,1,npchain)
title(main="Model capture Speed")
plotCaptureSpeedFit(datTurnVsStrikeSpeed_LL,draw_LF,2,npchain)
plotCaptureSpeedFit(datTurnVsStrikeSpeed_DL,draw_DF,3,npchain)

dev.off()

## Undershoot  ##
pdf(file= paste(strPlotExportPath,strUndershootDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotUndeshootClusterFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model undershoot on 1st turn to prey")
plotUndeshootClusterFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotUndeshootClusterFit(datTurnVsStrikeSpeed_DL,draw_DF,3)
dev.off()

## Distance ##
pdf(file= paste(strPlotExportPath,strDistanceDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotDistanceClustFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model distance from prey prior to capture")
plotDistanceClustFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotDistanceClustFit(datTurnVsStrikeSpeed_DL,draw_DF,3)

dev.off()
## plot 
##plot(xquant,dnorm(xquant,mean=tail(draw_NF$mu[2,2,,1],1),sd=tail(draw_NF$sigma[2,2,,1],1)),type='l',col=colourH[1],lty=1 )

### COVARIANCE PLOT ##
###Show covariance In the High Speed Capture Cluster ##
layout(matrix(c(1,2,3),1,3, byrow = FALSE))
## Plot the covariance Coefficients##
plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4),
     main="Undershoot-Speed", #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)
#lines(dALLb_rhoUS,col=colourLegL[4],lwd=3,lty=4)

plot(dNLb_rhoSD,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4),
     main="Speed-Distance", #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_rhoSD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoSD,col=colourLegL[3],lwd=3,lty=3)

##Speed TO Distance Covariance Coeff
plot(dNLb_rhoUD,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4),
     main="Undershoot-Distance", #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_rhoUD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUD,col=colourLegL[3],lwd=3,lty=3)



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

pdf(file= paste(strPlotExportPath,strDataPDFFileName,sep=""))
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

############# 3D

library(plot3D)
dev.off()
layout(matrix(c(1,2,3),1,3, byrow = FALSE))


pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/stat_UndershootSpeedDistance_3Dmodelplot_AB.pdf",sep=""))

scatter3D(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),grid=10,col = colourH[1],
          ticktype = "detailed",theta=0,phi=0,box=,type= c("shade", "wire", "dots"),zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),xlab="Undershoot", ylab="Distance",zlab="Speed")

scatter3D(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[2],
          zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)

scatter3D(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[3],
          zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)

dev.off()
###
#install.packages("plotly")

############# Plot Position Of Prey Prior Capture Bout 

pdf(file= paste(strPlotExportPath,"/PreyPositionPriorCapture_Validated.pdf",sep=""))

plotCaptureBoutPreyPositions()
dev.off()