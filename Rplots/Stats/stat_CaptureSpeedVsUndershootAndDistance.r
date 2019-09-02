### Kostas Lagogiannis 2019-06-24 
## 3D Gaussian Model for each group, to discover covariance structure in Undershoot to Distance/Speed
## I made this to complement the Clustering Method, so as to characterize the overall covariance structure

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
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourH[colourIdx],lty=2 )
  }
  ##Data
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
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,1,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,1,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
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
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,3,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,3,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
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

## Non Clustering Model, 3D Gaussian 
strmodel_capspeedVsUndershootAndDistance <- "
var x_rand[2,3];

model {

##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model  as determined by mod flag
  c[i,1:3] ~ dmnorm(mu[1,],prec[1, , ]) ## data in column 1 and 2
  mID[i] ~ dbern(0.5) ##Se Gaussian class membership randomly
  
}

##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture - Single Gaussian  Here -
for  (g in 1:1)
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


  ## Cluster's priors 
  mu[g,1] ~ dnorm(1, 0.00001)T(0.0,2) ##undershoot
  mu[g,2] ~ dnorm(15,0.00001)T(0,) ##cap speed
  mu[g,3] ~ dnorm(0.1,0.01)T(0,) ##Distance prey
  
  sigma[g,1] ~ dunif(0.0,0.30) ##undershoot prey - Keep it narrow within the expected limits
  sigma[g,2] ~ dunif(0.0,15) ## cap speed sigma 
  sigma[g,3] ~ dunif(0.0,0.3) ##dist prey - Keep it broad within the expected limits 

  ## Synthesize data from the distribution
  x_rand[g,] ~ dmnorm(mu[1,],prec[1,,])

}


} "


strModelPDFFileName <- "/stat/UndershootAnalysis/fig7S1-stat_modelCaptureSpeedVsUndershootAndDistance_Valid.pdf"
strDataPDFFileName <- "/stat/UndershootAnalysis/fig7-UndershootCaptureSpeedCV_scatter_Valid.pdf"
strCaptSpeedDensityPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelCaptureSpeed_Valid.pdf"
strUndershootDensityPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelUndershoot_Valid.pdf"
strDistanceDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelDistance_Valid.pdf"
strModelCovarPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelCaptureSpeedVsUndershootAndDistance_COVar.pdf"

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds","",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_Validated",".rds",sep="")) 

datTurnVsStrikeSpeed_NL <- data.frame( cbind(Undershoot=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],OnSetDistance=lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],Validated= lFirstBoutPoints$NL[,"Validated"] )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(Undershoot=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],OnSetDistance=lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(Undershoot=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],OnSetDistance=lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],Validated= lFirstBoutPoints$DL[,"Validated"] )

###Filter For Hunting Where Prey Is approached into strike distance, rather than Initial Prey Distance being within strike Distance
#datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$OnSetDistance > 0.6,]
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
steps <- 1000
nchains <- 5
nthin <- 2
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","rho","cov","sigma","x_rand","mID","mStrikeCount","pS") #Mixture Model
ldata_LF <- list(c=datTurnVsStrikeSpeed_LL,N=NROW(datTurnVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datTurnVsStrikeSpeed_NL,N=NROW(datTurnVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datTurnVsStrikeSpeed_DL,N=NROW(datTurnVsStrikeSpeed_DL)) ##Dry fed
ldata_ALL <- list(c=datTurnVsStrikeSpeed_ALL,N=NROW(datTurnVsStrikeSpeed_ALL)) ##Dry fed


jags_model_LF <- jags.model(textConnection(strmodel_capspeedVsUndershootAndDistance), data = ldata_LF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_LF$N))
update(jags_model_LF, 300)
draw_LF=jags.samples(jags_model_LF,steps,thin=nthin,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel_capspeedVsUndershootAndDistance), data = ldata_NF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_NF$N)) 
update(jags_model_NF,300)
draw_NF=jags.samples(jags_model_NF,steps,thin=nthin,variable.names=str_vars)

## Dry  Fed
jags_model_DF <- jags.model(textConnection(strmodel_capspeedVsUndershootAndDistance), data = ldata_DF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_DF$N))
update(jags_model_DF, 300)
draw_DF=jags.samples(jags_model_DF,steps,thin=nthin,variable.names=str_vars)

save(draw_NF,draw_LF,draw_DF,file = paste0(strDataExportDir,"stat_CaptSpeedVsUndershootAndDistance_RJags.RData"))
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
pBw   <- 0.02
## Strike Cluster Only (Fast speed) draw_NF$mu[,2,,]
## The Undershoot To Capt. Speed covar coefficient
dLLb_rhoUS<-density(tail(draw_LF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)  ## Undershoot-Speed Covar
dNLb_rhoUS<-density(tail(draw_NF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUS<-density(tail(draw_DF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)
#The Speed - Distance covar coefficient
dLLb_rhoSD<-density(tail(draw_LF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoSD<-density(tail(draw_NF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoSD<-density(tail(draw_DF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
##The UNdershoot Distance covar 
dLLb_rhoUD<-density(tail(draw_LF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoUD<-density(tail(draw_NF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUD<-density(tail(draw_DF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)

#dALLb_rho<-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=pBw)

save(dLLb_rhoSD,dNLb_rhoSD,dDLb_rhoSD,file = paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))
## ALL  


load(paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))
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




### MAIN COVARIANCE PLOT  (Fast Cluster)##
###Show covariance In the High Speed Capture Cluster ##

pdf(file= paste0(strPlotExportPath,strModelCovarPDFFileName),width=14,height=7,
    title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
nContours <- 5
### Show Speed Fit ###
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 2.7
lineXAxis = 3.0
lineTitle = 0.5

cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1

layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5),2,6, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,3.5,1))


  ## Plot the mean of the 2D Models ##
  ##Collect Draws from all chains
  plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1], xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
  points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
  points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
  #points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])
  
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Capture speed (mm/sec)  " ))
  mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
  
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
         pch=pchL, col=colourLegL,cex=cex)
  ###############
  
  
  ## Distance To Prey Vs Turn Ratio##
  plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0.5,1.5),ylim=c(0,0.5),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
  points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),col=colourHPoint[2],pch=pchL[2])
  points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),col=colourHPoint[3],pch=pchL[3])
  #points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])
  
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))
  
  contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
  contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
  contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
  
  mtext("B",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)



      ### COVARIANCES 
      ## Plot the covariance Coefficients##
      plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,4),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)
      #lines(dALLb_rhoUS,col=colourLegL[4],lwd=3,lty=4)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Capture speed and turn ratio "  ))
      mtext("C",at="topleft",outer=outer,side=2,col="black",font=2  ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
      
    
      plot(dNLb_rhoSD,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,4),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoSD,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoSD,col=colourLegL[3],lwd=3,lty=3)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Capture speed and distance"  ))
      
      mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
      
      ##Speed TO Distance Covariance Coeff
      plot(dNLb_rhoUD,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,4),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoUD,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoUD,col=colourLegL[3],lwd=3,lty=3)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Turn ratio and distance"  ))
      
      mtext("E",at="topleft",outer=F,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
      
dev.off()


### 3D density figure of Means ##
#library(MASS)
#library(plotly)
#den3d <- kde2d(x, y)
#persp(zLLD, box=FALSE)
## the new part:
#plot_ly(x=zLLD$x, y=zLLD$y, z=zLLD$z) %>% add_surface()

library( rgl )
library(plot3D)
# Static chart
ntail <- 150
datMu3D <-  data.frame( cbind.data.frame( TurnR=as.numeric(tail(draw_NF$mu[,1,,1],ntail)),CSpeed=tail(draw_NF$mu[,2,,1],ntail),Dist=tail(draw_NF$mu[,3,,1],ntail),col=colourHL[1])  )
datMu3D <- rbind(datMu3D,
                 data.frame( cbind.data.frame( TurnR=tail(draw_LF$mu[,1,,1],ntail),CSpeed=tail(draw_LF$mu[,2,,1],ntail),Dist=tail(draw_LF$mu[,3,,1],ntail),col=colourHL[2])  ))
datMu3D <- rbind(datMu3D,
                 data.frame( cbind.data.frame( TurnR=tail(draw_DF$mu[,1,,1],ntail),CSpeed=tail(draw_DF$mu[,2,,1],ntail),Dist=tail(draw_DF$mu[,3,,1],ntail),col=colourHL[3])  ))

rgl::plot3d( x=datMu3D$TurnR, y=datMu3D$CSpeed, z=datMu3D$Dist, col = datMu3D$col, type = "s", radius = 1.5,xlab="Turn Ratio",ylab="Capture Speed (mm/sec)",zlab="Distance to prey (mm)"
        )
rgl::rgl.postscript(paste0(strPlotExportPath,"fig7_Modelballs3D.pdf",fmt=pdf) )
decorate3d(xlim, ylim, zlim, 
           xlab = "x", ylab = "y", zlab = "z", 
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.03, 
           ...)



example(scatter3D)

####################################  Summary FIGURE Of Each Pair  #####
## PLot Model / Means and covariance ##
## Open Output PDF 
pdf(file= paste(strPlotExportPath,strModelPDFFileName,sep=""),width=14,height=7,
    title="A 3D statistical model for Capture Strike speed / Turn Ratio / Distance to Prey")

### Show Speed Fit ###
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 2.7
lineXAxis = 3.0
cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1


layout(matrix(c(1,2,3,4),2,2, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,2,1))

## Plot the mean of the 2D Models ##
##Collect Draws from all chains
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1], xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Capture Speed (mm/sec)  " ))
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)

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
       pch=pchL, col=colourLegL,cex=cex)
###############

## Distance To Prey Vs Speed ##
plot(tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA ,cex=cex,cex.axis=cex )
points(tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))
mtext(side = 2,cex=cex, line = lineAxis, expression(" Capture Speed (mm/sec)" ))

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
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)


## Distance To Prey Vs Turn Ratio##
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0.5,1.5),ylim=c(0,0.5),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))

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
mtext("C",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)



### Show Covar Of Undershoot to Distance  Membership
plot(dNLb_rhoUD,col=colourLegL[1],xlim=c(-1,1),ylim=c(0.4,10),lwd=4,lty=1,main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dLLb_rhoUD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUD,col=colourLegL[3],lwd=3,lty=3)
#lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density ") )
mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Turn ratio to distance covariance " ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)


dev.off()

##################################################
####################################################

### Show Covar Of Undershoot to Speed  Membership
plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-1,1),ylim=c(0.4,10),lwd=4,lty=1,main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)


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



## SLow Clust
clust <- 2
###Undershoot-Speed Covar
plot(density(draw_LF$sigma[clust,1,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,1,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,1,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,1,,1]),
     col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,1,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,1,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))

###Speed Distance
plot(density(draw_LF$sigma[clust,3,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,2,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,3,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,2,,1]),
      col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,3,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,2,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))


#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)

  ######################################################################
 ###         PLOT EMPIRICAL                              ##############
#####################################################################
###        UNdershoot Vs Capture speed               ###
densNL <-  kde2d(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,n=80)

covLL <- cov( 1/datTurnVsStrikeSpeed_LL$Undershoot,datTurnVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( 1/datTurnVsStrikeSpeed_DL$Undershoot,datTurnVsStrikeSpeed_DL$CaptureSpeed)
covNL  <- cov( 1/datTurnVsStrikeSpeed_NL$Undershoot,datTurnVsStrikeSpeed_NL$CaptureSpeed)


pdf(file= paste(strPlotExportPath,"distal",strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.5,1,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$CaptureSpeed ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex)  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,col=colourLegL[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$CaptureSpeed ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis-0.7, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,col=colourLegL[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$CaptureSpeed ~ datTurnVsStrikeSpeed_DL$Undershoot)
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

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$DistanceToPrey,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,1.0),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$DistanceToPrey ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$DistanceToPrey,col=colourLegL[2],
     ylim=c(0,1),xlim=c(0,2.0),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$DistanceToPrey ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
mtext(side = 2,cex=cex, line = 2.2, expression("Distance to prey  (mm) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$DistanceToPrey,col=colourLegL[3],
     ylim=c(0,1.0),xlim=c(0,2),   xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$DistanceToPrey ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex) 

dev.off()


### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig5S1-DetectionAngleDensity.pdf",sep=""))
  ### Show Speed Fit ###
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 2.7
  lineXAxis = 3.0
  cex = 1.4
  adj  = 3.5
  padj <- -8.0
  las <- 1
  
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.5,4.3,0.5,1))
  
  
  plot(density(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[1],xlim=c(-120.0,120),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA)
  lines(density(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[2],xlim=c(-120.0,120),lwd=3,lty=2)
  lines(density(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[3],xlim=c(-120.0,120),lwd=3,lty=3,main=NA)
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$NL))  ),
                    bquote(LF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$LL))  ),
                    bquote(DF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$DL))  )
         ), 
         col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)
  
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection  " ) )  )
  #mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=c
dev.off()


### Plot Initial Distance to Prey Vs Final Distance Scatter
plot(lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],lFirstBoutPoints$LL[,"DistanceToPrey"],col=colourLegL[2],pch=pchL[2])
points(lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],lFirstBoutPoints$NL[,"DistanceToPrey"],col=colourLegL[1],pch=pchL[1])
points(lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],lFirstBoutPoints$DL[,"DistanceToPrey"],col=colourLegL[3],pch=pchL[3])
############# 


### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig5S2-DetectionAngleVsDistance_scatter.pdf",sep=""))
  plot(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],col=colourLegL[2],pch=pchL[2],xlim=c(-120.0,120),lwd=2,lty=1,main=NA,xlab=NA,ylab=NA)
  points(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],col=colourLegL[1],pch=pchL[1],lwd=2)
  points(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],col=colourLegL[3],pch=pchL[3],lwd=2)
  mtext(side = 2,cex=cex, line = lineAxis, expression("Prey distance upon detection (mm)") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection (deg)  " ) )  )
dev.off()
##

# 
# library(plot3D)
# dev.off()
# layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# 
# 
# pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/stat_UndershootSpeedDistance_3Dmodelplot_AB.pdf",sep=""))
# 
# scatter3D(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),grid=10,col = colourH[1],
#           ticktype = "detailed",theta=0,phi=0,box=,type= c("shade", "wire", "dots"),zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),xlab="Undershoot", ylab="Distance",zlab="Speed")
# 
# scatter3D(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[2],
#           zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)
# 
# scatter3D(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[3],
#           zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)
# 
# dev.off()
###
#install.packages("plotly")

############# Plot Position Of Prey Prior Capture Bout 

pdf(file= paste(strPlotExportPath,"/PreyPositionPriorCapture_Validated.pdf",sep=""))

plotCaptureBoutPreyPositions()
dev.off()