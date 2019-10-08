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

##### Notes on Covariance Prior , I could have used a wishard :
## What I do know is that if I assume Omega is my prior guess for covariance
## matrix Sigma, then I use the following in JAGS:
#  invSigma ~ dwish(J*Omega,J) #where J is the degrees of freedom, Omega Prior Guess for
#and I get the correct prior.

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

## Fit Bernouli distribution on Number of Hunt |Events that have a high-speed strike 
## Probability of Strike Swim 
pS  ~ dnorm(sum(mID)/N,1000)T(0,1)
mStrikeCount ~ dbin(pS,N )

##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture (1 and 2)
for  (g in 1:2)
{
  prec[g,1:2,1:2] <- inverse(cov[g,,])
  
  cov[g,1,1] <- sigma[g,1]*sigma[g,1]
  cov[g,1,2] <- sigma[g,1]*sigma[g,2]*rho[g]
  cov[g,2,1] <- cov[g,1,2] ##sigma[g,1]*sigma[g,2]*rho[g]
  cov[g,2,2] <- sigma[g,2]*sigma[g,2]
  
  ## Priors 
  sigma[g,1] ~ dunif(0,1) ##dist prey - Keep it broad within the expected limits 
  
  rho[g] ~ dunif(-1,1) ##The covar coefficient
}
  ## Low Speed Captcha cluster
  mu[1,1] ~ dnorm(0.5,0.01)T(0.0,) ##Distance prey
  mu[1,2] ~ dnorm(5,1)T(0,) ##cap speed
  sigma[1,2] ~ dunif(0,2) ##the low cap speed sigma 

  ## High speed Capture Cluster
  mu[2,1] ~ dnorm(0.5,0.01)T(0.0,) ##Distance prey ##precision=1/sigma^2
  mu[2,2] ~ dnorm(35,1)T(mu[1,2],) ##cap speed
  sigma[2,2] ~ dunif(0,10) ##the high cap speed sigma 

## Synthesize data from the distribution
x_rand[1,] ~ dmnorm(mu[1,],prec[1,,])
x_rand[2,] ~ dmnorm(mu[2,],prec[2,,])

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
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM,cex=cex,cex.axis=cex 
       ,main=NA,xlab = NA,ylab=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=2 )
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
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"//huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated",".rds",sep="")) 

### Capture Speed vs Distance to prey ###
datDistanceVsStrikeSpeed_NL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"],RegistarIdx=lFirstBoutPoints$NL[,"RegistarIdx"],Validated= lFirstBoutPoints$NL[,"Validated"],CaptureDuration= (lFirstBoutPoints$NL[,"CaptureStrikeEndFrame"]-lFirstBoutPoints$NL[,"CaptureStrikeFrame"])/G_APPROXFPS ) )
datDistanceVsStrikeSpeed_LL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],Validated= lFirstBoutPoints$LL[,"Validated"],CaptureDuration= (lFirstBoutPoints$LL[,"CaptureStrikeEndFrame"]-lFirstBoutPoints$LL[,"CaptureStrikeFrame"])/G_APPROXFPS  )
datDistanceVsStrikeSpeed_DL <- data.frame( cbind(DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),RegistarIdx=lFirstBoutPoints$DL[,"RegistarIdx"],Validated= lFirstBoutPoints$DL[,"Validated"],CaptureDuration= (lFirstBoutPoints$DL[,"CaptureStrikeEndFrame"]-lFirstBoutPoints$DL[,"CaptureStrikeFrame"])/G_APPROXFPS  )

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
str_vars <- c("mu","rho","sigma","cov","x_rand","mID","mStrikeCount","pS","RegistarIdx")
ldata_LF <- list(c=datDistanceVsStrikeSpeed_LL,N=NROW(datDistanceVsStrikeSpeed_LL)) ##Live fed
ldata_NF <- list(c=datDistanceVsStrikeSpeed_NL,N=NROW(datDistanceVsStrikeSpeed_NL)) ##Not fed
ldata_DF <- list(c=datDistanceVsStrikeSpeed_DL,N=NROW(datDistanceVsStrikeSpeed_DL)) ##Dry fed
ldata_ALL <- list(c=datDistanceVsStrikeSpeed_ALL,N=NROW(datDistanceVsStrikeSpeed_ALL)) ##Dry fed



### RUN MODEL ###
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
#jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsDistance), data = ldata_ALL, 
#                            n.adapt = 500, n.chains = 3, quiet = F)
#update(jags_model_ALL, 500)
#draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)

save(draw_LF,draw_NF,draw_DF,file =paste(strDataExportDir,"stat_CaptSpeedVsDistance_RJags.RData",sep=""))






### Load Pre Calc Results
load(file =paste(strDataExportDir,"stat_CaptSpeedVsDistance_RJags.RData",sep=""))
#### Main Figure 4 - Show Distance Vs Capture speed clusters for all groups - and Prob Of Capture Strike###

## Load COvariance (dLLb_rhoSD) - Calculated by 3D model in stat_CaptureSpeedVsUndershootAndDistance ##
load(file = paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))


#######################################################
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

#########
## Denote Fast/Slow CLuster Membership of Data Points - 
##Make List For Mean Number of Times Strike Was Classed as fast (score likelihood this is a fast one), and the RegIDx and Plot Point type,
lClustScore_LF <- list(fastClustScore=apply(draw_LF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_LL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_LL)))
lClustScore_LF$pchL[lClustScore_LF$fastClustScore > 0.7] <- 16

lClustScore_NF <- list(fastClustScore=apply(draw_NF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_NL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_NL)))
lClustScore_NF$pchL[lClustScore_NF$fastClustScore > 0.7] <- 16

lClustScore_DF <- list(fastClustScore=apply(draw_DF$mID[, (900):1000,1][,],1,mean) ,RegistarIdx=datDistanceVsStrikeSpeed_DL$RegistarIdx,pchL=rep_len(1,NROW(datDistanceVsStrikeSpeed_DL)))
lClustScore_DF$pchL[lClustScore_DF$fastClustScore > 0.7] <- 16

##Make Distance Density Of Each Cluster
dens_dist_NF_fast <- density(datDistanceVsStrikeSpeed_NL$DistanceToPrey[lClustScore_NF$pchL == 16])
dens_dist_NF_slow <- density(datDistanceVsStrikeSpeed_NL$DistanceToPrey[lClustScore_NF$pchL == 1])

##Make Distance Density Of Each Cluster
dens_dist_LF_fast <- density(datDistanceVsStrikeSpeed_LL$DistanceToPrey[lClustScore_LF$pchL == 16])
dens_dist_LF_slow <- density(datDistanceVsStrikeSpeed_LL$DistanceToPrey[lClustScore_LF$pchL == 1])


plot(dens_dist_NF_fast)
lines(dens_dist_NF_slow)

plot(dens_dist_LF_fast)
lines(dens_dist_LF_slow)



### Estimate  densities  ###
nContours <- 6
ntail <-2000
pBw   <- 0.02 

zLL <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,2,,],ntail)),n=180)
zNL <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,2,,],ntail)),n=180)
zDL <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,2,,],ntail)),n=180)
#zALL <- kde2d(c(tail(draw_ALL$mu[,1,,1],ntail)), c(tail(draw_ALL$mu[,2,,1],ntail)),n=80)


## Check out the covar coeffient , compare estimated densities

dLLb_rho_slow <-density(tail(draw_LF$rho[1,,1],ntail),kernel="gaussian",bw=0.05)
dNLb_rho_slow <-density(tail(draw_NF$rho[1,,1],ntail),kernel="gaussian",bw=0.05)
dDLb_rho_slow <-density(tail(draw_DF$rho[1,,1],ntail),kernel="gaussian",bw=0.05)

dLLb_rho_fast <-density(tail(draw_LF$rho[2,,1],ntail),kernel="gaussian",bw=0.05)
dNLb_rho_fast <-density(tail(draw_NF$rho[2,,1],ntail),kernel="gaussian",bw=0.05)
dDLb_rho_fast <-density(tail(draw_DF$rho[2,,1],ntail),kernel="gaussian",bw=0.05)

#dALLb_rho <-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=0.05)
##dALLb_rho[[2]] <-density(tail(draw_ALL$rho[2,,1],ntail),kernel="gaussian",bw=0.05)

dLLb_rho[[2]]<-density(tail(draw_LF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
dNLb_rho[[2]]<-density(tail(draw_NF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
dDLb_rho[[2]]<-density(tail(draw_DF$rho[2,,1],ntail),kernel="gaussian",bw=0.1)
#dALLb_rho[[2]]<-density(tail(draw_ALL$rho[1,,1],ntail),kernel="gaussian",bw=0.1)


## Check out the dist to prey variance  , compare estimated densities
#dLLb_sigmaD <- list();dNLb_sigmaD<-list();dDLb_sigmaD<-list();dALLb_sigmaD<-list()
dLLb_sigmaD <-density(tail(draw_LF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
dNLb_sigmaD <-density(tail(draw_NF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
dDLb_sigmaD <-density(tail(draw_DF$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)
#dALLb_sigmaD <-density(tail(draw_ALL$sigma[,1,,1],ntail),kernel="gaussian",bw=pBw)


## Check out the dist to prey variance  , compare estimated densities
#dLLb_sigmaD[[2]] <-density(tail(draw_LF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dNLb_sigmaD[[2]]<-density(tail(draw_NF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dDLb_sigmaD[[2]]<-density(tail(draw_DF$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)
#dALLb_sigmaD[[2]]<-density(tail(draw_ALL$sigma[2,2,,1],ntail),kernel="gaussian",bw=pBw)


dLLb_sigmaC <- list();dNLb_sigmaC<- list();dDLb_sigmaC<-list();dALLb_sigmaC<-list()
dLLb_sigmaC<-density(tail(draw_LF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
dNLb_sigmaC<-density(tail(draw_NF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
dDLb_sigmaC<-density(tail(draw_DF$sigma[,2,,1],ntail),kernel="gaussian",bw=1)
#dALLb_sigmaC<-density(tail(draw_ALL$sigma[,2,,1],ntail),kernel="gaussian",bw=1)

pdf(file= paste(strPlotExportPath,strMainPDFFilename,sep=""),width=14,height=7,
    title="A Gaussian clustering statistical model for capture strike speed and distance to prey")

outer = FALSE
line = 2.8 ## SubFig Label Params
lineAxis = 2.7
lineTitle = 2.7
lineXAxis = 3.0
cex = 1.4
adj  = 1.0
padj <- -8.0
las <- 1
nContours <- 5
npchain<-3

layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6),2,6, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
#par(mar = c(5,4.5,3,1))
par(mar = c(4.5,4.7,2,1))

plotCaptureSpeedFit(datDistanceVsStrikeSpeed_NL,draw_NF,1,npchain)
mtext("B",at="topleft",outer=F,side=2,col="black",font=2,  las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
#title(main="Model capture Speed")
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_LL,draw_LF,2,npchain)
mtext("C",at="topleft",outer=F,side=2,col="black",font=2,  las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_DL,draw_DF,3,npchain)
mtext("D",at="topleft",outer=F,side=2,col="black",font=2,  las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)


#### ## Probability Density of Strike capture ####
plot(density(draw_NF$pS,pBw=0.05),col=colourLegL[1],xlim=c(0,1),ylim=c(0.4,10),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA,
     cex=cex,cex.axis=cex )
lines(density(draw_LF$pS),col=colourLegL[2],lwd=3,lty=2)
lines(density(draw_DF$pS),col=colourLegL[3],lwd=3,lty=3)
#lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)
mtext(side = 1,cex=cex, line = lineXAxis, expression(paste("Probability of high speed capture  ["~p["s"]~"]" ) ) ,cex.main=cex )
mtext(side = 2,cex=cex, line = lineAxis, expression("Density  function" ))

mtext("E",at="topleft",outer=F,side=2,col="black",font=2,  las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
#### ## Probability Density of Strike capture ####
legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
                  #, bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)



## Plot the mean of the 2D Models Cluster ##
ntail <- 2000
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1],
     xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex )
#points(tail(draw_NF$mu[2,1,,1],ntail),tail(draw_NF$mu[2,2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
#points(tail(draw_LF$mu[2,1,,1],ntail),tail(draw_LF$mu[2,2,,1],ntail),col=colourH[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_DF$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[3],pch=pchL[3])

#points(tail(draw_ALL$mu[,1,,1],ntail),tail(draw_ALL$mu[,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineXAxis, expression("Distance to prey  ["~delta~"] (mm)" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Capture speed (mm/sec)  " ))
mtext("F",at="topleft",outer=F,side=2,col="black",font=2,     las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[3],lty=2)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
                  #bquote(All ~ '#' ~ .(ldata_ALL$N)  )
                  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL, col=colourLegL,cex=cex)



## Plot COvariance - Calculated by 3D model in stat_CaptureSpeedVsUndershootAndDistance
plot(dNLb_rhoSD,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
lines(dLLb_rhoSD_fast,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoSD,col=colourLegL[3],lwd=3,lty=3)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
mtext(side = 3,cex=cex, line = lineTitle-3, expression("Capture distance and speed "  ))
mtext("G",at="topleft",outer=F,side=2,col="black",font=2,     las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)

dev.off()







#### FIG 4 / Capture Speed Only Model And Data ##
pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName ,sep=""),width=14)

par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),1,3, byrow = FALSE))
npchain<-3
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_NL,draw_NF,1,npchain)
#title(main="Model capture Speed")
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_LL,draw_LF,2,npchain)
plotCaptureSpeedFit(datDistanceVsStrikeSpeed_DL,draw_DF,3,npchain)


dev.off()
#embed_fonts(strCaptSpeedDensityPDFFileName)
















##Get the synthesized data:
#plot(tail((draw_NF$x_rand[,1,,1]) , ntail),tail((draw_NF$x_rand[,2,,1]) , ntail),col=colourH[1])
#points(tail((draw_LF$x_rand[,1,,1]) , ntail),tail((draw_LF$x_rand[,2,,1]) , ntail),col=colourH[2])
#points(tail((draw_DF$x_rand[,1,,1]) , ntail),tail((draw_DF$x_rand[,2,,1]) , ntail),col=colourH[3])
#points(tail((draw_ALL$x_rand[,2,,1]) , ntail),tail((draw_ALL$x_rand[,1,,1]) , ntail),col=colourH[4],pch=1,cex=1.6)

####################################
## PLot Model / Means and covariance ##
## Open Output PDF 
message(strModelVarPDFFilename)
pdf(file= paste(strPlotExportPath,strModelVarPDFFilename,sep=""),width=14,height=7,
    title="A statistical model for Capture Strike speed And Distance to prey")

outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 3.5
padj <- -28.0
las <- 1

layout(matrix(c(1,2,3,4),2,2, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))

## Plot the mean of the 2D Models ##
ntail <- 600
plot(tail(draw_NF$mu[1,1,,1],ntail),tail(draw_NF$mu[1,2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,0.6),ylim=c(10,60),ylab=NA,xlab=NA )
points(tail(draw_NF$mu[2,1,,1],ntail),tail(draw_NF$mu[2,2,,1],ntail),col=colourH[1],pch=pchL[1], xlim=c(0,0.6),ylim=c(10,60),ylab=NA,xlab=NA )

points(tail(draw_LF$mu[1,1,,1],ntail),tail(draw_LF$mu[1,2,,1],ntail),col=colourH[2],pch=pchL[2])
points(tail(draw_LF$mu[2,1,,1],ntail),tail(draw_LF$mu[2,2,,1],ntail),col=colourH[2],pch=pchL[2])

points(tail(draw_DF$mu[1,1,,1],ntail),tail(draw_DF$mu[1,2,,1],ntail),col=colourH[3],pch=pchL[3])
points(tail(draw_DF$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[3],pch=pchL[3])

#points(tail(draw_ALL$mu[,1,,1],ntail),tail(draw_ALL$mu[,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = 2.2, expression("Distance to Prey (mm) "~(delta) ))
mtext(side = 2,cex=cex, line = 2.2, expression("Capture Speed (mm/sec)  " ))

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  ),
                  bquote(All ~ '#' ~ .(ldata_ALL$N)  ) ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL, col=colourLegL)
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)


  
  ## Plot the COVARIANCE ##
  pdf(file= paste(strPlotExportPath,strModelCoVarPDFFilename,sep=""),width=14,height=7,
      title="A statistical model for Covariance of Capture speed to Distance to prey")
  
  layout(matrix(c(1,2),1,2, byrow = TRUE))
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(3.9,4.3,3,1))
  
  plot(dNLb_rho_slow,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,5),
       main=NA,cex=cex, #"Density Inference of Turn-To-Prey Slope ",
       xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
  lines(dLLb_rho_slow,col=colourLegL[2],lwd=3,lty=2)
  lines(dDLb_rho_slow,col=colourLegL[3],lwd=3,lty=3)
  #lines(dALLb_rho,col=colourLegL[4],lwd=3,lty=4)
  mtext(side = 2,cex=cex, line = lineAxis-0.5, expression("Density ") )
  mtext(side = 1,outer=F,cex=cex, line = lineXAxis,adj=0.5 ,expression(paste("Capture speed to prey distance  covariance ",(rho["s"]) ) ))
  mtext(side = 3,cex=cex, line = lineAxis-2, expression("Slow") )
  mtext("A",at="topleft",side=2,col="black",font=2,las=las,line=line,padj=-18,adj=adj,cex=cex)
  
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
                    bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
                    bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
                    #bquote(ALL ~ '#' ~ .(ldata_ALL$N)  ) 
         ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)
  
  
  plot(dNLb_rho_fast,col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,5),
       main=NA,cex=cex, #"Density Inference of Turn-To-Prey Slope ",
       xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
  lines(dLLb_rho_fast,col=colourLegL[2],lwd=3,lty=2)
  lines(dDLb_rho_fast,col=colourLegL[3],lwd=3,lty=3)
  #lines(d
  
  mtext(side = 1,outer=F,cex=cex, line = lineXAxis,adj=0.5 ,expression(paste("Capture speed to prey distance covariance ",(rho["f"]) ) ))
  mtext(side = 2,cex=cex, line = lineAxis-0.5, expression("Density ") )
  mtext(side = 3,cex=cex, line = lineAxis-2, expression("Fast") )
  mtext("B",at="topleft",side=2,col="black",font=2,las=las,line=line,padj=-18,adj=adj,cex=cex)
  
  dev.off()





### ADD DISTANCE TO PREY VARIANCE COMPARISON
plot(dNLb_sigmaD,col=colourLegL[1],xlim=c(0,0.5),lwd=3,lty=1,ylim=c(0,20),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_sigmaD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_sigmaD,col=colourLegL[3],lwd=3,lty=3)
#lines(dALLb_sigmaD,col=colourLegL[4],lwd=3,lty=4)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Variance Prey Distance  ",delta) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )

### PloT CAPT SPEED VARIANCE 

plot(dNLb_sigmaC,col=colourLegL[1],xlim=c(0.0,30),lwd=3,lty=1,ylim=c(0,0.3),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb_sigmaC,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_sigmaC,col=colourLegL[3],lwd=3,lty=3)
#lines(dALLb_sigmaC,col=colourLegL[4],lwd=3,lty=4)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("Variance Capture Speed  ") ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )

dev.off()




#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)

###    Plot Cluster Membership ratios , mean number of HUnt Events in every cluster #####
##Select the Cluster ID that has the fastest capture speed, as being the capture strike
idxCaptClust_LF <- ifelse ( mean(draw_LF$mu[2,2,,1]) > mean(draw_LF$mu[1,2,,1]),2,1 ) 
idxCaptClust_NF <- ifelse ( mean(draw_NF$mu[2,2,,1]) > mean(draw_NF$mu[1,2,,1]),2,1 ) 
idxCaptClust_DF <- ifelse ( mean(draw_DF$mu[2,2,,1]) > mean(draw_DF$mu[1,2,,1]),2,1 ) 


#' Count the numbe of huntevents classed as having a capture strike accoridng to the 2-gaussian fit
getClusterOccupancySamples <- function(draw_G,ntail)
{
  ##get mean membership per Hunt Event 
  vMembership <- vector()
  n <- NROW(draw_G$mu[1,2,,1])
  ntail <- min(n, ntail) ##check for ntail exceeding number of samples
  for (i in 1:ntail)
  {
    ##Get Last nSamples For each Hunt Event Memberhip vector / Count number of events in Fast cluster over each sample
    vMembership[i] <- sum(draw_G$mID[, (n-ntail):n,1][,i])
    
  }

  return(vMembership)  
}

vMembership_LF <- getClusterOccupancySamples(draw_LF,100);
vMembership_NF <- getClusterOccupancySamples(draw_NF,100);
Membership_DF <- getClusterOccupancySamples(draw_DF,100);
layout(matrix(c(1,2,3),3,1, byrow = TRUE))
hist(vMembership_LF/NROW(draw_LF$mID),xlim=c(0,1),main= paste("LF Strike Speed",prettyNum( mean(draw_LF$mu[idxCaptClust_LF,2,,1]),digits=4 )  ) )
hist(vMembership_NF/NROW(draw_NF$mID),xlim=c(0,1),main= paste("NF Strike Speed",prettyNum( mean(draw_NF$mu[idxCaptClust_NF,2,,1]),digits=4 )  ) )
hist(vMembership_DF/NROW(draw_DF$mID),xlim=c(0,1),main= paste("DF Strike Speed",prettyNum( mean(draw_DF$mu[idxCaptClust_DF,2,,1]),digits=4 )  ) )

layout(matrix(c(1,2,3),3,1, byrow = TRUE))
hist(draw_LF$mStrikeCount/NROW(draw_LF$mID),xlim=c(0,1),main= paste("LF Strike Speed",prettyNum( mean(draw_LF$mu[idxCaptClust_LF,2,,1]),digits=4 )  ) )
hist(draw_NF$mStrikeCount/NROW(draw_NF$mID),xlim=c(0,1),main= paste("NF Strike Speed",prettyNum( mean(draw_NF$mu[idxCaptClust_NF,2,,1]),digits=4 )  ) )
hist(draw_DF$mStrikeCount/NROW(draw_DF$mID),xlim=c(0,1),main= paste("DF Strike Speed",prettyNum( mean(draw_DF$mu[idxCaptClust_DF,2,,1]),digits=4 )  ) )


## Show Whther LF does proportionally more strike swimms than the other two
pdf(file= paste(strPlotExportPath,strClusterOccupancyPDFFileName,sep=""))
boxplot(vMembership_NF/NROW(draw_NF$mID),vMembership_LF/NROW(draw_LF$mID),vMembership_DF/NROW(draw_DF$mID),ylim=c(0,1.0),
        main="Estimated ratio of high speed captures ",names=c("NF","LF","DF") , col=colourH )
dev.off()
#table(draw_LF$mID[,,1])[idxCaptClust_LF]/table(draw_LF$mID[,,1])[1]
#table(draw_NF$mID[,,1])[idxCaptClust_LF]/table(draw_NF$mID[,,1])[1]
#table(draw_DF$mID[,,1])[idxCaptClust_LF]/table(draw_DF$mID[,,1])[1]
###




#### Plot Prey Location  ###########
## The Original list if the lFirstBout data from runHuntepisode analysis
source("plotTrackScatterAndDensities.r")
#plotCaptureBoutPreyPositions

##########################################################
########## GAPE TIMING - Capture Strike Duration #########
##########################################################


### Make Gape - Timing Estimates for Fast Captures
##These can be compared with Estimated time from frame durations of tracked hunt episodes- 
gape_timing_comp_LF <- datDistanceVsStrikeSpeed_LL$DistanceToPrey[lClustScore_LF$pchL == 16]/datDistanceVsStrikeSpeed_LL$CaptureSpeed[lClustScore_LF$pchL == 16]
gape_timing_comp_NF <- datDistanceVsStrikeSpeed_NL$DistanceToPrey[lClustScore_NF$pchL == 16]/datDistanceVsStrikeSpeed_NL$CaptureSpeed[lClustScore_NF$pchL == 16]
gape_timing_comp_DF <- datDistanceVsStrikeSpeed_DL$DistanceToPrey[lClustScore_DF$pchL == 16]/datDistanceVsStrikeSpeed_DL$CaptureSpeed[lClustScore_DF$pchL == 16]

dens_timing_LF_fast <- density(gape_timing_comp_LF,bw=0.002)
dens_timing_NF_fast <- density(gape_timing_comp_NF,bw=0.002)
dens_timing_DF_fast <- density(gape_timing_comp_DF,bw=0.002)

### Compare To Bout Duration ###

##Select the fast Capture/Strikes And Check duration 
dens_captiming_LF_fast <- density(datDistanceVsStrikeSpeed_LL$CaptureDuration[lClustScore_LF$pchL == 16],bw=0.1)
dens_captiming_NF_fast <- density(datDistanceVsStrikeSpeed_NL$CaptureDuration[lClustScore_NF$pchL == 16],bw=0.1)
dens_captiming_DF_fast <- density(datDistanceVsStrikeSpeed_DL$CaptureDuration[lClustScore_DF$pchL == 16],bw=0.1)

###
layout(matrix(c(1,2),2,1, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
#par(mar = c(5,4.5,3,1))
par(mar = c(4.5,4.7,2,1))


plot(dens_timing_LF_fast,col=colourLegL[2],lwd=2,main="Capture peak speed/distance to prey")
lines(dens_timing_NF_fast,col=colourLegL[1],lwd=2)
lines(dens_timing_DF_fast,col=colourLegL[3],lwd=2)

plot(dens_captiming_LF_fast,col=colourLegL[2],lwd=3,main="Capture swim Duration")
lines(dens_captiming_NF_fast,col=colourLegL[1],lwd=3)
lines(dens_captiming_DF_fast,col=colourLegL[3],lwd=3)


plot(datDistanceVsStrikeSpeed_LL$CaptureDuration,datDistanceVsStrikeSpeed_LL$CaptureSpeed)



plot( lFirstBoutPoints$LL[,"PeakSpeedDistance"][lClustScore_LF$pchL == 16],lFirstBoutPoints$LL[,"CaptureSpeed"][lClustScore_LF$pchL == 16],xlim=c(0,0.5),ylim=c(0,70),ylab="Speed",xlab="Dist Travelled to Peak Speed")
plot( lFirstBoutPoints$NL[,"PeakSpeedDistance"][lClustScore_NF$pchL == 16],lFirstBoutPoints$NL[,"CaptureSpeed"][lClustScore_NF$pchL == 16],xlim=c(0,0.5),ylim=c(0,70),ylab="Speed",xlab="Dist Travelled to Peak Speed")
plot( lFirstBoutPoints$DL[,"PeakSpeedDistance"][lClustScore_DF$pchL == 16],lFirstBoutPoints$DL[,"CaptureSpeed"][lClustScore_DF$pchL == 16],xlim=c(0,0.5),ylim=c(0,70),ylab="Speed",xlab="Dist Travelled to Peak Speed")


## I ve extracted MotionBoutData on Time to peak speed and distance travelled --
##Papers suggest peak speed hits close to prey - Here tho - SpeedDistance is displacent of as sum(FishSpeed)
plot( lFirstBoutPoints$LL[,"PeakSpeedDistance"][lClustScore_LF$pchL == 16],lFirstBoutPoints$LL[,"DistanceToPrey"][lClustScore_LF$pchL == 16],xlim=c(0,0.5),ylim=c(0,0.5),ylab="Dist To Prey",xlab="Dist Travelled to Peak Speed")
plot( lFirstBoutPoints$NL[,"PeakSpeedDistance"][lClustScore_NF$pchL == 16],lFirstBoutPoints$NL[,"DistanceToPrey"][lClustScore_NF$pchL == 16],xlim=c(0,0.5),ylim=c(0,0.5),ylab="Dist To Prey",xlab="Dist Travelled to Peak Speed")
plot( lFirstBoutPoints$DL[,"PeakSpeedDistance"][lClustScore_DF$pchL == 16],lFirstBoutPoints$DL[,"DistanceToPrey"][lClustScore_DF$pchL == 16],xlim=c(0,0.5),ylim=c(0,0.5),ylab="Dist To Prey",xlab="Dist Travelled to Peak Speed")

cor( lFirstBoutPoints$LL[,"PeakSpeedDistance"][lClustScore_LF$pchL == 16],lFirstBoutPoints$LL[,"DistanceToPrey"][lClustScore_LF$pchL == 16])
cor( lFirstBoutPoints$NL[,"PeakSpeedDistance"][lClustScore_NF$pchL == 16],lFirstBoutPoints$NL[,"DistanceToPrey"][lClustScore_NF$pchL == 16])
cor( lFirstBoutPoints$DL[,"PeakSpeedDistance"][lClustScore_DF$pchL == 16],lFirstBoutPoints$DL[,"DistanceToPrey"][lClustScore_DF$pchL == 16])


plot( lFirstBoutPoints$LL[,"CaptureSpeed"],lFirstBoutPoints$LL[,"NFramesToPeakSpeed"])
plot( lFirstBoutPoints$NL[,"CaptureSpeed"],lFirstBoutPoints$NL[,"NFramesToPeakSpeed"])


###Time Until Min Distance to Prey
plot(density((lFirstBoutPoints$LL[,"ColisionFrame"]-lFirstBoutPoints$LL[,"CaptureBoutStartFrame"])/G_APPROXFPS ),xlim=c(0,1),col=colourLegL[2],main="All captures")
lines(density((lFirstBoutPoints$NL[,"ColisionFrame"]-lFirstBoutPoints$NL[,"CaptureBoutStartFrame"])/G_APPROXFPS,na.rm=T),col=colourLegL[1])
lines(density((lFirstBoutPoints$DL[,"ColisionFrame"]-lFirstBoutPoints$DL[,"CaptureBoutStartFrame"])/G_APPROXFPS,na.rm=T),col=colourLegL[3])

###Time Until Min Distance to Prey
plot(density((lFirstBoutPoints$LL[,"ColisionFrame"][lClustScore_LF$pchL == 16]-lFirstBoutPoints$LL[,"CaptureBoutStartFrame"][lClustScore_LF$pchL == 16])/G_APPROXFPS ),xlim=c(0,1),col=colourLegL[2],main="fast")
lines(density((lFirstBoutPoints$NL[,"ColisionFrame"][lClustScore_NF$pchL == 16]-lFirstBoutPoints$NL[,"CaptureBoutStartFrame"][lClustScore_NF$pchL == 16])/G_APPROXFPS,na.rm=T),col=colourLegL[1])
lines(density((lFirstBoutPoints$DL[,"ColisionFrame"][lClustScore_DF$pchL == 16]-lFirstBoutPoints$DL[,"CaptureBoutStartFrame"][lClustScore_DF$pchL == 16])/G_APPROXFPS,na.rm=T),col=colourLegL[3])

##Fast CLuster Data points On Time to Hit Prey on capture swim (Time to hit prey - should not covary with distance for a fixed gape timing)
timeToHit_LF_fast <- (lFirstBoutPoints$LL[,"ColisionFrame"][lClustScore_LF$pchL == 16]-lFirstBoutPoints$LL[,"CaptureBoutStartFrame"][lClustScore_LF$pchL == 16])/G_APPROXFPS 
timeToHit_NF_fast <-(lFirstBoutPoints$NL[,"ColisionFrame"][lClustScore_NF$pchL == 16]-lFirstBoutPoints$NL[,"CaptureBoutStartFrame"][lClustScore_NF$pchL == 16])/G_APPROXFPS
timeToHit_DF_fast <- (lFirstBoutPoints$DL[,"ColisionFrame"][lClustScore_DF$pchL == 16]-lFirstBoutPoints$DL[,"CaptureBoutStartFrame"][lClustScore_DF$pchL == 16])/G_APPROXFPS


### Time to Hit Prey Should be invariant to Distance from prey in FAST captures - Gape Timing
###
layout(matrix(c(1,2,3),3,1, byrow = TRUE))
plot( lFirstBoutPoints$LL[,"DistanceToPrey"][lClustScore_LF$pchL == 16],timeToHit_LF_fast,xlim=c(0,0.5),pch=pchL[1],ylim=c(0,0.6))
plot( lFirstBoutPoints$NL[,"DistanceToPrey"][lClustScore_NF$pchL == 16],timeToHit_NF_fast,pch=pchL[2],xlim=c(0,0.5),ylim=c(0,0.6))
plot( lFirstBoutPoints$DL[,"DistanceToPrey"][lClustScore_DF$pchL == 16],timeToHit_DF_fast,pch=pchL[3],xlim=c(0,0.5),ylim=c(0,0.6))

sd( timeToHit_LF_fast,na.rm=T )
sd( timeToHit_NF_fast,na.rm=T ) 
sd( timeToHit_DF_fast,na.rm=T ) 
