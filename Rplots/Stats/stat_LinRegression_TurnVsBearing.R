##3-09-2018
###  Estimates the hidden function of Turn Vs Bearing To Prey - Linear Regression gives pdf of slope param.  - 
## Tried both a Non-parametric Gaussian Process with Bayesian Inference (Failed) But Also a simple Linear Model
## Requires the lFirstBoutPoints list of dataframes - which is constucted in 

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

#source("DataLabelling/labelHuntEvents_lib.r")
### GP Process Estimation Of Hunt Rate Vs Prey Density Using Bayesian Inference Model
myplot_res<- function(ind,qq=0.05){
  
  
  xplotLim <- c(-100,100)
  yplotLim <- c(-100,100)
  plot(bearingLL,turnsLL,col=colourP[1],
       main = "GP Regression Of Turn Vs Bearing To Prey ",
       ylab="Turn To Prey",
       xlab="Bearing To Prey",
       xlim = xplotLim,
       ylim = yplotLim,
       pch=16,
       sub=paste("GP tau:",format(mean(drawLL$tau),digits=4 ),
                 "rho:",format(mean(drawLL$rho),digits=4 ) )  
       )
  
  legend("topright",legend = c(paste("LL #",nDatLL),paste("NL #",nDatNL),paste("DL #",nDatDL)),fill=colourH)
  
  points(bearingNL,turnsNL,col=colourP[2],pch=16,xlim = xplotLim)
  points(bearingDL,turnsDL,col=colourP[3],pch=16,xlim = xplotLim)
  
  muLL=apply(drawLL$lambda[,(steps-ind):steps,1],1,mean)
  muNL=apply(drawNL$lambda[,(steps-ind):steps,1],1,mean)
  muDL=apply(drawDL$lambda[,(steps-ind):steps,1],1,mean)
  
  lines(bearingLL,muLL,col=colourH[1],lwd=4,xlim = xplotLim)
  lines(bearingNL,muNL,col=colourH[2],lwd=4,xlim = xplotLim)
  lines(bearingDL,muDL,col=colourH[3],lwd=4,xlim = xplotLim)
  
  band=apply(drawLL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(bearingLL,rev(bearingLL)),c(band[1,],rev(band[2,])),col=colourR[1])
  
  band=apply(drawNL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(bearingNL,rev(bearingNL)),c(band[1,],rev(band[2,])),col=colourR[2])
  
  band=apply(drawDL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(bearingDL,rev(bearingDL)),c(band[1,],rev(band[2,])),col=colourR[3])
  
}

## This Is the Linear Regression Model Used
modelLin <- "model{

  # Likelihood
  for(i in 1:N){
    turn[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- beta[1] + beta[2]*bearing[i] 
  }

  # Prior for beta
  for(j in 1:2){
    beta[j] ~ dnorm(0,0.0001)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  sigma     <- 1/sqrt(inv.var)

}"

####Select Subset Of Data To Analyse

datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
remove(lMotionBoutDat)
lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData.rds",sep="") ) #Processed Registry on which we add )
remove(lFirstBoutPoints) ##Load From File
lFirstBoutPoints <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData",".rds",sep="") ) #Processed Registry on which we add )

# datMotionBoutCombinedAll <-  data.frame( do.call(rbind,lMotionBoutDat ) )
# lFirstBoutPoints <- list() ##Add Dataframes Of 1st bout Turns for Each Group
# strGroupID <- levels(datTrackedEventsRegister$groupID)
# for (gp in strGroupID)
# {
#   groupID <- which(levels(datTrackedEventsRegister$groupID) == gp)
#   
#   datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm <- as.numeric(datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm)
#   datMotionBoutCombined <-datMotionBoutCombinedAll[datMotionBoutCombinedAll$groupID == as.numeric(groupID), ] #Select Group
#   
#   datMotionBoutCombined$boutRank <- as.numeric(datMotionBoutCombined$boutRank)
#   datMotionBoutTurnToPrey <- datMotionBoutCombined[abs(datMotionBoutCombined$OnSetAngleToPrey) >= abs(datMotionBoutCombined$OffSetAngleToPrey) , ]
#   datMotionBoutTurnToPrey <- datMotionBoutTurnToPrey[!is.na(datMotionBoutTurnToPrey$RegistarIdx),]
#   ## Punctuate 1st Turn To Prey
#   #lFirstBoutPoints[[gp]] <- cbind(OnSetAngleToPrey = datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey,
#   #                            Turn= datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey - datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1,]$OffSetAngleToPrey
#   #                            , RegistarIdx=datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$RegistarIdx)
#   lFirstBoutPoints[[gp]] <- cbind(OnSetAngleToPrey = datMotionBoutTurnToPrey[datMotionBoutTurnToPrey$turnSeq == 1 ,]$OnSetAngleToPrey,
#                                   Turn= datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$OnSetAngleToPrey - datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1,]$OffSetAngleToPrey
#                                   , RegistarIdx=datMotionBoutTurnToPrey[ datMotionBoutTurnToPrey$turnSeq == 1 ,]$RegistarIdx)
#   
#   
# }
#   

##Add The Empty Test Conditions
#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
#datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]

## Get Event Counts Within Range ##
datTurnVsPreyLL <- cbind(lFirstBoutPoints$LL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$LL[,"Turn"]),lFirstBoutPoints$LL[,"RegistarIdx"] )
datTurnVsPreyLL <- datTurnVsPreyLL[!is.na(datTurnVsPreyLL[,1]),]


datTurnVsPreyNL <- cbind(lFirstBoutPoints$NL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$NL[,"Turn"]),lFirstBoutPoints$NL[,"RegistarIdx"] )
datTurnVsPreyNL <- datTurnVsPreyNL[!is.na(datTurnVsPreyNL[,1]),]

datTurnVsPreyDL <- cbind(lFirstBoutPoints$DL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$DL[,"Turn"]),lFirstBoutPoints$DL[,"RegistarIdx"] )
datTurnVsPreyDL <- datTurnVsPreyDL[!is.na(datTurnVsPreyDL[,1]),]

##Outlier datTurnVsPreyDL[13,] <- NA

##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points]
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
pchL <- c(16,2,4)
#
#Thse RC params Work Well to Smooth LF And NF
tauRangeA =100000 #10000
rhoMaxA = 1000
Noise = 1 ##The Gaussian Noise Term

burn_in=10;
steps=100000;
thin=1;


##Larva Event Counts Slice
nDatLL <- NROW(datTurnVsPreyLL)
nDatNL <- NROW(datTurnVsPreyNL)
nDatDL <- NROW(datTurnVsPreyDL)

##Order Data in Bearing Sequence 
bearingLL=datTurnVsPreyLL[,1]
ordLL=order(bearingLL)
bearingLL=bearingLL[ordLL]
turnsLL=datTurnVsPreyLL[,2][ordLL]

##Order Data in Bearing Sequence 
bearingNL=datTurnVsPreyNL[,1]
ordNL=order(bearingNL)
bearingNL=bearingNL[ordNL]
turnsNL=datTurnVsPreyNL[,2][ordNL]


##Order Data in Bearing Sequence 
bearingDL=datTurnVsPreyDL[,1]
ordDL=order(bearingDL)
bearingDL=bearingDL[ordDL]
turnsDL=datTurnVsPreyDL[,2][ordDL]


dataLL=list(turn=turnsLL,bearing=bearingLL,N=nDatLL,tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataNL=list(turn=turnsNL,bearing=bearingNL,N=nDatNL,tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataDL=list(turn=turnsDL,bearing=bearingDL,N=nDatDL,tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);

varnames=c("tau","rho","alpha","lambda")
varnames=c("beta","sigma")

library(rjags)
fileConn=file("model.tmp")
#writeLines(modelGPV1,fileConn);
writeLines(modelLin,fileConn);
close(fileConn)

mLL=jags.model(file="model.tmp",data=dataLL);
mNL=jags.model(file="model.tmp",data=dataNL);
mDL=jags.model(file="model.tmp",data=dataDL);
#update(mLL,burn_in);update(mNL,burn_in);update(mDL,burn_in)


drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)



## Plot ### 
ind = 10000 ## Number of last sampled values
## Save the Mean Slope and intercept
##quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2]
muLLa=mean(drawLL$beta[,(steps-ind):steps,1][1,]) 
muLLb=mean(drawLL$beta[,(steps-ind):steps,1][2,])
muNLa=mean(drawNL$beta[,(steps-ind):steps,1][1,])
muNLb=mean(drawNL$beta[,(steps-ind):steps,1][2,])
muDLa=mean(drawDL$beta[,(steps-ind):steps,1][1,])
muDLb=mean(drawDL$beta[,(steps-ind):steps,1][2,])
sig=mean(drawLL$sigma[,(steps-ind):steps,1])
###Plot Density of Slope
dLLb<-density(drawLL$beta[,(steps-ind):steps,1][2,])
dNLb<-density(drawNL$beta[,(steps-ind):steps,1][2,])
dDLb<-density(drawDL$beta[,(steps-ind):steps,1][2,])


pdf(file= paste(strPlotExportPath,"/stat/stat_densityolinregressionslope.pdf",sep=""))
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),
     main="Density Inference of Turn-To-Prey Slope ",
     xlab=expression(paste("slope ",gamma) ) )
lines(dLLb,col=colourH[2],xlim=c(0.5,1.2),lwd=3,lty=2)
lines(dNLb,col=colourH[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

### PLot Scatter with regression lines with Conf intervals##
#X11()

pdf(file= paste(strPlotExportPath,"/stat/stat_TurnToPrey_LinearRegression.pdf",sep=""))
plot(lFirstBoutPoints[["DL"]][,1], lFirstBoutPoints[["DL"]][,2],
     main=paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[1] ,pch=pchL[1]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=1,lty=2);
#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(lFirstBoutPoints[["DL"]][,2] ~ lFirstBoutPoints[["DL"]][,1]),col=colourH[4],lwd=1.0) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(lFirstBoutPoints[["LL"]][,1], lFirstBoutPoints[["LL"]][,2],pch=pchL[2],col=colourP[2])
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(lFirstBoutPoints[["LL"]][,2] ~ lFirstBoutPoints[["LL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points(lFirstBoutPoints[["NL"]][,1], lFirstBoutPoints[["NL"]][,2],pch=pchL[3],col=colourP[3])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(lFirstBoutPoints[["NL"]][,2] ~ lFirstBoutPoints[["NL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muNLa,b=muNLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       , pch=pchL,col=colourL)

dev.off()



### Model Evidence ### 
## Compare Likelyhoods between models for undershooting (linear slope fit) 
## Establish if DryFed Data Belong to NF rather then LF
## 
## This was given by Giovanni , based on model evidence formula (see wikipedia Bayesian linear regression)
getParams <- function(data,a0=1,b0=1,sigma0=1){
  n=nrow(data)
  y=data[,2]
  X=cbind(1,data[,1])
  Lambda0 = diag(sigma0,2)
  Lambda  = t(X)%*%X+Lambda0 
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
  mu0 = c(1,0)               
  mu  = solve(t(X)%*%X+Lambda0)%*%(t(X)%*%X%*%beta_hat+Lambda0%*%mu0)
  a=a0+n/2                           
  b=b0+0.5*(t(y)%*%y+t(mu0)%*%Lambda0%*%mu0-t(mu)%*%Lambda%*%mu)
  
  return(list(n=n,a=a,b=b,mu=mu,lambda=Lambda))  
}

##Marginal Likelyhood 
MarginalLikelihood <- function(MLParams,a0,b0)
{
  return (1/(2*pi)^(MLParams$n/2))* sqrt( det(diag(sigma0,2))/det( MLParams$lambda))*b0/MLParams$b*gamma(MLParams$a)/gamma(a0)
}

b0=1
a0=1
MLparamsLL <- getParams( cbind(dataLL$turn,dataLL$bearing),a0,b0 )
MLparamsDL <- getParams( cbind(dataDL$turn,dataDL$bearing),a0,b0 )
MLparamsNL <- getParams( cbind(dataNL$turn,dataNL$bearing),a0,b0 )

dataNLDL <- rbind(cbind(dataNL$turn,dataNL$bearing),cbind(dataDL$turn,dataDL$bearing))
MLparamsNLDL <- getParams( dataNLDL,a0,b0 )

dataDLLL <- rbind(cbind(dataDL$turn,dataDL$bearing),cbind(dataLL$turn,dataLL$bearing))
MLparamsDLLL <- getParams( dataDLLL,a0,b0 )


ML_LL <- MarginalLikelihood(MLparamsLL,a0,b0)
ML_DL <- MarginalLikelihood(MLparamsDL,a0,b0)
ML_NL <- MarginalLikelihood(MLparamsNL,a0,b0)
ML_NLDL <- MarginalLikelihood(MLparamsNLDL,a0,b0)
ML_DLLL <- MarginalLikelihood(MLparamsDLLL,a0,b0)


##Now Compare ##
# A value of K > 1 means that M1 is more strongly supported by the data under consideration than M2.
ML_DL*ML_NL/(ML_NLDL)

##
ML_DL*ML_LL/(ML_DLLL)

mean(dataDL$turn/dataDL$bearing)
mean(dataNL$turn/dataNL$bearing)

#strPlotName <-  paste(strPlotExportPath,"/stat_TurnVsBearing_GPEstimate-tauMax",tauRangeA,"Rho",rhoMaxA,".pdf",sep="")
#pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 
#myplot_res(1000)
#dev.off()

# 
# X11()
# hist(drawLL$beta[1,,1],breaks=seq(0.91,1.15,length=00),col=colourH[1],
#      #xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
#      xlab=expression(paste("Turn to Prey Bearing ",lambda)),main=paste("Slope ") )
# hist(drawNL$beta[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[2],xlim=c(5,15))
# hist(drawDL$beta[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[3],xlim=c(5,15))
# 
# myplot_res(1000)
# 
# X11()
# hist(drawLL$beta[2,,1],breaks=seq(0.9,1.1,length=100),col=colourH[1],xlim=c(0.9,1.1),
#      #xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
#      xlab=expression(paste("Turn to Prey Bearing ",lambda)),main=paste("Slope ") )
# 
# 
# hist(drawNL$beta[2,,1],breaks=seq(0.9,1.1,length=100),col=colourH[2],xlim=c(0.9,1.1),add=T  )
# 
# 
# X11()
# hist(drawLL$beta[2,,1])


##Plot Densities Summary
#sampLL <- coda.samples(mLL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
#sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
#sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
#X11()
#plot(sampLL)
#X11()
#plot(sampNL)
#X11()
#plot(sampDL,main="DL")



# 
# modelGPV1="model {
#   # Likelihood
# 
# for(i in 1:N){
#   turn[i] ~ dnorm(lambda[i],eps)
#   #n[i] ~ dpois(lambda[i])
# }
# 
# eps~dexp(10)
# lambda ~ dmnorm(Mu, Sigma.inv)
# #n ~ dmnorm(Mu, Sigma.inv)
# Sigma.inv <- inverse(Sigma)
# 
# # Set up mean and covariance matrix
# for(i in 1:N) {
#   Mu[i] <- alpha
#   Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
# 
#   for(j in (i+1):N) {
#     Sigma[i,j] <- pow(tau,2) * exp( - rho * pow(bearing[i] - bearing[j], 2) )
#     Sigma[j,i] <- Sigma[i,j]
#   }
# }
# 
# alpha ~ dnorm(0,1e-4)T(0,) 
# #tau ~ dnorm(tauRange,1e-1)T(0,)
# #rho = rhoMax
# 
# tau0 ~ dgamma(tauRange,0.2) 
# tau  ~ dgamma(tauRange,0.2) 
# rho ~ dunif(0,rhoMax)
# 
# }"

#+ beta[3]*turn[i]


#modelLin2 <- "model {
#	for (i in 1:N){
# 		turn[i] ~ dnorm(turn.hat[i], tau)
# 		turn.hat[i] <- beta[1] + beta[2] * bearing[i]
# 	}
# 	beta[1] ~ dnorm(0, .0001)
# 	beta[2] ~ dnorm(0, .0001)
# 	tau <- pow(sigma, -2)
# 	sigma ~ dunif(0, 100)
# }"
