##3-09-2018
###  Estimates the hidden function of Turn Vs Bearing To Prey  - 
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


modelGPV1="model {
  # Likelihood

for(i in 1:N){
  turn[i] ~ dnorm(lambda[i],eps)
  #n[i] ~ dpois(lambda[i])
}

eps~dexp(10)
lambda ~ dmnorm(Mu, Sigma.inv)
#n ~ dmnorm(Mu, Sigma.inv)
Sigma.inv <- inverse(Sigma)

# Set up mean and covariance matrix
for(i in 1:N) {
  Mu[i] <- alpha
  Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)

  for(j in (i+1):N) {
    Sigma[i,j] <- pow(tau,2) * exp( - rho * pow(bearing[i] - bearing[j], 2) )
    Sigma[j,i] <- Sigma[i,j]
  }
}

alpha ~ dnorm(0,1e-4)T(0,) 
#tau ~ dnorm(tauRange,1e-1)T(0,)
#rho = rhoMax

tau0 ~ dgamma(tauRange,0.2) 
tau  ~ dgamma(tauRange,0.2) 
rho ~ dunif(0,rhoMax)

}"

#+ beta[3]*turn[i]

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
datMotionBoutCombinedAll <-  data.frame( do.call(rbind,lMotionBoutDat ) )
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
remove(lMotionBoutDat)
lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData.rds",sep="") ) #Processed Registry on which we add )


for (gp in strGroupID)
{
  groupID <- which(levels(datTrackedEventsRegister$groupID) == gp)
  
  datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm <- as.numeric(datMotionBoutCombinedAll$vMotionBoutDistanceToPrey_mm)
  datMotionBoutCombined <-datMotionBoutCombinedAll[datMotionBoutCombinedAll$groupID == as.numeric(groupID), ] #Select Group
  
  datMotionBoutCombined$boutRank <- as.numeric(datMotionBoutCombined$boutRank)
  ## Punctuate 1st Turn To Prey
  #lFirstBoutPoints[[gp]] <- cbind(OnSetAngleToPrey = datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey,
  #                            Turn= datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey - datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1,]$OffSetAngleToPrey
  #                            , RegistarIdx=datMotionBoutCombined[datMotionBoutCombined$turnSeq == 1 & datMotionBoutCombined$boutSeq == 1 ,]$RegistarIdx)
  lFirstBoutPoints[[gp]] <- cbind(OnSetAngleToPrey = datMotionBoutCombined[datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey,
                                  Turn= datMotionBoutCombined[ datMotionBoutCombined$boutSeq == 1 ,]$OnSetAngleToPrey - datMotionBoutCombined[ datMotionBoutCombined$boutSeq == 1,]$OffSetAngleToPrey
                                  , RegistarIdx=datMotionBoutCombined[ datMotionBoutCombined$boutSeq == 1 ,]$RegistarIdx)
}
  

##Add The Empty Test Conditions
#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
#datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]

## Get Event Counts Within Range ##
datTurnVsPreyLL <- cbind(lFirstBoutPoints$LL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$LL[,"Turn"]) )
datTurnVsPreyLL <- datTurnVsPreyLL[!is.na(datTurnVsPreyLL[,1]),]


datTurnVsPreyNL <- cbind(lFirstBoutPoints$NL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$NL[,"Turn"]) )
datTurnVsPreyNL <- datTurnVsPreyNL[!is.na(datTurnVsPreyNL[,1]),]

datTurnVsPreyDL <- cbind(lFirstBoutPoints$DL[,"OnSetAngleToPrey"] , as.numeric(lFirstBoutPoints$DL[,"Turn"]) )
datTurnVsPreyDL <- datTurnVsPreyDL[!is.na(datTurnVsPreyDL[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
colourH <- c(rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.01,0.01,0.9,0.8),rgb(0.00,0.00,0.0,1.0))
colourP <- c(rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.01,0.01,0.8,0.5),rgb(0.00,0.00,0.0,1.0))
colourR <- c(rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.01,0.01,0.9,0.4),rgb(0.00,0.00,0.0,1.0))
##Thse RC params Work Well to Smooth LF And NF
tauRangeA =100000 #10000
rhoMaxA = 1000
Noise = 1 ##The Gaussian Noise Term

burn_in=10;
steps=1000;
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


dataLL=list(turn=turnsLL,bearing=turnsLL,N=nDatLL,tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataNL=list(turn=turnsNL,bearing=turnsNL,N=nDatNL,tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
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

strPlotName <-  paste(strPlotExportPath,"/stat_TurnVsBearing_GPEstimate-tauMax",tauRangeA,"Rho",rhoMaxA,".pdf",sep="")
pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 
myplot_res(1000)
dev.off()



hist(drawLL$beta[2,,1],breaks=seq(0,2,length=100),col=colourH[1],xlim=c(-1,15),
     #xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
     xlab=expression(paste("Turn to Prey Bearing ",lambda)),main=paste("Slope ") )
hist(drawNL2$q[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[2],xlim=c(5,15))
hist(drawDL2$q[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[3],xlim=c(5,15))

myplot_res(1000)

X11()
hist(drawLL$beta[2,,1],breaks=seq(0.9,1.1,length=100),col=colourH[1],xlim=c(0.9,1.1),
     #xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
     xlab=expression(paste("Turn to Prey Bearing ",lambda)),main=paste("Slope ") )

hist(drawNL$beta[2,,1],breaks=seq(0.9,1.1,length=100),col=colourH[2],xlim=c(0.9,1.1),add=T  )


X11()
hist(drawLL$beta[2,,1])

ind = 100
muLL=mean(drawLL$beta[,(steps-ind):steps,1][2,])
muNL=mean(drawNL$beta[,(steps-ind):steps,1][2,])
sig=mean(drawLL$sigma[,(steps-ind):steps,1])

##Plot Densities Summary
sampLL <- coda.samples(mLL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL,main="LL")
X11()
plot(sampNL,main="NL")
X11()
plot(sampDL,main="DL")