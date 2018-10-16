##15-10-2018
### Model how Eye Angle reports distance from prey 
## Used to compare differences in distance estimation between rearing groups 


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

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

modelLin2 <- "model {
	for (i in 1:N){
		turn[i] ~ dnorm(turn.hat[i], tau)
		turn.hat[i] <- beta[1] + beta[2] * bearing[i]
	}
	beta[1] ~ dnorm(0, .0001)
	beta[2] ~ dnorm(0, .0001)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}"

##The Eye Angle Vs Distance Model
## Regression of an exponential Function for Eye Distance
modelExp  <- "model{
  phi_0 ~ dnorm(10,2) # Idle Eye Position
  phi_max ~ dnorm(35,5) # Max Eye Vergence Angle
  lambda ~ dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
  u1 ~ dunif(0, 4) ## End Hunt Distance - Close to prey
  u0 ~ dunif(u1, 5) ##Start Hunt Distance -Far 
  
  # Likelihood
  for(i in 1:N){
    ##Make indicator if hunt event is within sampled Range 
    #if (u1 < distP[i]  & distP[i] < u0) {
    s[i] <- step(u1 - distP[i])*step(distP[i] - u0) 

    phi_hat[i] <- phi_0 + s[i] * phi_max* (1-exp(-lambda*(distMax[i] - distP[i] ) )) 
    phi[i] ~ dnorm(phi_hat[i],sigma[s[i]+1]) ##choose sigma 

  }

  # Prior Sigma On Eye Angle when  In Or Out of hunt region 
  for(j in 1:2){
    sigma[j] ~ dgamma(0.01, 0.01) ##Draw 
  }

}"

####Select Subset Of Data To Analyse

strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 

lEyeMotionDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData.rds",sep="") ) #Processed Registry on which we add )

datEyeVsPreyCombinedAll <-  data.frame( do.call(rbind,lEyeMotionDat ) )

strGroupID <- levels(datTrackedEventsRegister$groupID)


##Add The Empty Test Conditions
#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
#datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]

## Get Event Counts Within Range ##
datLEyePointsLL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$LEyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$DistToPrey),
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$DistToPreyInit ))
datREyePointsLL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$REyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$DistToPrey),
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),]$DistToPreyInit ))

datLEyePointsLL <- datLEyePointsLL[!is.na(datLEyePointsLL[,2]),]
datREyePointsLL <- datREyePointsLL[!is.na(datLEyePointsLL[,2]),]


datLEyePointsNL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL"),]$LEyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL"),]$DistToPrey) )
datREyePointsNL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL"),]$REyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL"),]$DistToPrey) )
datREyePointsNL <- datREyePointsNL[!is.na(datREyePointsNL[,2]),]

datLEyePointsDL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL"),]$LEyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL"),]$DistToPrey) )
datREyePointsDL <- cbind(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL"),]$REyeAngle,
                         as.numeric(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL"),]$DistToPrey) )


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points]
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
pchL <- c(16,2,4)
#
#Thse RC params Work Well to Smooth LF And NF
burn_in=10;
steps=5000;
thin=1;


##Larva Event Counts Slice
nDatLL <- NROW(datLEyePointsLL)
nDatNL <- NROW(datLEyePointsNL)
nDatDL <- NROW(datLEyePointsDL)

##Test limit data

vsamples <- sample (nDatLL,size=2000)
dataLL=list(phi=datLEyePointsLL[vsamples,1],distP=datLEyePointsLL[vsamples,2],N=NROW(vsamples),distMax=datLEyePointsLL[vsamples,3] );
dataNL=list(phi=datLEyePointsNL[,1],distP=datLEyePointsNL[,2],N=nDatNL);
dataDL=list(phi=datLEyePointsDL[,1],distP=datLEyePointsDL[,2],N=nDatDL);


varnames=c("u0","u1","phi_0","phi_max","lambda","sigma")

library(rjags)
fileConn=file("model.tmp")
#writeLines(modelGPV1,fileConn);
writeLines(modelExp,fileConn);
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


X11()
hist(drawLL$lambda[1,,1],breaks=100,col=colourH[1],
     xlab=paste("Turn to Prey Bearing "),main=paste("Slope ") )

## Plot the infered function
X11()
plot((seq(0,5,by=0.01)), mean(drawLL$phi_max )*(1-exp(- mean(drawLL$lambda)*(5-seq(0,5,by=0.01)) )) 
      ,type="l")

X11()
hist(drawLL$sigma[2,,1],breaks=100,col=colourH[1],
     xlab=paste(""),main=paste("During hunt Sigma  ") )

X11()
hist(drawLL$sigma[1,,1],breaks=100,col=colourH[1],
     xlab=paste(" "),main=paste("Outside hunt Sigma ") )


X11()
hist(drawLL$phi_max[1,,1])

X11()
hist(drawLL$phi_0[1,,1])

X11()
hist(drawLL$u0[1,,1],breaks=100)


X11()
hist(drawLL$u1[1,,1],breaks=100)

hist(drawNL$beta[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[2],xlim=c(5,15))
hist(drawDL$beta[1,,1],breaks=seq(0,30,length=200),add=T,col=colourH[3],xlim=c(5,15))


X11()
hist(drawLL$u0[1,,1],breaks=seq(0.9,1.1,length=100),col=colourH[1],xlim=c(0.9,1.1),
     #xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
     xlab=expression(paste("Turn to Prey Bearing ",lambda)),main=paste("Slope ") )


hist(drawNL$beta[2,,1],breaks=seq(0.9,1.1,length=100),col=colourH[2],xlim=c(0.9,1.1),add=T  )


X11()
hist(drawLL$beta[2,,1])

ind = 100
##Save the Mean Slope and intercept
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
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),main="Density Inference of Turn-To-Prey Slope ")
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





##Plot Densities Summary
sampLL <- coda.samples(mLL,                      variable.names=varnames,                      n.iter=20000, progress.bar="none")
sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL)
X11()
plot(sampNL)
X11()
plot(sampDL,main="DL")