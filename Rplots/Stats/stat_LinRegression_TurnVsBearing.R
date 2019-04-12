##3-09-2018
###  Estimates the hidden function of Turn Vs Bearing To Prey - Linear Regression gives pdf of slope param.  - 
## Tried both a Non-parametric Gaussian Process with Bayesian Inference (Failed) But Also a simple Linear Model
## Requires the lFirstBoutPoints list of dataframes - which is constucted in 

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetC",".rds",sep="") #Processed Registry on which we add 

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
  beta[1] ~ dnorm(0,4)
  beta[2] ~ dnorm(10,4)
  

  # Prior for the inverse variance
  inv.var   ~ dgamma(10, 205)
  sigma     <- 1/sqrt(inv.var)

}"

####Select Subset Of Data To Analyse
plot(dgamma(1:100,            shape=10,scale=205           ))
plot(dnorm(1:100,   mean=1,sd=4           ))

datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
remove(lMotionBoutDat)
lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #Processed Registry on which we add )
remove(lFirstBoutPoints) ##Load From File
lFirstBoutPoints <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_SetC",".rds",sep="") ) #Processed Registry on which we add )

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
## Added Flag whether Capture Strike move Detected 
flagWithCaptureStrike <- 1
randomSubset <- 30
datTurnVsPreyLL <- cbind(OnSetAngleToPrey=lFirstBoutPoints$LL[,"OnSetAngleToPrey"] , Turn= as.numeric(lFirstBoutPoints$LL[,"Turn"]),RegistarIdx=lFirstBoutPoints$LL[,"RegistarIdx"],
                         CaptureStrikeDetected=datTrackedEventsRegister[lFirstBoutPoints$LL[,"RegistarIdx"],"CaptureStrikeDetected"] )
#datTurnVsPreyLL <- datTurnVsPreyLL[ sample(NROW(datTurnVsPreyLL),size=randomSubset), ]
#datTurnVsPreyLL <- datTurnVsPreyLL[!is.na(datTurnVsPreyLL[,1]) & datTurnVsPreyLL[,4] == flagWithCaptureStrike,]

datTurnVsPreyNL <- cbind(OnSetAngleToPrey=lFirstBoutPoints$NL[,"OnSetAngleToPrey"] ,Turn= as.numeric(lFirstBoutPoints$NL[,"Turn"]),lFirstBoutPoints$NL[,"RegistarIdx"],CaptureStrikeDetected=datTrackedEventsRegister[lFirstBoutPoints$NL[,"RegistarIdx"],"CaptureStrikeDetected"] )
#datTurnVsPreyNL <- datTurnVsPreyNL[ sample(NROW(datTurnVsPreyNL),size=randomSubset), ]
#datTurnVsPreyNL <- datTurnVsPreyNL[!is.na(datTurnVsPreyNL[,1]) & datTurnVsPreyNL[,4] == flagWithCaptureStrike,]


datTurnVsPreyDL <- cbind(OnSetAngleToPrey=lFirstBoutPoints$DL[,"OnSetAngleToPrey"] , Turn=  as.numeric(lFirstBoutPoints$DL[,"Turn"]),lFirstBoutPoints$DL[,"RegistarIdx"],CaptureStrikeDetected=datTrackedEventsRegister[lFirstBoutPoints$DL[,"RegistarIdx"],"CaptureStrikeDetected"] )
#datTurnVsPreyDL <- datTurnVsPreyDL[ sample(NROW(datTurnVsPreyDL),size=randomSubset), ]
#datTurnVsPreyDL <- datTurnVsPreyDL[!is.na(datTurnVsPreyDL[,1]) & datTurnVsPreyDL[,4] == flagWithCaptureStrike,]

##Outlier datTurnVsPreyDL[13,] <- NA

##For the 3 Groups 
#colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
#colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points]
#colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
#pchL <- c(16,2,4)
#
#Thse RC params Work Well to Smooth LF And NF
tauRangeA =100000 #10000
rhoMaxA = 1000
Noise = 1 ##The Gaussian Noise Term

burn_in=10;
steps=100000;
thin=1;
nchains <- 3



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

mLL=jags.model(file="model.tmp",data=dataLL,n.chains = nchains);
mNL=jags.model(file="model.tmp",data=dataNL,n.chains = nchains);
mDL=jags.model(file="model.tmp",data=dataDL,n.chains = nchains);
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
muNLb=mean(drawNL$beta[,(steps-ind):steps,1][2,]) #Slope
muDLa=mean(drawDL$beta[,(steps-ind):steps,1][1,])
muDLb=mean(drawDL$beta[,(steps-ind):steps,1][2,])
sig=mean(drawLL$sigma[,(steps-ind):steps,1])
###Plot Density of Slope
pBw <- 0.01
dLLb<-density(drawLL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)
dNLb<-density(drawNL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)
dDLb<-density(drawDL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)

##Open Output PDF 
#pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootLinRegressions_Cap",G_THRES_CAPTURE_SPEED,"Strike",flagWithCaptureStrike,"_SetC2.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")
pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootLinRegressions_SetC2.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")

outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 3.5
padj <- -23.0
las <- 1

layout(matrix(c(1,2),1,2, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,1,1))
plot( datTurnVsPreyDL[,"OnSetAngleToPrey"],datTurnVsPreyDL[,"Turn"],
     main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
     ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
     xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[3] ,pch=pchL[3]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
mtext(side = 1,cex=0.8, line = 2.2, expression("Bearing To Prey Prior Turn "~(phi^degree) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Bearing To Prey After Turn "~(theta^degree) ))

#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(datTurnVsPreyDL[,"Turn"] ~ datTurnVsPreyDL[,"OnSetAngleToPrey"]),col=colourH[3],lwd=2.0,lty=2) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
#abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
#abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(datTurnVsPreyLL[,"OnSetAngleToPrey"],datTurnVsPreyLL[,"Turn"],pch=pchL[2],col=colourP[2],cex=1.2)
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(datTurnVsPreyLL[,"Turn"] ~ datTurnVsPreyLL[,"OnSetAngleToPrey"]),col=colourH[2],lwd=2.0,lty=2)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5,lty=1.5) ##Fit Line / Regression
#abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
#abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points( datTurnVsPreyNL[,"OnSetAngleToPrey"],datTurnVsPreyNL[,"Turn"],pch=pchL[1],col=colourP[1])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(datTurnVsPreyNL[,"Turn"] ~ datTurnVsPreyNL[,"OnSetAngleToPrey"]),col=colourH[1],lwd=2.0,lty=2)
abline(a=muNLa,b=muNLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
#abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)

##Guiding lines / Draw 0 Vertical Line
segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);

legend("topleft",
       legend=c(  expression (),
                         bquote(NF["e"] ~ '#' ~ .(NROW(datTurnVsPreyNL[,"Turn"]))  ),
                         bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  ),
                         bquote(DF["e"] ~ '#' ~ .(NROW(datTurnVsPreyDL[,"Turn"]))  )  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL, col=colourLegL)
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)


##Density Estimation
plot(dNLb,col=colourLegL[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA) #expression(paste("slope ",gamma) ) )
lines(dLLb,col=colourLegL[2],xlim=c(0.5,1.2),lwd=3,lty=2)
lines(dDLb,col=colourLegL[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topright",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(NROW(datTurnVsPreyNL[,"Turn"]))  ),
                  bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  ),
                  bquote(DF["e"] ~ '#' ~ .(NROW(datTurnVsPreyDL[,"Turn"]))  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3),lwd=3)
mtext(side = 1,cex=0.8, line = 2.2, expression(paste("slope ",gamma) ))
mtext(side = 2,cex=0.8, line = 2.2, expression("Density ") )
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

### PLot Scatter with regression lines with Conf intervals##
dev.off()

pdf(file= paste(strPlotExportPath,"/stat/boxplot_UndershootRatio_Cap",G_THRES_CAPTURE_SPEED,"Strike",flagWithCaptureStrike,"_SetC.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")
#pdf(file= paste(strPlotExportPath,"/stat/boxplot_UndershootRatio_RandSub_SetC.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")

boxplot(datTurnVsPreyNL[,"Turn"]/datTurnVsPreyNL[,"OnSetAngleToPrey"],
        datTurnVsPreyLL[,"Turn"]/datTurnVsPreyLL[,"OnSetAngleToPrey"],
        datTurnVsPreyDL[,"Turn"]/datTurnVsPreyDL[,"OnSetAngleToPrey"],
        names=c("NL","LL","DL"),main="Undershoot slope ",col=colourH)
dev.off()
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
