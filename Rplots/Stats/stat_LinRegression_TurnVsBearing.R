## kostasl 3-09-2018
###  Estimates the hidden function of Turn Vs Bearing To Prey - Linear Regression gives pdf of slope param.  - 
## Tried both a Non-parametric Gaussian Process with Bayesian Inference (Failed) But Also a simple Linear Model
## Requires the lFirstBoutPoints list of dataframes - which is constucted in 
## 18-03-19 Update : Add confidence 5% intervals 
source("config_lib.R")
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
  beta[1] ~ dnorm(0,2)
  beta[2] ~  dnorm(1,1/sqrt(sigmaU))T(0.0,2) ##undershoot
  sigmaU ~ dgamma(1, 1) ##dunif(0.0,0.50)

  # Prior for the inverse variance
  inv.var   ~  dgamma(5, 2)
  sigma     <- 1/sqrt(inv.var)

}"

####Select Subset Of Data To Analyse
plot(dgamma(1:100,            shape=5,scale=2           ))

plot(dnorm(1:100,   mean=1,sd=1           ))

datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
remove(lMotionBoutDat)
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #Processed Registry on which we add )
remove(lFirstBoutPoints) ##Load From File
#lFirstBoutPoints <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_SetC",".rds",sep="") ) #Processed Registry on which we add )

## The Original list if the lFirstBout data from runHuntepisode analysis
datMotionBoutsToValidate <-readRDS(file=paste0(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_ToValidate.rds") ) 
lFirstBoutPoints <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated.rds",sep="") ) #Processed Registry on which we add )


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
                         CaptureStrikeDetected=lFirstBoutPoints$LL[,"doesCaptureStrike"] ) ##datTrackedEventsRegister[lFirstBoutPoints$NL[,"RegistarIdx"],"CaptureStrikeDetected"]
#datTurnVsPreyLL <- datTurnVsPreyLL[ sample(NROW(datTurnVsPreyLL),size=randomSubset), ]
#datTurnVsPreyLL <- datTurnVsPreyLL[!is.na(datTurnVsPreyLL[,1]) & datTurnVsPreyLL[,4] == flagWithCaptureStrike,]

datTurnVsPreyNL <- cbind(OnSetAngleToPrey=lFirstBoutPoints$NL[,"OnSetAngleToPrey"] ,Turn= as.numeric(lFirstBoutPoints$NL[,"Turn"]),
                         lFirstBoutPoints$NL[,"RegistarIdx"],
                         CaptureStrikeDetected=lFirstBoutPoints$NL[,"doesCaptureStrike"] )
#datTurnVsPreyNL <- datTurnVsPreyNL[ sample(NROW(datTurnVsPreyNL),size=randomSubset), ]
#datTurnVsPreyNL <- datTurnVsPreyNL[!is.na(datTurnVsPreyNL[,1]) & datTurnVsPreyNL[,4] == flagWithCaptureStrike,]


datTurnVsPreyDL <- cbind(OnSetAngleToPrey=lFirstBoutPoints$DL[,"OnSetAngleToPrey"] , Turn=  as.numeric(lFirstBoutPoints$DL[,"Turn"]),lFirstBoutPoints$DL[,"RegistarIdx"],
                         CaptureStrikeDetected=lFirstBoutPoints$DL[,"doesCaptureStrike"] )
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
steps=50000;
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
varnames=c("beta","sigma","sigmaU")

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
ind = steps*0.10 ## Number of last sampled values
## Save the Mean Slope and intercept
##quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2]
muLLa=mean(drawLL$beta[,(steps-ind):steps,1][1,]) 
muLLb=mean(drawLL$beta[,(steps-ind):steps,1][2,])
muNLa=mean(drawNL$beta[,(steps-ind):steps,1][1,])
muNLb=mean(drawNL$beta[,(steps-ind):steps,1][2,]) #Slope
muDLa=mean(drawDL$beta[,(steps-ind):steps,1][1,])
muDLb=mean(drawDL$beta[,(steps-ind):steps,1][2,])
#sig=mean(drawLL$sigma[,(steps-ind):steps,1])
sigLL=mean(drawLL$sigmaU[,(steps-ind):steps,1])

###Plot Density of Slope
pBw <- 0.01
dLLb<-density(drawLL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)
dNLb<-density(drawNL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)
dDLb<-density(drawDL$beta[,(steps-ind):steps,1][2,],kernel="gaussian",bw=pBw)

## Compare Turn Ration Between Groups 
##Add Density That LF undershoot More than NF/DF
turnRatio_LFvsNF <- (drawLL$beta[,(steps-ind):steps,1][2,])-drawNL$beta[,(steps-ind):steps,1][2,]
turnRatio_LFvsDF <- (drawLL$beta[,(steps-ind):steps,1][2,])-drawDL$beta[,(steps-ind):steps,1][2,]
turnRatio_NFvsDF <- (drawNL$beta[,(steps-ind):steps,1][2,])-drawDL$beta[,(steps-ind):steps,1][2,]

dLLbVsNF <- density(turnRatio_LFvsNF,kernel="gaussian",bw=pBw)
dLLbVsDF <- density(turnRatio_LFvsDF,kernel="gaussian",bw=pBw)
dNLbVsDF <- density(turnRatio_NFvsDF,kernel="gaussian",bw=pBw)
dNLbVsNF <- density(drawNL$beta[,(steps-ind):steps,1][2,]-sample(drawNL$beta[,(steps-ind):steps,1][2,]),kernel="gaussian",bw=pBw)

PUndershoot_LFgtNF <- length(turnRatio_LFvsNF[turnRatio_LFvsNF < 0])/length(turnRatio_LFvsNF) ## ProbValLessThan(dLLbVsNF,0)
PUndershoot_LFgtDF <- length(turnRatio_LFvsDF[turnRatio_LFvsDF < 0])/length(turnRatio_LFvsDF) ## ProbValLessThan(dLLbVsNF,0)
PUndershoot_NFgtDF <- length(turnRatio_NFvsDF[turnRatio_NFvsDF < 0])/length(turnRatio_NFvsDF) ## ProbValLessThan(dLLbVsNF,0)

print( paste("Prob that LF pooled data have stronger undershoot than NF:", PUndershoot_LFgtNF ) )
print( paste("Prob that LF pooled data have lower undershoot than DF:", PUndershoot_LFgtDF ) )
print( paste("Prob that NF pooled data have lower undershoot than DF:", PUndershoot_NFgtDF ) )
print( paste("(Control-Validation)Prob that NF pooled data have lower undershoot than DN:", ProbValLessThan(dNLbVsNF,0) ))

###Density of STD Dev on TurnRatio
dsigLL=density(drawLL$sigmaU[,(steps-ind):steps,1])  
dsigDL=density(drawDL$sigmaU[,(steps-ind):steps,1])  
dsigNL=density(drawNL$sigmaU[,(steps-ind):steps,1])  

gammaLL <- dataLL$turn/ dataLL$bearing
gammaNL <- dataNL$turn/ dataNL$bearing
gammaDL <- dataNL$turn/ dataNL$bearing
###Plot DATA Density of Slope
pBw <- 0.2
dDatLLb<-density(gammaLL,kernel="gaussian",bw=pBw)
dDatNLb<-density(gammaNL,kernel="gaussian",bw=pBw)
dDatDLb<-density(gammaDL,kernel="gaussian",bw=pBw)



  ##### ######################
  ### MAIN FIGURE ############
  
  ################################  
pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_SetC2.pdf",sep="")
    ,width=7,height=7,
    title="First Turn To prey / Undershoot Ratio")
  
  outer = FALSE
  line = 1 ## SubFig Label Params
  line <- 2.6 ## SubFig Label Params
  lineAxis = line## 3.2    
  cex = 1.4
  adj  = 3.5
  padj <- -16.5
  las <- 1
  
  #layout(matrix(c(1,2),1,2, byrow = FALSE))
  ##Margin: (Bottom,Left,Top,Right )
  #par(mar = c(3.95,4.75,1,1))
  par(mar = c(4.2,4.8,1.1,1))
  plot( datTurnVsPreyDL[,"OnSetAngleToPrey"],datTurnVsPreyDL[,"Turn"],
       main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
       xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
       ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
       xlim=c(-100,100),
       ylim=c(-100,100),
       col=colourP[3] ,pch=pchL[3],cex=cex,cex.axis=cex) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
  mtext(side = 1,cex=cex, line = lineAxis, expression("Prey azimuth prior to turn ("~theta^degree~")" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Turn to prey ("~phi^degree~")" ))
  
  #text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
  ##abline(lm(datTurnVsPreyDL[,"Turn"] ~ datTurnVsPreyDL[,"OnSetAngleToPrey"]),col=colourLegL[3],lwd=1.5,lty=2) ##Fit Line / Regression
  abline(a=muDLa,b=muDLb,col=colourLegL[3],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[3],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[3],lwd=4.0) ##Fit Line / Regression
  
  #abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
  ##LL
  points(datTurnVsPreyLL[,"OnSetAngleToPrey"],datTurnVsPreyLL[,"Turn"],pch=pchL[2],col=colourP[2],cex=1.2)
  #text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
  #abline(lm(datTurnVsPreyLL[,"Turn"] ~ datTurnVsPreyLL[,"OnSetAngleToPrey"]),col=colourLegL[2],lwd=1.5,lty=2)
  abline(a=muLLa,b=muLLb,col=colourLegL[2],lwd=4.0,lty=1.5) ##Fit Line / Regression
  abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[2],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[2],lwd=4.0) ##Fit Line / Regression
  
  #abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
  ##NL
  points( datTurnVsPreyNL[,"OnSetAngleToPrey"],datTurnVsPreyNL[,"Turn"],pch=pchL[1],col=colourP[1])
  #text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
  #abline(lm(datTurnVsPreyNL[,"Turn"] ~ datTurnVsPreyNL[,"OnSetAngleToPrey"]),col=colourLegL[1],lwd=1.5,lty=2)
  abline(a=muNLa,b=muNLb,col=colourLegL[1],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[1],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[1],lwd=4.0) ##Fit Line / Regression
  
  ##Guiding lines / Draw 0 Vertical Line
  segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);
  
  legend("topleft",
         legend=c(  expression (),
                           bquote(NF["e"] ~ '#' ~ .(NROW(datTurnVsPreyNL[,"Turn"]))  ),
                           bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  ),
                           bquote(DF["e"] ~ '#' ~ .(NROW(datTurnVsPreyDL[,"Turn"]))  )  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         pch=pchL, col=colourLegL,cex=cex)
  #mtext("C",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
  
dev.off()  

pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_Density.pdf",sep=""),
    width=7,height=7,
    title="First Turn To prey / Undershoot Ratio - Posterior Density")
  par(mar = c(4.2,4.8,1.1,1))
  ##Density Estimation
  plot(dNLb,col=colourLegL[1],xlim=c(0.0,2),lwd=4,lty=1,ylim=c(0,18),
       main=NA, #"Density Inference of Turn-To-Prey Slope ",
       xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
  lines(dLLb,col=colourLegL[2],xlim=c(0.5,1.2),lwd=4,lty=2)
  lines(dDLb,col=colourLegL[3],xlim=c(0.5,1.2),lwd=4,lty=3)
  
  #lines(dDatNLb,col=colourLegL[1],xlim=c(0.5,1.2),lwd=2,lty=1,ylim=c(0,20),
  #     main=NA, #"Density Inference of Turn-To-Prey Slope ",
  #     xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
  #lines(dDatLLb,col=colourLegL[2],xlim=c(0.5,1.2),lwd=2,lty=2)
  #lines(dDatDLb,col=colourLegL[3],xlim=c(0.5,1.2),lwd=2,lty=3)
  
  legend("topright",
         legend=c(  expression (),
                    bquote(NF["e"] ~ '#' ~ .(NROW(datTurnVsPreyNL[,"Turn"]))  ),
                    bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  ),
                    bquote(DF["e"] ~ '#' ~ .(NROW(datTurnVsPreyDL[,"Turn"]))  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL,lty=c(1,2,3),lwd=3,cex=cex)
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Estimated mean turn-ratio (",gamma,")" ) ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
 # mtext("D",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
  
  print(paste("LF has different mean with Prob "))
  ProbValLessThan(dLLb,0.7)
  ### PLot Scatter with regression lines with Conf intervals##
dev.off()

##############################################
### Individually Plotted Turn Ratios With   ##
## With COFINDENCE Intervals 5%             ##
##############################################
pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_DF.pdf",sep="")
    ,width=7,height=7,
    title="First Turn To prey / Undershoot Ratio Dry Fed Group")
  
  outer = FALSE
  line = 1 ## SubFig Label Params
  line <- 2.6 ## SubFig Label Params
  lineAxis = line## 3.2    
  cex = 1.4
  adj  = 3.5
  padj <- -16.5
  las <- 1
  
  par(mar = c(4.2,4.8,1.1,1))
  plot( datTurnVsPreyDL[,"OnSetAngleToPrey"],datTurnVsPreyDL[,"Turn"],
        main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
        xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
        ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
        xlim=c(-100,100),
        ylim=c(-100,100),
        col=colourP[3] ,pch=pchL[3],cex=cex,cex.axis=cex) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
  
  #text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
  ##abline(lm(datTurnVsPreyDL[,"Turn"] ~ datTurnVsPreyDL[,"OnSetAngleToPrey"]),col=colourLegL[3],lwd=1.5,lty=2) ##Fit Line / Regression
  abline(a=muDLa,b=muDLb,col=colourLegL[3],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[3],lwd=4.0) ##Fit Line / Regression
  abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[3],lwd=4.0) ##Fit Line / Regression
  ##Guiding lines / Draw 0 Vertical Line
  segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(DF["e"] ~ '#' ~ .(NROW(datTurnVsPreyDL[,"Turn"]))  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL[3],lty=NA,pch=pchL[3],lwd=3,cex=cex)
  
  mtext(side = 1,cex=cex, line = lineAxis, expression("Prey azimuth prior to turn ("~theta^degree~")" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Turn to prey ("~phi^degree~")" ))
  
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_NF.pdf",sep="")
    ,width=7,height=7,
    title="First Turn To prey / Undershoot Ratio NOT Fed Group")

outer = FALSE
line = 1 ## SubFig Label Params
line <- 2.6 ## SubFig Label Params
lineAxis = line## 3.2    
cex = 1.4
adj  = 3.5
padj <- -16.5
las <- 1

par(mar = c(4.2,4.8,1.1,1))
plot( datTurnVsPreyNL[,"OnSetAngleToPrey"],datTurnVsPreyNL[,"Turn"],pch=pchL[1],col=colourP[1],
      main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
      xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
      ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
      xlim=c(-100,100),
      ylim=c(-100,100),
      cex=cex,cex.axis=cex  )

#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
#abline(lm(datTurnVsPreyNL[,"Turn"] ~ datTurnVsPreyNL[,"OnSetAngleToPrey"]),col=colourLegL[1],lwd=1.5,lty=2)
abline(a=muNLa,b=muNLb,col=colourLegL[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[1],lwd=4.0) ##Fit Line / Regression

##Guiding lines / Draw 0 Vertical Line
segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);

legend("topleft",
       legend=c(  expression (),
                  bquote(NF["e"] ~ '#' ~ .(NROW(datTurnVsPreyNL[,"Turn"]) ) ) ) , ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL[1],lty=NA,pch=pchL[1],lwd=3,cex=cex)

mtext(side = 1,cex=cex, line = lineAxis, expression("Prey azimuth prior to turn ("~theta^degree~")" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Turn to prey ("~phi^degree~")" ))

dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_LF.pdf",sep="")
    ,width=7,height=7,
    title="First Turn To prey / Undershoot Ratio Live Fed Group")
  
  outer = FALSE
  line = 1 ## SubFig Label Params
  line <- 2.6 ## SubFig Label Params
  lineAxis = line## 3.2    
  cex = 1.4
  adj  = 3.5
  padj <- -16.5
  las <- 1
  
  par(mar = c(4.2,4.8,1.1,1))
  
  plot(datTurnVsPreyLL[,"OnSetAngleToPrey"],datTurnVsPreyLL[,"Turn"],pch=pchL[2],col=colourP[2],
         main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
         xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
         ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
         xlim=c(-100,100),
         ylim=c(-100,100),
         cex=cex,cex.axis=cex )
    #text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
    #abline(lm(datTurnVsPreyLL[,"Turn"] ~ datTurnVsPreyLL[,"OnSetAngleToPrey"]),col=colourLegL[2],lwd=1.5,lty=2)
    abline(a=muLLa,b=muLLb,col=colourLegL[2],lwd=4.0,lty=1.5) ##Fit Line / Regression
    abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[1],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[1],col=colourLegE[2],lwd=4.0) ##Fit Line / Regression
    abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,],p=c(0.05, 0.95))[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,],p=c(0.05, 0.95))[2],col=colourLegE[2],lwd=4.0) ##Fit Line / Regression
  
  ##Guiding lines / Draw 0 Vertical Line
  segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL[2],lty=NA,pch=pchL[2],lwd=3,cex=cex)
  
  mtext(side = 1,cex=cex, line = lineAxis, expression("Prey azimuth prior to turn ("~theta^degree~")" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Turn to prey ("~phi^degree~")" ))

dev.off()


# DATA Turn Ratio
################################  
pdf(file= paste(strPlotExportPath,"/stat/fig4S3_stat_TurnRatioHistogram.pdf",sep=""),width=7,height=7,title="First Turn To prey / Turn Ratio")
  
 #### PLOT Histogram from DATA
  layout(matrix(c(1,2,3),3,1, byrow = FALSE))
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(3.95,4.75,1,1))
  
  hist(dataLL$turn/ dataLL$bearing,breaks=seq(0,2,0.2),col=colourR[2],main="LF",xlab=NA,cex=cex,cex.axis=cex)
  hist(dataNL$turn/ dataNL$bearing,breaks=seq(0,2,0.2),col=colourR[1],main="NF",xlab=NA,cex=cex,cex.axis=cex)
  hist(dataDL$turn/ dataDL$bearing,breaks=seq(0,2,0.2),col=colourR[3],main="DF",xlab=NA,cex=cex,cex.axis=cex)
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio") )
dev.off()
  
###################### Further Turn Analysis #########################
##### Check out Undershoot Consistency - Ratio of Angle to Prey At start and end of bout ####
  layout(matrix(c(1,2,3),3,1, byrow = FALSE))
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.5,4.3,0.5,1))

  idxReg <- row.names(datTrackedEventsRegister[datTrackedEventsRegister$groupID == "NL",])
  ### Calc Turn Ratio
  boutTurnRatio_NL <- (datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey -
                        datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OffSetAngleToPrey )/ datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey
  
  boutTurnRatio_NL <- boutTurnRatio_NL[which(boutTurnRatio_NL > -2 & boutTurnRatio_NL < 2)]
  #plot(boutTurnRatioF,ylab="On/Off",)
  hist(boutTurnRatio_NL,xlim=c(-3,3),breaks=50)
  
  idxReg <- row.names(datTrackedEventsRegister[datTrackedEventsRegister$groupID == "LL",])
  boutTurnRatio_LL <-  (datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey -
                          datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OffSetAngleToPrey )/ datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey
  
  boutTurnRatio_LL <- boutTurnRatio_LL[which(boutTurnRatio_LL > -2 & boutTurnRatio_LL < 2)]
  #plot(boutTurnRatioF,ylab="On/Off",)
  hist(boutTurnRatio_LL,xlim=c(-3,3),breaks=50)
  
    
  idxReg <- row.names(datTrackedEventsRegister[datTrackedEventsRegister$groupID == "DL",])
  boutTurnRatio_DL <-  (datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey -
                          datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OffSetAngleToPrey )/ datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx %in% idxReg,]$OnSetAngleToPrey
  
  boutTurnRatio_DL <- boutTurnRatio_DL[which(boutTurnRatio_DL > -2 & boutTurnRatio_DL < 2)]
  #plot(boutTurnRatioF,ylab="On/Off",)
  hist(boutTurnRatio_DL,xlim=c(-3,3),breaks=50)
  
  
  
  boxplot(boutTurnRatio_NL,boutTurnRatio_LL,boutTurnRatio_DL)
  
  plot(datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx==idxReg,]$boutRank,
       datMotionBoutsToValidate[datMotionBoutsToValidate$RegistarIdx==idxReg,]$OffSetAngleToPrey)
  



pdf(file= paste(strPlotExportPath,"/stat/boxplot_UndershootRatio_Cap",G_THRES_CAPTURE_SPEED,"Strike",flagWithCaptureStrike,"_SetC.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")
#pdf(file= paste(strPlotExportPath,"/stat/boxplot_UndershootRatio_RandSub_SetC.pdf",sep=""),width=14,height=7,title="First Turn To prey / Undershoot Ratio")

boxplot(datTurnVsPreyNL[,"Turn"]/datTurnVsPreyNL[,"OnSetAngleToPrey"],
        datTurnVsPreyLL[,"Turn"]/datTurnVsPreyLL[,"OnSetAngleToPrey"],
        datTurnVsPreyDL[,"Turn"]/datTurnVsPreyDL[,"OnSetAngleToPrey"],
        names=c("NL","LL","DL"),main="Undershoot slope ",col=colourH,ylim=c(0,2.0) )
dev.off()




### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig4S1-DetectionAngleDensity_LF.pdf",sep=""))
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


#plot(density(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[1],xlim=c(-120.0,120),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA)
plot(density(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[2],xlim=c(-120.0,120),lwd=3,lty=2,main=NA,xlab=NA,ylab=NA)
#lines(density(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[3],xlim=c(-120.0,120),lwd=3,lty=3,main=NA)

legend("topleft",
       legend=c(  expression (),
                  #bquote(NF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$NL))  ),
                  bquote(LF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$LL))  ) ),
                  #bquote(DF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$DL))  )
       
       col=colourLegL[2],lty=c(2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection  " ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=c
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











####
## LF Only for MM grant proposal 

pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_SetC2_LF.pdf",sep="")
    ,width=7,height=7,
    title="First Turn To prey / Undershoot Ratio")

outer = FALSE
line = 1 ## SubFig Label Params
line <- 2.6 ## SubFig Label Params
lineAxis = line## 3.2    
cex = 1.4
adj  = 3.5
padj <- -16.5
las <- 1

#layout(matrix(c(1,2),1,2, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.95,4.75,1,1))
par(mar = c(4.2,4.8,1.1,1))
plot( datTurnVsPreyLL[,"OnSetAngleToPrey"],datTurnVsPreyLL[,"Turn"],
      main=NA,#paste("Turn Size Vs Bearing To Prey ", sep=""),
      xlab=NA,#expression("Bearing To Prey Prior Turn "~(phi^degree) ),
      ylab=NA,#expression("Bearing To Prey After Turn "~(theta^degree) ),
      xlim=c(-100,100),
      ylim=c(-100,100),
      col=colourP[2] ,pch=pchL[2],cex=cex,cex.axis=cex) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
mtext(side = 1,cex=cex, line = lineAxis, expression("Prey azimuth prior to turn ("~theta^degree~")" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Turn to prey ("~phi^degree~")" ))

abline(lm(datTurnVsPreyLL[,"Turn"] ~ datTurnVsPreyLL[,"OnSetAngleToPrey"]),col=colourLegL[2],lwd=1.5,lty=2)
abline(a=muLLa,b=muLLb,col=colourLegL[2],lwd=4.0,lty=1.5) ##Fit Line / Regression

##Guiding lines / Draw 0 Vertical Line
segments(0,-100,0,100); segments(-100,0,100,0); segments(-100,-100,100,100,lwd=1,lty=2);

legend("topleft",
       legend=c(  expression (),
                  
                  bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  )),
                     #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       pch=pchL[2], col=colourLegL[2],cex=cex)
#mtext("C",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)

dev.off()  

pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_UndershootLinRegressions_Density_LF.pdf",sep=""),
    width=7,height=7,
    title="First Turn To prey / Undershoot Ratio - Posterior Density")
par(mar = c(4.2,4.8,1.1,1))
##Density Estimation
#plot(dNLb,col=colourLegL[1],xlim=c(0.0,2),lwd=4,lty=1,ylim=c(0,18),
     main=NA, #"Density Inference of Turn-To-Prey Slope ",
     xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )

plot(dLLb,col=colourLegL[2],lwd=4,lty=2,xlim=c(0.0,2),ylim=c(0,18),
      main=NA, #"Density Inference of Turn-To-Prey Slope ",
      xlab=NA,ylab=NA,cex=cex,cex.axis=cex)

#lines(dDatNLb,col=colourLegL[1],xlim=c(0.5,1.2),lwd=2,lty=1,ylim=c(0,20),
#     main=NA, #"Density Inference of Turn-To-Prey Slope ",
#     xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
#lines(dDatLLb,col=colourLegL[2],xlim=c(0.5,1.2),lwd=2,lty=2)
#lines(dDatDLb,col=colourLegL[3],xlim=c(0.5,1.2),lwd=2,lty=3)

legend("topright",
       legend=c(  expression (),
                  bquote(LF["e"] ~ '#' ~ .(NROW(datTurnVsPreyLL[,"Turn"]))  )
                  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL[2],lty=c(2),lwd=3,cex=cex)
mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Estimated turn ratio ",(gamma) ) ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
# mtext("D",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)

### PLot Scatter with regression lines with Conf intervals##
dev.off()





