
###  Estimates the hidden function of Hunt Rate Vs Prey Desnsity for each group - Uses Non-parametric Gaussian Process with Bayesian Inference
### This File Generates the Gaussian Processs Plot with the correct uncesrtainty region ##
###
## 24/08/2018 Recalculated after Adding the D19 Set (Unlabelled Yet)- With the Experimental Repetitions of LL - 
## Note : I fix some obvious tracker Overestimates of Prey Density in
## the Empty Groups by caping to max 5 Prey (In reality these are debri sitting on the bottom of the dish)

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

#source("DataLabelling/labelHuntEvents_lib.r")
### GP Process Estimation Of Hunt Rate Vs Prey Density Using Bayesian Inference Model
myplot_res<- function(ind,qq=0.05){
  
  
  ### Show Speed Fit ###
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 3.2
  lineXAxis = 3.0
  cex = 1.4
  adj  = 3.5
  padj <- -8.0
  las <- 1
  
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(3.9,4.7,2,1))
  
  
  xplotLim <- c(0,60)
  yplotLim <- c(0,80)
  plot(foodlevelsLL,countsLL,col=colourLegL[2],cex=cex,cex.axis=cex,
       main = NA,
       ylab=NA,
       xlab=NA,
       xlim = xplotLim,
       ylim = yplotLim,
       pch=pchL[2],
       sub=paste("GP tau:",format(mean(drawLL$tau),digits=4 ),
                 "rho:",format(mean(drawLL$rho),digits=4 ) )  
       )
  
  
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Initial prey count in ROI (Tracker estimate)" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression(" Number of hunt events in 10 min" ))
  
  
  legend("topright",pch=pchL,cex=cex,
         legend = c(paste("NF #",nDatNL),paste("LF #",nDatLL),paste("DF #",nDatDL)),col=colourLegL)
  
  points(foodlevelsNL,countsNL,col=colourLegL[1],pch=pchL[1],xlim = xplotLim,cex=cex)
  points(foodlevelsDL,countsDL,col=colourLegL[3],pch=pchL[3],xlim = xplotLim,cex=cex)
  
  muLL=apply(apply(drawLL$lambda[,,1],1,tail,ind),2,mean)
  muNL=apply(apply(drawNL$lambda[,,1],1,tail,ind),2,mean)
  muDL=apply(apply(drawDL$lambda[,,1],1,tail,ind),2,mean)
  
  lines(foodlevelsLL,muLL,col=colourH[1],lwd=4,xlim = xplotLim)
  lines(foodlevelsNL,muNL,col=colourH[2],lwd=4,xlim = xplotLim)
  lines(foodlevelsDL,muDL,col=colourH[3],lwd=4,xlim = xplotLim)
  
  band=apply(apply(drawLL$lambda[,,1],1,tail,ind),2,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsLL,rev(foodlevelsLL)),c(band[1,],rev(band[2,])),col=colourR[1])
  
  band=apply(apply(drawNL$lambda[,,1],1,tail,ind),2,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsNL,rev(foodlevelsNL)),c(band[1,],rev(band[2,])),col=colourR[2])
  
  band=apply(apply(drawDL$lambda[,,1],1,tail,ind),2,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsDL,rev(foodlevelsDL)),c(band[1,],rev(band[2,])),col=colourR[3])
  
}


modelGPV1="model {
  # Likelihood

for(i in 1:N){
  n[i] ~ dnorm(lambda[i],eps)
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
    Sigma[i,j] <- pow(tau,2) * exp( - rho * pow(food[i] - food[j], 2) )
    #Sigma[i,j] <- pow(tau,2) * exp( - 0.5*rho^2 * pow(food[i] - food[j], 2) )
    #Sigma[i,j] <-  exp( - 0.5* pow((food[i] - food[j])*rho, 2) )
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

##Fixed Rho Model - From Gio 19-05-18 - We fixed rho since estimated both tau and rho made it very sensitive and noisy/non-smooth#
##He uses this to estimate Params for GP, and then constructs a classic GP with the estimated parameters
modelGPFixedRho="model {
  # Likelihood
  
  n ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:N) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
  
    for(j in (i+1):N) {
      Sigma[i,j] <- pow(tau,2) * exp( - 0.5*rho^2 * pow(food[i] - food[j], 2) )
      Sigma[j,i] <- Sigma[i,j]
    }
  }
 
  alpha=0 
  tau0 ~ dgamma(2,0.2) 
  tau  ~ dgamma(2,.2) 
  rho = 0.1
  
}"



#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

#load("out/setn-12-D-5-16-datHuntStat.RData")


#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStat <- makeHuntStat(datHuntLabelledEventsKL)

#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
#strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated-Merged" ##Has the Empty-Test Condition Fish Too (unlabelled)
#strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged" ##Warning Set Includes Repeated Test For some LF fish - One In Different Food Density

message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#datHuntLabelledEventsSBMerged <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
datHuntLabelledEventsSBMerged <- getLabelledHuntEventsSet()
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged)

##Add The Empty Test Conditions
#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
#datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]

## Get Event Counts Within Range ##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE) )
datHuntVsPreyLE[datHuntVsPreyLE[,1] > 5,1] <- 5
datHuntVsPreyL <- rbind(datHuntVsPreyLL,datHuntVsPreyLE)

datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE) )
datHuntVsPreyNE[datHuntVsPreyNE[,1] > 5,1] <- 5
datHuntVsPreyN <- rbind(datHuntVsPreyNL,datHuntVsPreyNE)

datHuntVsPreyN <- datHuntVsPreyN[!is.na(datHuntVsPreyN[,1]),]


datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE) )
datHuntVsPreyDE[datHuntVsPreyDE[,1] > 5,1] <- 5 ##Cap 
datHuntVsPreyD <- rbind(datHuntVsPreyDL,datHuntVsPreyDE)
##Remove NA 
datHuntVsPreyD <- datHuntVsPreyD[!is.na(datHuntVsPreyD[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,320)
colourH <- c(rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.01,0.01,0.9,0.8),rgb(0.00,0.00,0.0,1.0))
colourP <- c(rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.01,0.01,0.8,0.5),rgb(0.00,0.00,0.0,1.0))
colourR <- c(rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.01,0.01,0.9,0.4),rgb(0.00,0.00,0.0,1.0))
##Thse RC params Work Well to Smooth LF And NF
tauRangeA =50 #10000
rhoMaxA = 0.8
Noise = 1 ##The Gaussian Noise Term

burn_in=10;
steps=10000;
thin=2;


##Larva Event Counts Slice
datSliceLL <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2], ]
datSliceNL <- datHuntVsPreyN[datHuntVsPreyN[,1] >= preyCntRange[1] & datHuntVsPreyN[,1] <= preyCntRange[2], ]
datSliceDL <- datHuntVsPreyD[datHuntVsPreyD[,1] >= preyCntRange[1] & datHuntVsPreyD[,1] <= preyCntRange[2], ]
nDatLL <- length(datSliceLL[,1])
nDatNL <- length(datSliceNL[,1])
nDatDL <- length(datSliceDL[,1])


nEventsLL2=as.numeric(datSliceLL[seq(1,nDatLL),2])
nEventsNL2=as.numeric(datSliceNL[seq(1,nDatNL),2])
nEventsDL2=as.numeric(datSliceDL[seq(1,nDatDL),2])

nFoodLL2=as.numeric(datSliceLL[seq(1,nDatLL),1])
nFoodNL2=as.numeric(datSliceNL[seq(1,nDatNL),1])
nFoodDL2=as.numeric(datSliceDL[seq(1,nDatDL),1])


dataLL2=list(n=nDatLL,N=length(nDatLL),food=as.integer(nFoodLL2) );
dataNL2=list(n=nDatNL,N=length(nDatNL),food=as.integer(nFoodNL2) ) ;
dataDL2=list(n=nDatDL,N=length(nDatDL),food=as.integer(nFoodDL2) ) ;


foodlevelsLL=dataLL2[dataLL2$food >= preyCntRange[1] & dataLL2$food <= preyCntRange[2]]$food
ordLL=order(foodlevelsLL)
foodlevelsLL=foodlevelsLL[ordLL]
countsLL=nEventsLL2[ordLL]

foodlevelsNL=dataNL2[dataNL2$food >= preyCntRange[1] & dataNL2$food <= preyCntRange[2]]$food
ordNL=order(foodlevelsNL)
foodlevelsNL=foodlevelsNL[ordNL]
countsNL=nEventsNL2[ordNL]

foodlevelsDL=dataDL2[dataDL2$food >= preyCntRange[1] & dataDL2$food <= preyCntRange[2]]$food
ordDL=order(foodlevelsDL)
foodlevelsDL=foodlevelsDL[ordDL]
countsDL=nEventsDL2[ordDL]

dfoodlevelsNL <- density(foodlevelsNL,bw=0.5)
dfoodlevelsLL <- density(foodlevelsLL,bw=0.5)
dfoodlevelsDL <- density(foodlevelsDL,bw=0.5)


dataLL=list(n=countsLL,food=foodlevelsLL,N=length(countsLL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataNL=list(n=countsNL,food=foodlevelsNL,N=length(countsNL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataDL=list(n=countsDL,food=foodlevelsDL,N=length(countsDL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);

varnames=c("tau","rho","alpha","lambda")


library(rjags)
fileConn=file("model.tmp")
writeLines(modelGPV1,fileConn);
close(fileConn)

mLL=jags.model(file="model.tmp",data=dataLL);
mNL=jags.model(file="model.tmp",data=dataNL);
mDL=jags.model(file="model.tmp",data=dataDL);
update(mLL,burn_in);#update(mNL,burn_in);update(mDL,burn_in)


drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)

save(drawLL,drawNL,drawDL,file = paste0(strDataExportDir,"stat_GPProcessHuntRateVsPreyDensity_RJags.RData"))

strPlotName <-  paste(strPlotExportPath,"/stat/fig2S1-stat_HuntEventRateLabelledT30V50VsPrey_GPEstimate2-tauMax",tauRangeA,"Rho",rhoMaxA,".pdf",sep="")
pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 
  myplot_res(2500)

dev.off()

cex <- 1.4
plot(dfoodlevelsNL,lwd=4,lty=1,xlim=c(0,60),col=colourLegL[1],main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dfoodlevelsLL,lwd=4,lty=1,col=colourLegL[2])
lines(dfoodlevelsDL,lwd=4,lty=1,col=colourLegL[3])

## Compare Prey Density ###
par(mar = c(3.9,4.7,2,1))
strCumPlotName <-  paste(strPlotExportPath,"/stat/fig2S2-InitPreyCount_CDF.pdf",sep="")
pdf(strCumPlotName,width=8,height=8,title="Compare prey density testing conditions between groups") 
  plot(ecdf(foodlevelsNL),xlim=c(0,60),lwd=4,lty=1,col=colourLegL[1],main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex,pch=pchL[1])
  lines(ecdf(foodlevelsLL),xlim=c(0,60),lwd=4,lty=2,pch=pchL[2],col=colourLegL[2],cex=cex)
  lines(ecdf(foodlevelsDL),xlim=c(0,60),lwd=4,lty=3,pch=pchL[3],col=colourLegL[3],cex=cex)
  mtext(side = 1,cex=cex, line = 2.7, expression("Initial prey count in ROI (Tracker estimate)" ))
  mtext(side = 2,cex=cex, line = 2.2, expression(" Cumulative distribution " ))
  
  legend("bottomright",pch=pchL,cex=cex,
         legend = c(paste("NF #",nDatNL),paste("LF #",nDatLL),paste("DF #",nDatDL)),col=colourLegL)
dev.off()

###

X11()
myplot_res(4000)

