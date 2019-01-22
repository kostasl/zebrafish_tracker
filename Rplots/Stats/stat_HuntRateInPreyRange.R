### Makes Bayssian inference on  hunt events rate - within a given rannge of prey counts##

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


##### Do Jags Statistics ###
modelI="model { 
qq ~ dnorm(10,0.001)T(0,400)
tt ~ dnorm(0,0.001)T(0,100)

for(j in 1:NTOT){
q[j] ~ dnorm(qq,tt) 
n[j] ~ dpois(q[j])
}
}"

modelG="model { 
q ~ dnorm(15,0.0001)T(0,400)

for(j in 1:NTOT){
n[j] ~ dpois(q)
}


}"

fromchain=1000

#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

colourH <- c(rgb(0.01,0.7,0.01,0.2),rgb(0.9,0.01,0.01,0.2),rgb(0.01,0.01,0.9,0.2),rgb(0.00,0.00,0.0,1.0))

strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged" ##Warning Set Includes Repeated Test For some LF fish - One In Different Food Density

message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSBMerged <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged)


## Get Event Counts Within Range ##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE) )
datHuntVsPreyL <- rbind(datHuntVsPreyLL,datHuntVsPreyLE)

datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE) )
datHuntVsPreyN <- rbind(datHuntVsPreyNL,datHuntVsPreyNE)

datHuntVsPreyN <- datHuntVsPreyN[!is.na(datHuntVsPreyN[,1]),]


datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE) )
datHuntVsPreyD <- rbind(datHuntVsPreyDL,datHuntVsPreyDE)
##Remove NA 
datHuntVsPreyD <- datHuntVsPreyD[!is.na(datHuntVsPreyD[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,10)

##Larva Event Counts Slice
datSliceLL <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2],2 ]
datSliceNL <- datHuntVsPreyN[datHuntVsPreyN[,1] >= preyCntRange[1] & datHuntVsPreyN[,1] <= preyCntRange[2],2 ]
datSliceDL <- datHuntVsPreyD[datHuntVsPreyD[,1] >= preyCntRange[1] & datHuntVsPreyD[,1] <= preyCntRange[2],2 ]
nDatLL <- length(datSliceLL)
nDatNL <- length(datSliceNL)
nDatDL <- length(datSliceDL)

nLL2=as.numeric(datSliceLL[seq(1,nDatLL,1)])
nNL2=as.numeric(datSliceNL[seq(1,nDatNL,1)])
nDL2=as.numeric(datSliceDL[seq(1,nDatDL,1)])

dataLL2=list(n=nLL2,NTOT=length(nLL2),food=as.integer(datHuntVsPreyL[,1]));
dataNL2=list(n=nNL2,NTOT=length(nNL2),food=as.integer(datHuntVsPreyN[,1]));
dataDL2=list(n=nDL2,NTOT=length(nDL2),food=as.integer(datHuntVsPreyD[,1]));

varnames1=c("n","q")
burn_in=1000;
steps=100000;
thin=2;

library(rjags)
strModelName = "modelG.tmp"
fileConn=file(strModelName)
writeLines(modelG,fileConn);
close(fileConn)

mLL2=jags.model(file=strModelName,data=dataLL2);
mNL2=jags.model(file=strModelName,data=dataNL2);
mDL2=jags.model(file=strModelName,data=dataDL2);

update(mLL2,burn_in)
update(mNL2,burn_in)
update(mDL2,burn_in)

drawLL2=jags.samples(mLL2,steps,thin=thin,variable.names=varnames1)
drawNL2=jags.samples(mNL2,steps,thin=thin,variable.names=varnames1)
drawDL2=jags.samples(mDL2,steps,thin=thin,variable.names=varnames1)

##q[idx,sampleID,chainID]
##PLot The Param Mean Distributions
#strPlotName <- paste("plots/stat_HuntEventRateParamPreyRange",preyCntRange[1],preyCntRange[2], ".pdf",sep="-")
#pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 
X11()

hist(drawLL2$q[1,,1],breaks=seq(0,15,length=100),col=colourH[1],xlab="Hunt Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
hist(drawNL2$q[1,,1],breaks=seq(0,15,length=100),add=T,col=colourH[2])
hist(drawDL2$q[1,,1],breaks=seq(0,15,length=100),add=T,col=colourH[3])

legend("topright",legend = c(paste("LL #",nDatLL),paste("NL #",nDatNL),paste("DL #",nDatDL)),fill=colourH)

#dev.off()
##hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
##hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
##hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

### Do Linear Rate Vs Prey Model Plot -- Add Quantile Plots of Uncertainty in Parameters##
X11()
strDensityplotFileName <- paste("plots/stat_HuntRateVsInitPreyCountModel.pdf",collapse=NULL,sep="");
#pdf(strDensityplotFileName,width=8,height=8)

ind=length(drawLL2$mu[1,,])
food_level = seq(0,80,0.2)

mu=apply(drawLL2$mu[,(ind-1000):ind,],1,mean)
qmu1 = drawLL2$mu[1,(ind-1000):ind,]
qmu2 = drawLL2$mu[2,(ind-1000):ind,]
lowup=matrix(NA,2,length(food_level))
for(i in 1:length(food_level)){
  lowup[,i]=quantile(qmu1+qmu2*food_level[i],probs=c(0.05,0.95))
}
plot(food_level,mu[1]+mu[2]*food_level,main="Hunt Rate Estimation Vs Initial Prey Count",ylim=c(0,40),t='l',col=colourH[1],cex=6.5,lwd=4)
polygon(c(food_level,rev(food_level)),c(lowup[1,],rev(lowup[2,])),ylim=c(0,40),col=colourH[1],cex=6.5)

points(dataLL2$food,dataLL2$n,col=colourH[4] ,t='p')
