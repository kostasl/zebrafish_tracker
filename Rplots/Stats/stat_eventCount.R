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

nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

## Get Event Counts Within Range ##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )

### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)
nDat <- nrow(datHuntVsPreyLL)
datSliceLL <- datHuntVsPreyLL[datHuntVsPreyLL[,1] > preyCntRange[1] & datHuntVsPreyLL[,1] < preyCntRange[2],2 ]
datSliceNL <- datHuntVsPreyNL[datHuntVsPreyNL[,1] > preyCntRange[1] & datHuntVsPreyNL[,1] < preyCntRange[2],2 ]
datSliceDL <- datHuntVsPreyDL[datHuntVsPreyDL[,1] > preyCntRange[1] & datHuntVsPreyDL[,1] < preyCntRange[2],2 ]


nLL2=as.numeric(datSliceLL[seq(1,nDat,1)])
nNL2=as.numeric(datSliceNL[seq(1,nDat,1)])
nDL2=as.numeric(datSliceDL[seq(1,nDat,1)])

dataLL2=list(n=nLL2,NTOT=length(nLL2));
dataNL2=list(n=nNL2,NTOT=length(nNL2));
dataDL2=list(n=nDL2,NTOT=length(nDL2));

varnames1=c("q")
burn_in=1000;
steps=100000;
thin=10;

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

strPlotName <- paste("plots/stat_HuntEventRateParamPreyRange",preyCntRange[1],preyCntRange[2], ".pdf",sep="-")
pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 
X11()

hist(drawLL2$q[1,,1],breaks=seq(5,35,length=100),col=colourH[1],xlab="Rate Parameter",main=paste("Comparison using Poisson fit, to H.Events with  (",preyCntRange[1],"-",preyCntRange[2],") prey") )
hist(drawNL2$q[1,,1],breaks=seq(5,35,length=100),add=T,col=colourH[2])
hist(drawDL2$q[1,,1],breaks=seq(5,35,length=100),add=T,col=colourH[3])

legend(5,500,legend = c(paste("LL #",length(nLL2)),paste("NL #",length(nNL2)),paste("DL #",length(nDL2))),fill=colourH)

dev.off()
#hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
#hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
#hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

