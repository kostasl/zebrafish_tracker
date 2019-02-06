### Makes Bayssian inference on  hunt events rate - within a given rannge of prey counts ##
## Used to plot spontaneous eye vergence events count in the EMPTY Test conditions ##
### TODO Also Plot Duration Of Eye Vergence Frames ###


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


## Get Event Counts Within Range  - Along With Total Number of Hunting frames for each Larva##

datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LL ) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LE ) )
datHuntVsPreyL <- datHuntVsPreyLE#rbind(datHuntVsPreyLL,datHuntVsPreyLE)

datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NL ) )
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NE ) )
                         
datHuntVsPreyN <- datHuntVsPreyNE #rbind(datHuntVsPreyNL,datHuntVsPreyNE)

datHuntVsPreyN <- datHuntVsPreyN[!is.na(datHuntVsPreyN[,1]),]


datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DL ) )
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DE ) )
datHuntVsPreyD <-datHuntVsPreyDE #rbind(datHuntVsPreyDL,datHuntVsPreyDE)
##Remove NA 
datHuntVsPreyD <- datHuntVsPreyD[!is.na(datHuntVsPreyD[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,1)

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
plotsamples = 10000
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
strPlotName <- paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventRateParamPreyRange",preyCntRange[1],"-",preyCntRange[2], ".pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 

hist(tail(drawLL2$q[1,,1],plotsamples),breaks=seq(0,15,length=100),col=colourH[1],xlab=" Event Rate (r)",
     main=paste("Spontaneous Eye Vergence Events")) #(",preyCntRange[1],"-",preyCntRange[2],") prey")
hist(tail(drawNL2$q[1,,1],plotsamples),breaks=seq(0,15,length=100),add=T,col=colourH[2])
hist(tail(drawDL2$q[1,,1],plotsamples),breaks=seq(0,15,length=100),add=T,col=colourH[3])

legend("topright",legend = c(paste("LF #",nDatLL),paste("NF #",nDatNL),paste("DF #",nDatDL)),fill=colourH)

dev.off()


#### Plot Density ###
###Plot Density of Slope
dLLb<-density(tail(drawLL2$q[,,1],plotsamples)   )
dNLb<-density(tail(drawNL2$q[,,1] ,plotsamples)   )
dDLb<-density(tail(drawDL2$q[,,1],plotsamples)  )


pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDensityPreyRange",preyCntRange[1],"-",preyCntRange[2], ".pdf",sep=""))
plot(dDLb,col=colourL[1],lwd=3,lty=1, xlim=c(0.1,10), ylim=c(0,2),
     main="Rate of Spontaneous Eye Vergence Events ",
     xlab=expression(paste("Rate  ",lambda) ) )
lines(dLLb,col=colourL[2],xlim=c(0.5,1.2),lwd=3,lty=2)
lines(dNLb,col=colourL[3],xlim=c(0.5,1.2),lwd=3,lty=3)

legend("topleft",legend=paste(c("DF #","LF #","NF #"),c(nDatDL,nDatLL ,nDatNL ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

##Raw Histograms
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventCounts",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHuntVsPreyL[,2],breaks=seq(0,35,3),col=colourR[2],main="LE",xlab="",ylim=c(0,50))
hist(datHuntVsPreyN[,2],breaks=seq(0,35,3),col=colourR[3],main="NE",xlab="",ylim=c(0,50))
hist(datHuntVsPreyD[,2],breaks=seq(0,35,3),col=colourR[1],main="DE",xlab="# Eye Vergence Events ",ylim=c(0,50))
dev.off()
#### Also Plot Duration Of Eye Vergence Frames ###


##Raw Histograms
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHuntVsPreyL[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[2],main="LE",xlab="",ylim=c(0,50))
hist(datHuntVsPreyN[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[3],main="NE",xlab="",ylim=c(0,50))
hist(datHuntVsPreyD[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[1],main="DE",xlab="# Eye Vergence Events ",ylim=c(0,50))
dev.off()
#### Also Plo

boxplot(datHuntVsPreyL[,3]/G_APPROXFPS,datHuntVsPreyD[,3]/G_APPROXFPS,datHuntVsPreyN[,3]/G_APPROXFPS,
        main="Spontaneous Hunt Events Duration per Larva ",notch=TRUE,names=c("LE","DE","NE"), ylab="(sec)" )

##hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
##hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
##hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

### Do Linear Rate Vs Prey Model Plot -- Add Quantile Plots of Uncertainty in Parameters##
