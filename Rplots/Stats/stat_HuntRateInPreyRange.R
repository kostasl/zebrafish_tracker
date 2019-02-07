### Makes Bayssian inference on  hunt events rate and duration - within a given rannge of prey counts ##
## Used to plot spontaneous eye vergence events count in the EMPTY Test conditions ##
## Using an Exp distribution for Event Count, as
#  Poisson Distributions would be justified as a model of random occurance of huntevents through the 10 min recording time. Then at the group level the sum of individual poissons would still give
# render a poisson, thus modelling the group rate. 
# ** Yet the  empirical distribution does not look like poisson- but rather EXP like, heavy on the low event counts, and with a long tail. A Poisson For the group would imply an underlying 
# 
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

modelGEventRatePois="model { 
q ~ dnorm(15,0.0001)T(0,400)

for(j in 1:NTOT){
n[j] ~ dpois(q)
}


}"

modelGEventRateExp="model { 
q ~ dnorm(10,0.0001)T(0,100)

for(j in 1:NTOT){
n[j] ~ dexp(q)
}


}"



library(rjags)
strModelName = "modelGroupEventRate.tmp"
fileConn=file(strModelName)
writeLines(modelGEventRateExp,fileConn);
close(fileConn)


modelGEventDuration="model { 
s ~ dnorm(5,0.0001)T(0,100)
r ~ dnorm(5,0.0001)T(0,100)

for(j in 1:NTOT){
d[j] ~ dgamma(s,r)
}

}"

library(rjags)
strModelName = "modelGroupEventDuration.tmp"
fileConn=file(strModelName)
writeLines(modelGEventDuration,fileConn);
close(fileConn)



fromchain=1000

#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

colourH <- c(rgb(0.01,0.7,0.01,0.2),rgb(0.9,0.01,0.01,0.2),rgb(0.01,0.01,0.9,0.2),rgb(0.00,0.00,0.0,1.0))

strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged" ##Warning Set Includes Repeated Test For some LF fish - One In Different Food Density

message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSBMerged <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

##Remove Dublicates - Choose Labels - Duration Needs To be > 5ms
datHuntLabelledEventsSBMerged_filtered <- datHuntLabelledEventsSBMerged [
                with(datHuntLabelledEventsSBMerged, ( convertToScoreLabel(huntScore) != "Not_HuntMode/Delete" &
                                                                 convertToScoreLabel(huntScore) != "Duplicate/Overlapping" &
                                                                  (endFrame - startFrame) > 200 ) |  ## limit min event dur to 5ms
                                                                   eventID == 0), ] ## Add the 0 Event, In Case Larva Produced No Events
                                                                   
datHuntLabelledEventsSBMerged_filtered <- datHuntLabelledEventsSBMerged_filtered[!is.na(datHuntLabelledEventsSBMerged_filtered$groupID),]
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged_filtered)

## Get Event Counts Within Range  - Along With Total Number of Hunting frames for each Larva##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LL ) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LE ) )
datHuntVsPreyL <- datHuntVsPreyLL#rbind(datHuntVsPreyLL,datHuntVsPreyLE)
datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NL ) )
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NE ) )
datHuntVsPreyN <- datHuntVsPreyNL #rbind(datHuntVsPreyNL,datHuntVsPreyNE)
datHuntVsPreyN <- datHuntVsPreyN[!is.na(datHuntVsPreyN[,1]),]

datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DL ) )
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DE ) )
datHuntVsPreyD <-datHuntVsPreyDL #rbind(datHuntVsPreyDL,datHuntVsPreyDE)
datHuntVsPreyD <- datHuntVsPreyD[!is.na(datHuntVsPreyD[,1]),] ##Remove NA 


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)

##Larva Event Counts Slice
datSliceLF <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2], ]
datSliceNF <- datHuntVsPreyN[datHuntVsPreyN[,1] >= preyCntRange[1] & datHuntVsPreyN[,1] <= preyCntRange[2], ]
datSliceDF <- datHuntVsPreyD[datHuntVsPreyD[,1] >= preyCntRange[1] & datHuntVsPreyD[,1] <= preyCntRange[2], ]

##Select the event Count Column 2
nEventsLF=as.numeric(datSliceLF[,2])
nEventsNF=as.numeric(datSliceNF[,2])
nEventsDF=as.numeric(datSliceDF[,2])
nDatDF = length(nEventsDF);nDatNF = length(nEventsNF);nDatLF = length(nEventsLF);
dataLL2=list(n=nEventsLF,NTOT=nDatLF,food=as.integer(datSliceLF[,1]));
dataNL2=list(n=nEventsNF,NTOT=nDatNF,food=as.integer(datSliceNF[,1]));
dataDL2=list(n=nEventsDF,NTOT=nDatDF,food=as.integer(datSliceDF[,1]));

varnames1=c("n","q")
burn_in=1000;
steps=100000;
plotsamples = 10000
thin=2;
chains = 3


mLL2=jags.model(file="modelGroupEventRate.tmp",data=dataLL2,n.chains=chains);
mNL2=jags.model(file="modelGroupEventRate.tmp",data=dataNL2,n.chains=chains);
mDL2=jags.model(file="modelGroupEventRate.tmp",data=dataDL2,n.chains=chains);

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

hist(tail(drawLL2$q[1,,1],plotsamples),breaks=seq(0,1,length=100),col=colourH[1],xlab=" Event Rate (r)",
     main=paste("Spontaneous Eye Vergence Events")) #(",preyCntRange[1],"-",preyCntRange[2],") prey")
hist(tail(drawNL2$q[1,,1],plotsamples),breaks=seq(0,1,length=100),add=T,col=colourH[2])
hist(tail(drawDL2$q[1,,1],plotsamples),breaks=seq(0,1,length=100),add=T,col=colourH[3])

legend("topright",legend = c(paste("LF #",nDatLF),paste("NF #",nDatNF),paste("DF #",nDatDF)),fill=colourH)

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

legend("topleft",legend=paste(c("DF #","LF #","NF #"),c(nDatDF,nDatLF ,nDatNF ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

##Raw Histograms
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventCounts",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
Plim <- max(range(datHuntVsPreyL[,2])[2],range(datHuntVsPreyD[,2])[2],range(datHuntVsPreyN[,2])[2])
hist(datHuntVsPreyL[,2],breaks=seq(0,Plim,1),col=colourR[2],main="LE",xlab="",xlim=c(0,Plim),ylim=c(0,20))
hist(datHuntVsPreyN[,2],breaks=seq(0,Plim,1),col=colourR[3],main="NE",xlab="",xlim=c(0,Plim),ylim=c(0,20))
hist(datHuntVsPreyD[,2],breaks=seq(0,Plim,1),col=colourR[1],main="DE",xlab="# Eye Vergence Events ",xlim=c(0,Plim),ylim=c(0,20))
dev.off()

##Now Plot Infered Distributions
x <- seq(0,Plim,1)
N <- 1000
dev.off()

hist(datHuntVsPreyL[,2],breaks=seq(0,Plim,1),col=colourR[2],main="LE",xlab="",xlim=c(0,Plim),ylim=c(0,50))
for (r in tail(drawLL2$q[,,1],N))
  lines(x,100*dexp(x,r),type="p",pch=16,cex=0.4,main="LE",xlim=c(0,Plim),col=colourR[2] ) 
for (r in tail(drawDL2$q[,,1],N))
  lines(x,100*dexp(x,r),type="p",pch=16,cex=0.4,main="LE",xlim=c(0,Plim),col=colourR[1] ) 
for (r in tail(drawNL2$q[,,1],N))
  lines(x,100*dexp(x,r),type="p",pch=16,cex=0.4,main="NE",xlim=c(0,Plim),col=colourR[3] ) 



#### Also Plot Duration Of Eye Vergence Frames ###
###### Hunt Event Duration ############

## Raw Histograms for Total Hunt Duration Per Larva
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousTotalHuntDurationPerLarva",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHuntVsPreyL[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[2],main="LE",xlab="",xlim=c(0,40),ylim=c(0,40))
hist(datHuntVsPreyN[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[3],main="NE",xlab="",xlim=c(0,40),ylim=c(0,40))
hist(datHuntVsPreyD[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[1],main="DE",xlab=" Total Duration per Larva spent in Spontaneous Hunt Events ",xlim=c(0,40),ylim=c(0,40))
dev.off()
#### Also Plo

boxplot(datHuntVsPreyL[,3]/G_APPROXFPS,datHuntVsPreyD[,3]/G_APPROXFPS,datHuntVsPreyN[,3]/G_APPROXFPS,
        main="Spontaneous Hunt Events Duration per Larva ",notch=TRUE,names=c("LE","DE","NE"), ylab="(sec)" )

####### Function Returns Hunt Event Durations for Group ID, excluding events 0 (Food Count Event) 
getdatHuntDuration <- function(strGroupID)
{
  return (with(datHuntLabelledEventsSBMerged_filtered,
               data.frame(DurationFrames=endFrame[groupID == strGroupID & eventID != 0]-startFrame[groupID == strGroupID & eventID != 0],
                          expID=expID[groupID == strGroupID & eventID != 0]) ))
}

##Random Init Of Chain 
initDurfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(r=rgamma(N,7,0.5), ##Base Line Vergence Prior to HuntOn
                                     s=rgamma(N,7,0.5) ),
                        simplify=FALSE)
  
  return(initlist)
}


## Make DataStruct With Durations for each group - Exclude Event 0 ###
datHDuration_LE <- getdatHuntDuration("LE")
datHDuration_NE <- getdatHuntDuration("NE")
datHDuration_DE <- getdatHuntDuration("DE")

## Plot Histogram Of Durations in approx SEC
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHDuration_LE$DurationFrames /G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[2],xlab="",ylim=c(0,60),main="LE")
hist(datHDuration_NE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[3],main="NE",xlab="",ylim=c(0,60))
hist(datHDuration_DE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[1],main="DE",xlab="Duration of Spontaneous Hunt Events (sec) ",ylim=c(0,60))
dev.off()

## Baysian Inference Fiting a Gamma distribution to the Hunt Event Duration Data ##
##Setup Data Structure To Pass To RJAgs
varnames1=c("d","s","r")
burn_in=1000;
steps=10000;
plotsamples = 2000
thin=2;
chains = 3

datJagsLE=list(d=datHDuration_LE$DurationFrames,NTOT=NROW(datHDuration_LE));
datJagsNE=list(d=datHDuration_NE$DurationFrames,NTOT=NROW(datHDuration_NE));
datJagsDE=list(d=datHDuration_DE$DurationFrames,NTOT=NROW(datHDuration_DE));

mDurLE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsLE,n.chains=chains,inits=initDurfunct(chains,1));
mDurNE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsNE,n.chains=chains,inits=initDurfunct(chains,1));
mDurDE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsDE,n.chains=chains,inits=initDurfunct(chains,1));

update(mDurLE,burn_in)
update(mDurNE,burn_in)
update(mDurDE,burn_in)

drawDurLE=jags.samples(mDurLE,steps,thin=thin,variable.names=varnames1)
drawDurNE=jags.samples(mDurNE,steps,thin=thin,variable.names=varnames1)
drawDurDE=jags.samples(mDurDE,steps,thin=thin,variable.names=varnames1)

hist(drawDurLE$r,xlim=c(0,0.01))
hist(drawDurNE$r,xlim=c(0,0.01))
hist(drawDurDE$r,xlim=c(0,0.01))

layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
plot(drawDurLE$r,drawDurLE$s,xlim=c(0,0.01),ylim=c(0,8))
plot(drawDurNE$r,drawDurNE$s,xlim=c(0,0.01),ylim=c(0,8))
plot(drawDurDE$r,drawDurDE$s,xlim=c(0,0.01),ylim=c(0,8))

x <- seq(0,2700,1)
ns <- 1000
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
plot(x/G_APPROXFPS,dgamma(x,rate=drawDurLE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurLE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="LE",xlim=c(0,6),col=colourR[2] ) 
#hist(drawDurLE$d/G_APPROXFPS,breaks=100,xlim=c(0,6))
points(x/G_APPROXFPS,dgamma(x,rate=drawDurDE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurDE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="DE",xlim=c(0,6),col=colourR[1] ) 
points(x/G_APPROXFPS,dgamma(x,rate=drawDurNE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurNE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="NE",xlim=c(0,6),col=colourR[3] ) 
