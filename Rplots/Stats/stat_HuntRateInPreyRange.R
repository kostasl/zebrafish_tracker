### Makes Bayssian inference on  hunt events rate and duration - within a given rannge of prey counts ##
## Used to plot spontaneous eye vergence events count in the EMPTY Test conditions ##
## Using an Exp distribution for Event Count, as
#  Poisson Distributions would be justified as a model of random occurance of huntevents through the 10 min recording time. Then at the group level the sum of individual poissons would still give
# render a poisson, thus modelling the group rate. 
# ** Yet the  empirical distribution does not look like poisson- but rather EXP like, heavy on the low event counts, and with a long tail. A Poisson For the group would imply an underlying 
# 
### TODO Also Plot Duration Of Eye Vergence Frames ###

library(fitdistrplus)

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
q ~ dgamma(1,0.001)


for(j in 1:NTOT){
n[j] ~ dpois(q)
}


}"

modelGEventRateExp="model { 
q ~ dgamma(1,0.2) #SHape , Rate

for(j in 1:NTOT){
  n[j] ~  dexp(q)


}
}"

##For a Changing Rate in Low/ High Regime - 
modelGEventRateWeib="model { 
q ~ dnorm(10,0.0001)T(0,100)

for(j in 1:NTOT){
  n[j] ~  dweib(q,rate,shape)
  }
}"


## Discrete - Geometric Cause Mixture of rates - assuming rates drawn from most informative Prior distribution (EXP)
## Give geometric
modelGEventRateGeom="model { 
q ~ dunif(0.0,1)
r ~ dgamma(1,1)

for(j in 1:NTOT){
  n[j] ~  dnegbin(q,r) ##R=1 for Geometric
  }
}"



library(rjags)
strModelName = "modelGroupEventRate.tmp"
fileConn=file(strModelName)
writeLines(modelGEventRateGeom,fileConn);
close(fileConn)


#####
##Models Each Larva in the population of group Individually 
modelGEventDuration="model { 

for(j in 1:NTOT){
  d[j] ~ dgamma(s[hidx[j]],r[hidx[j]])
}

## Init Prior Per Larva ##
for (l in 1:max(hidx))
{
  s[l] ~ dnorm(5,0.0001)T(0,1000)
  r[l] ~ dnorm(5,0.0001)T(0,1000)
}

}"

library(rjags)
strModelName = "modelGroupEventDuration.tmp"
fileConn=file(strModelName)
writeLines(modelGEventDuration,fileConn);
close(fileConn)


## Run Baysian Inference on Model for Hunt Event Counts IN A Group/ Test Condition
## Return Samples Drawn structure
mcmc_drawEventCountModels <- function(datHuntVsPrey,preyCountRange,strModelFilename)
{
  varnames1=c("n","q","r")
  burn_in=1000;
  steps=100000;
  plotsamples = 10000
  thin=2;
  chains = 3
  
  
  ##Larva Event Counts Slice
  datSliceF <- datHuntVsPrey[datHuntVsPrey[,1] >= preyCntRange[1] & datHuntVsPrey[,1] <= preyCntRange[2], ]
  
  ##Select the event Count Column 2
  nEvents=as.numeric(datSliceF[,2])
  nDat = length(nEvents);
  dataG=list(n=nEvents,NTOT=nDat,food=as.integer(datSliceF[,1]));
  
  model=jags.model(file=strModelFilename,data=dataG,n.chains=chains);
  
  update(model,burn_in)
  
  drawSamples=jags.samples(model,steps,thin=thin,variable.names=varnames1)
  
  return(drawSamples) 
}



## Combines a the Gamma Fit distributions with a histogram and an empirical density estimate of the data (with a given BW)
## datHDuration the raw Duration Dataframe, drawDur the MCMC Sampled values for the gamma distr. params
## The histogram is scaled down fit the range of 
plotDurationDensityFitComparison <- function(datHDuration,drawDur,lcolour,HLim,nplotSamples)
{
  XLim <- 6
  
  x <- seq(0,HLim,1) ##In Seconds
  N <- 1000
  
  
  densDur<-density(datHDuration[,1],bw=pBW)   
  histDur <- hist(datHDuration[,1]/G_APPROXFPS,breaks=seq(0,HLim/G_APPROXFPS,0.1),plot=FALSE)
  
  
  YScale <- range(densDur$y)[2]*1.05
  hist_scale <- max(histDur$counts)/YScale ##Scale Histogram Relative To Density Peak
  YLim <- 0.0015
  
  par(cex.lab = 1.1)
  par(cex = 0.8)
  par(cex.axis = 1.1)
  
  par(mar = c(5,5,2,5))
  plot(densDur$x/G_APPROXFPS, densDur$y,type='l',xlim=c(0,XLim),ylim=c(0,YLim),lty=2,lwd=4,ylab=NA,xlab=NA)
  
  
  pBW <- 120 ## Estimation BandWidth
  ## Live Fed
  for (c in 1:NROW(drawDur$r[1,1,]) )
    for (i in (NROW(drawDur$r[,,c])-nplotSamples):NROW(drawDur$r[,,c]) )
      lines(x/G_APPROXFPS,dgamma(x,rate=drawDur$r[,i,c],shape=drawDur$s[,i,c]),type="p",pch=16,cex=1.4,xlim=c(0,XLim),ylim=c(0,YLim), col=lcolour ) 
  
  par(new=T)
  plot(histDur$breaks[1:NROW(histDur$counts)],histDur$counts, axes=F, xlab=NA, ylab=NA,cex=1.1,
       xlim=c(0.0,XLim),ylim=c(0,max(histDur$counts)*1.10 ) ,lwd=2,pch=21,col="#000000CC")
  axis(side = 4)
  mtext(side = 4, line = 2.1, 'Counts')
  mtext(side = 2, line = 2.1, 'P(s)')
  
  
} ##Plot Function


## Compare Model TO Data Using CDF ##
plotEventCountDistribution_cdf <- function(datHEventCount,drawHEvent,lcolour,lpch,lty,Plim,nplotSamples=100,newPlot = FALSE)
{
  XLim <- 50
  x <- seq(0,XLim,1)

  cdfD_N <- ecdf(datHEventCount[,2])

  plot(cdfD_N,col=lcolour,pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim),ylim=c(0,1),cex=1.5,cex.lab=1.5,add=!newPlot)
  ##Construct CDF of Model by Sampling randomly from Model distribution for exp rate parameter
  for (c in 1:NROW(drawHEvent$q[1,1,])) {
    for (q in  tail(drawHEvent$q[,,c],nplotSamples) )
    {
      cdfM <- dnbinom(x,size=1,prob=q)##1-exp(-q*x) ##ecdf(  dexp( x, q  ) )
      lines(cumsum(cdfM),col=lcolour,lty=lty) #add=TRUE,
    }
  }
  plot(cdfD_N,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim),ylim=c(0,1),cex=1.5,cex.lab=1.5,add=TRUE)
  #axis(side = 4)
  #mtext(side = 4, line = 2.1, 'Counts')
  ### Draw Distribution oF Hunt Rates - 
  ##(z= p/(1-p))
#  hist( (1-tail(drawHEvent$q[,,c],nplotSamples))/tail(drawHEvent$q[,,c],nplotSamples)  )
}


## Plot the Histogram, Empirical Density And Regressor Fit For A group Showing DENSITY VS Histogram #
plotEventCountDistribution_hist <- function(datHEventCount,drawHEvent,lcolour,HLim,nplotSamples)
{
  XLim <- 15
  x <- seq(0,Plim,0.1)
  N <- nplotSamples
  pBW <-1.2
  
  
  densHEventS <-density(datHEventCount[,2],bw=pBW)   
  histHEvent <- hist(datHEventCount[,2],breaks=seq(0,Plim,1),plot=FALSE)
  
  YLim <- 0.25 ##range(densHEventS$y)[2]*1.05
  par(mar = c(5,5,2,5))
  
  plot(densHEventS$x, densHEventS$y,type='l',xlim=c(0,XLim),ylim=c(0,YLim*1.20),xlab=NA, ylab=NA,lty=2,lwd=4)
  for (c in 1:NROW(drawHEvent$q[1,1,])) 
    for (r in (NROW(drawHEvent$q[,,c])-N):NROW(drawHEvent$q[,,c]))
      lines(x,dexp(x,rate=drawHEvent$q[,r,c] ),type="l",pch=16,cex=0.4,xlim=c(0,HLim),col=lcolour) 
  
  par(new=T)

  plot(histHEvent$breaks[1:NROW(histHEvent$counts)],histHEvent$counts, cex=1.1,ylab=NA,xlab=NA,axes=F,
       xlim=c(0.0,XLim),ylim=c(0,30 ) ,lwd=2,pch=21,col="#000000CC")
  
  axis(side = 4)
  mtext(side = 4, line = 2.1, 'Counts')
  mtext(side = 2, line = 2.1, 'P(s)')
  
}

fromchain=1000

#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

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
## Added Larva ID to Check for Correlation Through Time of Day - Surrogate as LarvaID;s increased through the day of the experiment from 1-4
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LL ),datHuntStat[,"vIDLookupTable"]$LL$larvaID )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LE ),datHuntStat[,"vIDLookupTable"]$LE$larvaID  )
datHuntVsPreyL <- datHuntVsPreyLL#rbind(datHuntVsPreyLL,datHuntVsPreyLE)
datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]) & datHuntVsPreyL[,2] < 25,]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NL),datHuntStat[,"vIDLookupTable"]$NL$larvaID )
datHuntVsPreyNL <- datHuntVsPreyNL[!is.na(datHuntVsPreyN[,1]) & datHuntVsPreyNL[,2] < 25,]
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NE),datHuntStat[,"vIDLookupTable"]$NE$larvaID  )
datHuntVsPreyNE <- datHuntVsPreyNE[!is.na(datHuntVsPreyN[,1]) & datHuntVsPreyNE[,2] < 25,]

datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DL ),datHuntStat[,"vIDLookupTable"]$DL$larvaID  )
datHuntVsPreyDL <- datHuntVsPreyDL[!is.na(datHuntVsPreyDL[,1]) & datHuntVsPreyDL[,2] < 25,] ##Remove NA And High Fliers
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DE ),datHuntStat[,"vIDLookupTable"]$DE$larvaID  )
datHuntVsPreyDE <- datHuntVsPreyDE[!is.na(datHuntVsPreyDE[,1]) & datHuntVsPreyDE[,2] < 25,] ##Remove NA And High Fliers


## Check Number of Hunt Events For A LarvaID 
vEventCount <- vector()
datHuntVsPreyC <- datHuntVsPreyL ##Check Which Group
for (i in 1:4)
{
  vEventCount[i] <- ( sum(datHuntVsPreyC[ datHuntVsPreyC[,4] == i & !is.na(datHuntVsPreyC[,4]) ,2] ) )
}
print(vEventCount)

### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)


drawLL2 <- mcmc_drawEventCountModels(datHuntVsPreyLL,preyCntRange,"modelGroupEventRate.tmp")
drawNL2 <- mcmc_drawEventCountModels(datHuntVsPreyNL,preyCntRange,"modelGroupEventRate.tmp")
drawDL2 <- mcmc_drawEventCountModels(datHuntVsPreyDL,preyCntRange,"modelGroupEventRate.tmp")
drawLE2 <- mcmc_drawEventCountModels(datHuntVsPreyLE,preyCntRange,"modelGroupEventRate.tmp")
drawNE2 <- mcmc_drawEventCountModels(datHuntVsPreyNE,preyCntRange,"modelGroupEventRate.tmp")
drawDE2 <- mcmc_drawEventCountModels(datHuntVsPreyDE,preyCntRange,"modelGroupEventRate.tmp")

#### DEBUG CODE - Fitting An EXP Distribution TEST ###
f1L <- fitdist( datHuntVsPreyLL[,2],"exp")
plot(f1L,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f1D <- fitdist( datHuntVsPreyDL[,2],"exp");
plot(f1D,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f1N <- fitdist( datHuntVsPreyNL[,2],"exp");
plot(f1N,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
##### ### 
f2D <- fitdist( datHuntVsPreyDL[,2],"nbinom",lower = c(0, 0))
plot(f2D) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f2N <- fitdist( datHuntVsPreyNL[,2],"nbinom",lower = c(0, 0)) ##,start = list(scale = 1, shape = 1)
plot(f2N) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f2L <- fitdist( datHuntVsPreyLL[,2],"nbinom",lower = c(0, 0)) #,start = list(scale = 1, shape = 1)
plot(f2L) #,start = list(scale = 1, shape = 1) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
## ########


#### HUNT EVENT PER LARVA PLOT #####
## Comprehensive Plot On Number of Hunt Events
pdf(file= paste(strPlotExportPath,"/stat/stat_LiveFoodTestHuntEventCounts",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
##Now Plot Infered Distributions
Plim <- max(range(datHuntVsPreyLL[,2])[2],range(datHuntVsPreyDL[,2])[2],range(datHuntVsPreyNL[,2])[2])
x <- seq(0,Plim,0.1)
##Show Alignment with Empirical Distribution of HuntEvent Numbers
## Number of Hunt Events Per Larva
plotsamples <- 70
layout(matrix(c(1,2,3,4,5,6), 3,2, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,3.01,1,1))
lineTypeL[1] <- 1
plotEventCountDistribution_cdf(datHuntVsPreyNE,drawNE2,colourHE[1],pchL[1],lineTypeL[1],Plim,plotsamples,newPlot=TRUE )
plotEventCountDistribution_cdf(datHuntVsPreyNL,drawNL2,colourHL[1],pchL[3],lineTypeL[1],Plim,plotsamples,newPlot=FALSE )
legend("bottomright",legend = c(paste("Data NE #",NROW(datHuntVsPreyNE) ),paste("Model NE "),
                                paste("Data NL #",NROW(datHuntVsPreyNL) ),paste("Model NL ")), 
       col=c(colourP[4], colourLegE[1],colourP[4],colourLegL[1]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
plotEventCountDistribution_cdf(datHuntVsPreyLE,drawLE2,colourHE[2],pchL[1],lineTypeL[1],Plim,plotsamples,newPlot=TRUE )
plotEventCountDistribution_cdf(datHuntVsPreyLL,drawLL2,colourHL[2],pchL[3],lineTypeL[1],Plim,plotsamples,newPlot=FALSE )
legend("bottomright",legend = c(paste("Data LE #",NROW(datHuntVsPreyLE) ),paste("Model LE "),
                                paste("Data LL #",NROW(datHuntVsPreyLL) ),paste("Model LL ")), 
      col=c(colourP[4], colourLegE[2],colourP[4],colourLegL[2]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
plotEventCountDistribution_cdf(datHuntVsPreyDE,drawDE2,colourHE[3],pchL[1],lineTypeL[1],Plim,plotsamples,newPlot=TRUE  )
plotEventCountDistribution_cdf(datHuntVsPreyDL,drawDL2,colourHL[3],pchL[3],lineTypeL[1],Plim,plotsamples,newPlot=FALSE  )
legend("bottomright",legend = c(paste("Data DE #",NROW(datHuntVsPreyLE) ),paste("Model DE "),
                                paste("Data DL #",NROW(datHuntVsPreyLL) ),paste("Model DL ")), 
       col=c(colourP[4], colourLegE[3],colourP[4],colourLegL[3]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
mtext(side = 1,cex=0.8, line = 2.2, "Hunt Event Counts (N)")
mtext(side = 2,cex=0.8, line = 2.2, " F(x < N) ")

## Compare Distrib Of Model Params ##
schain <- 1:3
Ylim <- 50 
pBW <- 0.001
densHEventSampled_LL <-density(drawLL2$q[1,,schain],bw=pBW);densHEventSampled_DL <-density(drawDL2$q[1,,schain],bw=pBW);densHEventSampled_NL <-density(drawNL2$q[1,,schain],bw=pBW)
densHEventSampled_LE <-density(drawLL2$q[1,,schain],bw=pBW);densHEventSampled_DE <-density(drawDL2$q[1,,schain],bw=pBW);densHEventSampled_NE <-density(drawNL2$q[1,,schain],bw=pBW)   

plot(densHEventSampled_NL$x, densHEventSampled_NL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[1],col=colourL[1],lwd=4,ylab=NA,xlab=NA)
lines(densHEventSampled_LL$x, densHEventSampled_LL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[2],col=colourL[2],lwd=4,ylab=NA,xlab=NA)
lines(densHEventSampled_DL$x, densHEventSampled_DL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[3],col=colourL[3],lwd=4,ylab=NA,xlab=NA)
mtext(side = 1,cex=0.8, line =2.2, expression(paste("Model  Parameter (p",") ") ) )
mtext(side = 2,cex=0.8, line =2.2, paste("Density Bw:",densHEventSampled_LL$bw))
legend("topright",legend = c(paste("NF " ),paste("LF ") ,paste("DF ") ) ,
       col=colourL,lty=lineTypeL,lwd=3,seg.len=3,cex=1.1)

### Draw Distribution oF Hunt Rates - 
##(z= p/(1-p))
densHEventHuntRate_LE <-density((1-tail(drawLE2$q[,,schain],nplotSamples))/tail(drawLE2$q[,,schain],nplotSamples));
densHEventHuntRate_DE <- density((1-tail(drawDE2$q[,,schain],nplotSamples))/tail(drawDE2$q[,,schain],nplotSamples));   
densHEventHuntRate_NE <- density((1-tail(drawNE2$q[,,schain],nplotSamples))/tail(drawNE2$q[,,schain],nplotSamples));   

plot(densHEventHuntRate_LE$x, densHEventHuntRate_LE$y,type='l',xlim=c(0,12),ylim=c(0,1),lty=lineTypeL[2],col=colourL[2],lwd=4,ylab=NA,xlab=NA)
lines(densHEventHuntRate_NE$x,densHEventHuntRate_NE$y,type='l',xlim=c(0,12),ylim=c(0,1),lty=lineTypeL[1],col=colourL[1],lwd=4,ylab=NA,xlab=NA)
lines(densHEventHuntRate_DE$x, densHEventHuntRate_DE$y,type='l',xlim=c(0,12),ylim=c(0,1),lty=lineTypeL[3],col=colourL[3],lwd=4,ylab=NA,xlab=NA)
legend("topright",legend = c(paste("NF " ),paste("LF ") ,paste("DF ") ) ,
       col=colourL,lty=lineTypeL,lwd=3,seg.len=3,cex=1.1)
mtext(side = 1,cex=0.8, line =2.2, expression(paste(" Hunt Rate  (",lambda,") ") ) )

##BoxPlot
boxplot(datHuntVsPreyN[,2],datHuntVsPreyL[,2] ,datHuntVsPreyD[,2],
        main=NA,notch=TRUE,col=colourD,names=c("NF","LF","DF"),ylim=c(0,50)  )
mtext(side = 2,cex=0.8, line =2.2, "Hunt Event Counts (N)")

dev.off() 
################## ############# ### # # 


############## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ########################
#### Plot Duration Of Eye Vergence Frames Per Larva - Time Spent Hunting Per Larva ###
## Note: I do not have enough data points per larvae, so as to make a model for both 
## the Duration Per Larva, and the overall group - So best to do them separatelly 
## Statistics of overall duration per larva, and statistic of individual 
######                Hunt Event Duration                                 ############

####### Function Returns Hunt Event Durations for Group ID, excluding events 0 (Food Count Event) 
getdatHuntDuration <- function(strGroupID)
{
  
  datDurationPerEpisodePerLarva <- (with(datHuntLabelledEventsSBMerged_filtered,
               data.frame(DurationFrames=endFrame[groupID == strGroupID & eventID != 0]-startFrame[groupID == strGroupID & eventID != 0],
                          expID=expID[groupID == strGroupID & eventID != 0],
                          hidx=as.numeric(factor(expID[groupID == strGroupID & eventID != 0]))) )) ##Add hidx to use for correct prior Init in RJags
  
  datDurationPerEpisodePerLarva
  
  return (datDurationPerEpisodePerLarva)
}

##Random Init Of Chain 
initDurfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(r=abs(rnorm(N,3,1)), ##Base Line Vergence Prior to HuntOn
                                     s=abs(rnorm(N,3,1)) ),
                        simplify=FALSE)
  
  return(initlist)
}


## Make DataStruct With Durations for each group - Exclude Event 0 ###
datHDuration_LE <- getdatHuntDuration("LE")
datHDuration_NE <- getdatHuntDuration("NE")
datHDuration_DE <- getdatHuntDuration("DE")

datHDuration_L <- datHDuration_LE
datHDuration_N <- datHDuration_NE
datHDuration_D <- datHDuration_DE

## Baysian Inference Fiting a Gamma distribution to the Hunt Event Duration Data ##
##Setup Data Structure To Pass To RJAgs
varnames1=c("d","s","r","hidx")
burn_in=1000;
steps=10000;
plotsamples = 2000
thin=2;
chains = 3

datJagsLE=list(d=datHDuration_LE$DurationFrames/G_APPROXFPS,NTOT=NROW(datHDuration_LE),hidx=datHDuration_LE$hidx);
datJagsNE=list(d=datHDuration_NE$DurationFrames/G_APPROXFPS,NTOT=NROW(datHDuration_NE),hidx=datHDuration_NE$hidx);
datJagsDE=list(d=datHDuration_DE$DurationFrames/G_APPROXFPS,NTOT=NROW(datHDuration_DE),hidx=datHDuration_DE$hidx);

mDurLE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsLE,n.chains=chains,inits=initDurfunct(chains,max(unique(datJagsLE$hidx) )));
mDurNE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsNE,n.chains=chains,inits=initDurfunct(chains,max(unique(datJagsNE$hidx) )));
mDurDE=jags.model(file="modelGroupEventDuration.tmp",data=datJagsDE,n.chains=chains,inits=initDurfunct(chains,max(unique(datJagsDE$hidx) )));

update(mDurLE,burn_in)
update(mDurNE,burn_in)
update(mDurDE,burn_in)

drawDurL=jags.samples(mDurLE,steps,thin=thin,variable.names=varnames1)
drawDurN=jags.samples(mDurNE,steps,thin=thin,variable.names=varnames1)
drawDurD=jags.samples(mDurDE,steps,thin=thin,variable.names=varnames1)

layout(matrix(c(1,2,3,4,5,6), 3,2, byrow = FALSE))
hist(drawDurL$r)
hist(drawDurN$r)
hist(drawDurD$r)

hist(drawDurL$s)
hist(drawDurN$s)
hist(drawDurD$s)

plot(tail(drawDurL$r,1000), tail(drawDurL$s,1000),ylim=c(0,200))
plot(tail(drawDurD$r,1000), tail(drawDurD$s,1000),ylim=c(0,200))
plot(tail(drawDurN$r,1000), tail(drawDurN$s,1000),ylim=c(0,200))

######## Plot Comparison Of Duration Data ##########
##Now Plot Infered Distributions


## Plot Histogram Of Durations in approx SEC
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHDuration_LE$DurationFrames /G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[2],xlab="",ylim=c(0,60),main="LE")
hist(datHDuration_NE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[3],main="NE",xlab="",ylim=c(0,60))
hist(datHDuration_DE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[1],main="DE",xlab="Duration of Spontaneous Hunt Events (sec) ",ylim=c(0,60))
dev.off()

## Raw Histograms for Total Hunt Duration Per Larva
pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousTotalHuntDurationPerLarva",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
hist(datHuntVsPreyL[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[2],main="LE",xlab="",xlim=c(0,40),ylim=c(0,40))
hist(datHuntVsPreyN[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[3],main="NE",xlab="",xlim=c(0,40),ylim=c(0,40))
hist(datHuntVsPreyD[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[1],main="DE",xlab=" Total Duration per Larva spent in Spontaneous Hunt Events ",xlim=c(0,40),ylim=c(0,40))
dev.off()
#### Also Plo

## Box plot Of Total Duration Per Larva
boxplot(datHuntVsPreyN[,3]/G_APPROXFPS,datHuntVsPreyL[,3]/G_APPROXFPS,datHuntVsPreyD[,3]/G_APPROXFPS,
        main=NA,notch=TRUE,names=c("NE","LE","DE"),ylim=c(0,40), ylab="(sec)",col=colourH )




#### Plot Density ###
###Plot Density of Slope

### Make Plot Of Histograms and Gamma Fit 
strPlotName <- paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration_p",preyCntRange[1],"-",preyCntRange[2], ".pdf",sep="")
pdf(strPlotName,width=16,height=10,title="Comparing the Duration of spontaneous hunt events " ) 


Plim <- max( round(range(datHDuration_L[,1])[2] ),round(range(datHDuration_D[,1])[2] ),round(range(datHDuration_N[,1])[2] ))*1.1   
layout(matrix(c(1,2,3,4,4,5), 3,2, byrow = FALSE))
plotDurationDensityFitComparison(datHDuration_L,drawDurL,colourH[2],Plim,10) #colourR[2]
legend("topright",legend=paste(c("Empirical Density ","Baysian Gamma Fit","Histogram ") )
       ,pch=c(NA,NA,21),lwd=c(2,1,2),lty=c(2,1,NA),col=c("black",colourL[2],colourH[4]) )

plotDurationDensityFitComparison(datHDuration_D,drawDurD,colourH[1],Plim,10)
plotDurationDensityFitComparison(datHDuration_N,drawDurN,colourH[3],Plim,10)
 ## Plot Distrib Of Params
ns <- 200
plot(drawDurL$r[,(steps/thin-ns):(steps/thin),],drawDurL$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[2],xlim=c(0,0.01),ylim=c(0,6),
     xlab="Rate Parameter",ylab="Shape Parameter") 
points(drawDurD$r[,(steps/thin-ns):(steps/thin),],drawDurD$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[1] ) 
points(drawDurN$r[,(steps/thin-ns):(steps/thin),],drawDurN$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[3] ) 

legend("topright",legend=paste(c("DF #","LF #","NF #"),c(nDatDF,nDatLF ,nDatNF ) )
       ,pch=16,col=colourL)

boxplot(datHDuration_L$DurationFrames/G_APPROXFPS,datHDuration_D$DurationFrames/G_APPROXFPS,datHDuration_N$DurationFrames/G_APPROXFPS,
        main="Hunt Duration per Larva ",notch=TRUE,names=c("LE","DE","NE"),ylim=c(0,6), ylab="(sec)" )


dev.off()

#points(x/G_APPROXFPS,dgamma(x,rate=drawDurDE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurDE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="DE",xlim=c(0,6),col=colourR[1] ) 
#points(x/G_APPROXFPS,dgamma(x,rate=drawDurNE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurNE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="NE",xlim=c(0,6),col=colourR[3] ) 
