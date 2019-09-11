### KOstasl 2019

### Produces Figure that compares the statistics of total Larva HUnt Duration for each experiment in a group ###
### 
source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")



## Apply the Same Logic As Counting Hunt Events, to Counting Hunt Frames : ie Hunt Frames occur with some Rate Lambda, drawn from 
## From a (set of) Gamma Distribution 
## Discrete - Geometric Cause Mixture of rates - assuming rates drawn from most informative Prior distribution (EXP)
## Assuming nbinom(r,p) Poisson(L|a,b) Gamma(a,b) then r=a, p=1/(b+1) -> b=(1-p)/p
## Give Neg Binomial
modelLarvaHuntDuration="model { 
q ~ dunif(0.0,1)
r ~ dgamma(1,1) #dgamma(r - shape , lambda-Rate) λ^rx^(r−1)*exp(−λx)/Gamma(r)

for(j in 1:NTOT){
  d[j] ~  dnegbin(q,r) ##Number Of Hunt Frames Per Larva
  #d[j] ~ dweib(r, mu) ##dweib(v, lambda)
  }
}"


library(rjags)
strModelName = "modelLarvaHuntDuration.tmp"
fileConn=file(strModelName)
writeLines(modelLarvaHuntDuration,fileConn);
close(fileConn)

###Alternativelly We May Probe the duration of individual Events ... 
##Models Each Event in the population of larvaer Individually 
modelEventDuration="model { 

for(j in 1:NTOT){
d[j] ~ dgamma(s[hidx[j]],r[hidx[j]]) #dgamma( r-shape , lambda-Rate) λ^rx^(r−1)*exp(−λx)/Gamma(r)
}

## Init Prior Per Larva ##
for (l in 1:max(hidx))
{
  s[l] ~ dnorm(500,0.0001)T(0,1000)
  r[l] ~ dnorm(1,0.0001)T(0,1000)
}

}"


library(rjags)
strModelName = "modelLarvaEventDuration.tmp"
fileConn=file(strModelName)
writeLines(modelEventDuration,fileConn);
close(fileConn)


############## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ########################
#### Plot Duration Of Eye Vergence Frames Per Larva - Time Spent Hunting Per Larva ###
## Note: I do not have enough data points per larvae, so as to make a model for both 
## the Duration Per Larva, and the overall group - So best to do them separatelly 
## Statistics of overall duration per larva, and statistic of individual 
######                Hunt Event Duration                                 ############

## Run Baysian Inference on Model for Hunt Event Counts IN A Group/ Test Condition
## Return Samples Drawn structure
mcmc_drawHuntDurationModels <- function(datHuntVsPrey,preyCountRange,strModelFilename)
{
  varnames1=c("d","q","r")
  burn_in=1000;
  steps=100000;
  plotsamples = 10000
  thin=2;
  chains = 3
  
  
  ##Larva Event Counts Slice
  datSliceF <- datHuntVsPrey[datHuntVsPrey[,1] >= preyCntRange[1] & datHuntVsPrey[,1] <= preyCntRange[2], ]
  
  datJags=list(d=datHuntVsPrey[,3],NTOT=NROW(datHuntVsPrey));
  
  nDat = NROW(datHuntVsPrey);
  dataG=list(d=datHuntVsPrey[,3],NTOT=nDat,food=as.integer(datSliceF[,1]));
  
  model=jags.model(file=strModelFilename,data=dataG,n.chains=chains);
  
  update(model,burn_in)
  
  drawSamples=jags.samples(model,steps,thin=thin,variable.names=varnames1)
  
  return(drawSamples) 
}



## Compare Model TO Data Using CDF ##
## Comverting frames to Duration using Approx FPS - 
## Note : plot at intervals otherwise plot is slow and very large
plotHuntDurationDistribution_cdf <- function(datHDuration,drawHEvent,lcolour,lpch,lty,Plim,nplotSamples=100,newPlot = FALSE)
{
  XLim <- G_APPROXFPS*180
  XStep <- round(G_APPROXFPS)
  x <- seq(XStep,XLim,XStep)
  x_s <- seq(0,XLim,1)
  yDat <- datHDuration[,3]/G_APPROXFPS
  
  cdfD_N <- ecdf(yDat)
  
  plot(cdfD_N,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0, XLim/G_APPROXFPS),
       ylim=c(0,1),cex=cex,cex.axis=cex,cex.lab=cex,add=!newPlot)
  ##Construct CDF of Model by Sampling randomly from Model distribution for exp rate parameter
  for (c in 1:NROW(drawHEvent$q[1,1,])) {
    for (j in (NROW(drawHEvent$q[,,c])-nplotSamples):NROW(drawHEvent$q[,,c]) )
    {##seq(0,XLim,1)
     ## Step-wise approx
      #cdfM <- cumsum((XStep  * dnbinom(x,size=drawHEvent$r[,j,c], prob=  drawHEvent$q[,j,c]  ) ) )##1-exp(-q*x) ##ecdf(  dexp( x, q  ) )
      #lines(x/G_APPROXFPS,cdfM,col=lcolour,lty=lty) #add=TRUE,
      
      ##Complete on per frame duration data - cummulative distribution
      cdfM <- cumsum(dnbinom(x_s,size=drawHEvent$r[,j,c], prob=  drawHEvent$q[,j,c]  ))##1-exp(-q*x) ##ecdf(  dexp( x, q  ) )
      ##Sample CDF at key points 
      lines(x[1:NROW(cdfM[x])]/G_APPROXFPS,cdfM[x],col=lcolour,lty=lty) #add=TRUE,
    }
  }
  plot(cdfD_N ,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim/G_APPROXFPS),ylim=c(0,1),cex=cex,cex.axis=cex,cex.lab=cex,add=TRUE)
  
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
      lines(x/G_APPROXFPS,dgamma(x,rate=drawDur$r[,i,c],shape=drawDur$s[,i,c]),
            type="p",pch=16,cex=1.4,xlim=c(0,XLim),ylim=c(0,YLim), col=lcolour ) 
  
  par(new=T)
  plot(histDur$breaks[1:NROW(histDur$counts)],histDur$counts, axes=F, xlab=NA, ylab=NA,cex=cex,
       xlim=c(0.0,XLim),ylim=c(0,max(histDur$counts)*1.10 ) ,lwd=2,pch=21,col="#000000CC")
  axis(side = 4)
  mtext(side = 4, cex=cex,line = 2.1, 'Counts')
  mtext(side = 2,cex=cex, line = 2.1, 'P(s)')
  
  
} ##Plot Function


## Box Plot Assist To Connect Individual Hunt Event Counts Between Empty (Spontaneous) and Live Test Conditions (Evoked) 
## For Each Larva Experiment
plotConnectedHuntDuration <- function(datHuntStat,vDat,strCondTags)
{
  
  ### Plot Connected Larva Event Counts - To Show Individual Behaviour In Spontaneous Vs Evoked Activity
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    vRegL <- datHuntStat[,"vIDLookupTable"][[gL]]
    vRegE <- datHuntStat[,"vIDLookupTable"][[gE]]
    
    ## Fix Missing LarvaID: REMOVE FActor Field / Set NAs Which Are really ID 5
    vRegL[,"larvaID"] <- as.numeric(vRegL[,"larvaID"])
    vRegE[,"larvaID"] <- as.numeric(vRegE[,"larvaID"])
    vRegL[is.na(vRegL$larvaID),"larvaID"] <- 5 
    vRegE[is.na(vRegE$larvaID),"larvaID"] <- 5
    
    ## Drop Factors To Skip Errors 
    vRegL[,"dataSetID"] <- as.numeric(vRegL[,"dataSetID"])
    vRegE[,"dataSetID"] <- as.numeric(vRegE[,"dataSetID"])
    
    datSetID <- levels(vIDTable[[gE]]$dataSetID) #[ vIDTable[[gE]]$dataSetID[  vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID &
    #vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID] ]
    for (k in 1:NROW(vRegE) )
    {
      e <- vRegE[k,]
      
      ptSrc  <- vDat[[gE]][which(names(vDat[[gE]]) == e$expID) ]##Get Idx for Event Count from Specificed Experiment, and retrieve Event Count At Empty 
      ##Find the Live Test for this Larva, Via the Larva ID dataset ID combination
      LivePairExp <-(vRegL[vRegL$dataSetID == e$dataSetID & vRegL$larvaID == e$larvaID,])
      if (NROW(LivePairExp) == 0)
        next() ##Skip If Matched LiveTest Larva Is not Found
      
      message(e$expID)
      
      ptDest <- vDat[[gL]][which(names(vDat[[gL]]) %in% LivePairExp$expID) ]
      
      
      ##Plot The Lines Connect Each Empty Tested Larva With Itself In THe Live Fed Conditions 
      #Do not Plot NA connecting line
      idxSrc  <- match(gE,strCondTags) ##Bar Center Idx for Each Condition E. Fed
      idxDest <- match(gL ,strCondTags)           
      
      points(idxSrc,log10(ptSrc+1),pch=1,
             col=colourP[4],cex=1 )
      for (pIDest in ptDest) ##Sometimes there are Multiple Live Tests for an Empty One
        points(idxDest,log10(pIDest+1),pch=1, col=colourP[4],cex=1 )
      
      segments(gIdx,log10(ptSrc+1),gIdx+1,log10(ptDest+1), col=colourP[4])
      
      
    }##For Each Experiment
    
  } ## Go Through Pairs Of Conditions ##
}##End Of Function ConnectEventPoints



##Random Init Of Chain ## For nChains, and N priors 
initLarvaHuntDurfunct <- function(nchains,N)
{
  
  initlist <- replicate(nchains,list(r=abs(rgamma(N,100,0.1)), ##Base Line Vergence Prior to HuntOn
                                     q=abs(runif(N)) ),
                        simplify=FALSE)
  
  return(initlist)
}


####### Function Returns Hunt Event Durations for Group ID, excluding events 0 (Food Count Event) 
getHuntEventDuration <- function(datLabelled,strGroupID)
{
  
  datDurationPerEpisodePerLarva <- (with(datLabelled,
                                         data.frame(DurationFrames=endFrame[groupID == strGroupID & eventID != 0]-startFrame[groupID == strGroupID & eventID != 0],
                                                    expID=expID[groupID == strGroupID & eventID != 0],
                                                    hidx=as.numeric(factor(expID[groupID == strGroupID & eventID != 0]))) )) ##Add hidx to use for correct prior Init in RJags
  
  datDurationPerEpisodePerLarva
  
  return (datDurationPerEpisodePerLarva)
}


##Random Init Of Chain For nChains, and N priors 
initEventDurfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(r=abs(rnorm(N,3,1)), ##Base Line Vergence Prior to HuntOn
                                     s=abs(rnorm(N,3,1)) ),
                        simplify=FALSE)
  
  return(initlist)
}


mcmc_drawEventDurationModels <- function(datHuntVsPrey,preyCountRange,strModelFilename)
{
  varnames1=c("d","s","r","hidx")
  burn_in=1000;
  steps=10000;
  thin=2;
  chains = 3
  
  
  ##Larva Event Counts Slice
  datSliceF <- datHuntVsPrey[datHuntVsPrey[,1] >= preyCntRange[1] & datHuntVsPrey[,1] <= preyCntRange[2], ]
  
  datJags=list(d=datHuntVsPrey[,3],NTOT=NROW(datHuntVsPrey));
  
  nDat = NROW(datHuntVsPrey);
  dataG=list(d=datHuntVsPrey$DurationFrames,NTOT=nDat,hidx=datHuntVsPrey$hidx,food=as.integer(datSliceF[,1]));
  
  model=jags.model(file=strModelFilename,data=dataG,n.chains=chains,inits=initEventDurfunct(chains,NROW(unique(datHuntVsPrey$hidx)) ) );
  
  update(model,burn_in)
  
  drawSamples=jags.samples(model,steps,thin=thin,variable.names=varnames1)
  
  return(drawSamples) 
}

## Plots A collection of Gamma Distributions given the vector of parameters Shape and Rate
plotDurationGammaSamples <- function (gammaShape,gammaRate,lcolour)
{
  plot(1:2000/G_APPROXFPS,dgamma(1:2000,shape=gammaShape[1],rate=gammaRate[1]),ylim=c(0,0.05),type="l",col=lcolour)
  for (i in 1:NROW(gammaShape))
    lines(1:2000/G_APPROXFPS,dgamma(1:2000,shape=gammaShape[i],rate=gammaRate[i]),col=lcolour)
}

GammaShapeSamples <- function(drawHuntDuration,plotsamples)
{
  vGammaShape <- vector()#tail(drawHD_NE$r[1,,schain],plotsamples)  
  for (i in 1:NROW(drawHuntDuration$r))
    vGammaShape <- append(vGammaShape,tail(drawHuntDuration$s[i,,schain],plotsamples))
  
  return(vGammaShape)
}
GammaScaleSamples <- function(drawHuntDuration,plotsamples)
{
  vGammaScale <- vector()#1/tail(drawHD_NE$r[1,,schain],plotsamples)
  for (i in 1:NROW(drawHuntDuration$r))
    vGammaScale <- append(vGammaScale,1/tail(drawHuntDuration$r[i,,schain],plotsamples))
  
  return(vGammaScale)
}


plotDurationGammaSamples <- function (gammaShape,gammaRate,lcolour)
{
  plot(1:3000/G_APPROXFPS,dgamma(1:3000,shape=gammaShape[1],rate=gammaRate[1]),ylim=c(0,0.05),type="l",col=lcolour)
  for (i in 1:NROW(gammaShape))
    lines(1:3000/G_APPROXFPS,dgamma(1:3000,shape=gammaShape[i],rate=gammaRate[i]),col=lcolour)
}


############ LOAD DATA #################

############ LOAD EVENTS LIst and Fix ####
## Warning Set Includes Repeated Test For some LF fish - One In Different Food Density
## Merged2 Contains the Fixed, Remerged EventID 0 files, so event Counts appear for all larvae recorded.
#strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged3" 


## Load From Central Function
## Removeshese Are Double/2nd Trials on LL, or Simply LL unpaired to any LE (Was checking Rates)
datHuntLabelledEventsSBMerged_fixed <- getLabelledHuntEventsSet()


################# # ## # # 
## Get Summarized Hunt Results Per Larva ####
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged_fixed)

## Get Event Duration Per Group ###
datHEvent_LE <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"LE")
datHEvent_NE <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"NE")
datHEvent_DE <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"DE")
datHEvent_LL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"LL")
datHEvent_NL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"NL")
datHEvent_DL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_fixed,"DL")

### Load Prior RJags Sampled Values ###
load(file =paste(strDataExportDir,"stat_HuntDurationInPreyRange_nbinomRJags.RData",sep="")) ## Total Hint Duration per Larvae
load(file =paste(strDataExportDir,"stat_HEventDurationInPreyRange_nbinomRJags.RData",sep="") ) ## Hunt Episode Duration


## Baysian Inference Fiting a Gamma distribution to the Hunt Event Duration Data ##
##Setup Data Structure To Pass To RJAgs

## Get Event Counts Within Range  - Along With Total Number of Hunting frames for each Larva##
## Added Larva ID to Check for Correlation Through Time of Day - Surrogate as LarvaID;s increase through the day of the experiment from 1 morning -4 evening
## No Effect Of Time Of Data is evident Found #

##***## NOTE Remove Rec Without Any Hunt Events / ie 0 Duration ### 
###  This is because The Hunt Duration Is Discontinuous if we include the Larvae With 0 events, and this creates modelling problems 
## Thus when Looking for Hunt Duration per larva we need to exclude the ones that did not hunt.
datHuntDurVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LL ),datHuntStat[,"vIDLookupTable"]$LL$larvaID )
datHuntDurVsPreyLL <- datHuntDurVsPreyLL[!is.na(datHuntDurVsPreyLL[,1]) & datHuntDurVsPreyLL[,2] > 0 ,]
datHuntDurVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LE ),datHuntStat[,"vIDLookupTable"]$LE$larvaID  )
datHuntDurVsPreyLE <- datHuntDurVsPreyLE[!is.na(datHuntDurVsPreyLE[,1]) & datHuntDurVsPreyLE[,2] > 0,] ##& datHuntDurVsPreyLE[,2] > 0

datHuntDurVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NL),datHuntStat[,"vIDLookupTable"]$NL$larvaID )
datHuntDurVsPreyNL <- datHuntDurVsPreyNL[!is.na(datHuntDurVsPreyNL[,1]) & datHuntDurVsPreyNL[,2] > 0,]
datHuntDurVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NE),datHuntStat[,"vIDLookupTable"]$NE$larvaID  )
datHuntDurVsPreyNE <- datHuntDurVsPreyNE[!is.na(datHuntDurVsPreyNE[,1]) & datHuntDurVsPreyNE[,2] > 0  ,]

datHuntDurVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DL ),datHuntStat[,"vIDLookupTable"]$DL$larvaID  )
datHuntDurVsPreyDL <- datHuntDurVsPreyDL[!is.na(datHuntDurVsPreyDL[,1]) & datHuntDurVsPreyDL[,2] > 0,] ##Remove datHuntDurVsPreyDL[,2] NA And High Fliers
datHuntDurVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DE ),datHuntStat[,"vIDLookupTable"]$DE$larvaID  )
datHuntDurVsPreyDE <- datHuntDurVsPreyDE[!is.na(datHuntDurVsPreyDE[,1]) & datHuntDurVsPreyDE[,2] > 0,] ##Remove NA And High Fliers


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)
plotsamples <- 200
schain <-1:3

## Run Sampler On Larva Total Hunt Duration 
drawDurLE <- mcmc_drawHuntDurationModels(datHuntDurVsPreyLE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurNE <- mcmc_drawHuntDurationModels(datHuntDurVsPreyNE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurDE <- mcmc_drawHuntDurationModels(datHuntDurVsPreyDE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurLL <- mcmc_drawHuntDurationModels(datHuntDurVsPreyLL,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurDL <- mcmc_drawHuntDurationModels(datHuntDurVsPreyDL,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurNL <- mcmc_drawHuntDurationModels(datHuntDurVsPreyNL,preyCntRange,"modelLarvaHuntDuration.tmp" )
##
save(drawDurLE,drawDurNE,drawDurDE,drawDurLL,drawDurNL,drawDurDL,
     file =paste(strDataExportDir,"stat_HuntDurationInPreyRange_nbinomRJags.RData",sep=""))


### Run Sampler On Hunt Event Duration , Modelling Duration distribution Per Larva
drawHD_LE <-mcmc_drawEventDurationModels (datHEvent_LE,preyCntRange, "modelLarvaEventDuration.tmp")
drawHD_NE <-mcmc_drawEventDurationModels (datHEvent_NE,preyCntRange, "modelLarvaEventDuration.tmp")
drawHD_DE <-mcmc_drawEventDurationModels (datHEvent_DE,preyCntRange, "modelLarvaEventDuration.tmp")
drawHD_NL <-mcmc_drawEventDurationModels (datHEvent_NL,preyCntRange, "modelLarvaEventDuration.tmp")
drawHD_LL <-mcmc_drawEventDurationModels (datHEvent_LL,preyCntRange, "modelLarvaEventDuration.tmp")
drawHD_DL <-mcmc_drawEventDurationModels (datHEvent_DL,preyCntRange, "modelLarvaEventDuration.tmp")
##
save(drawHD_LE,drawHD_NE,drawHD_DE,drawHD_NL,drawHD_LL,drawHD_DL,
     file =paste(strDataExportDir,"stat_HEventDurationInPreyRange_nbinomRJags.RData",sep=""))


## Make Derived Duration Results ##
### The Prob Of Success p from NegBinom translates to Gamma Rate p/(1-p), or scale: (1-p)/p
## This Rate is the Hunt Frames Producing rate now ###
HLarvaHuntGammaRate_LE <-tail(drawDurLE$q[,,schain],plotsamples)/(1-tail(drawDurLE$q[,,schain],plotsamples));
HLarvaHuntGammaRate_LL <-tail(drawDurLL$q[,,schain],plotsamples)/(1-tail(drawDurLL$q[,,schain],plotsamples));
HLarvaHuntGammaRate_DE <- tail(drawDurDE$q[,,schain],plotsamples)/(1-tail(drawDurDE$q[,,schain],plotsamples));
HLarvaHuntGammaRate_DL <- (tail(drawDurDL$q[,,schain],plotsamples)/(1-tail(drawDurDL$q[,,schain],plotsamples)));      
HLarvaHuntGammaRate_NE <- (tail(drawDurNE$q[,,schain],plotsamples)/(1-tail(drawDurNE$q[,,schain],plotsamples)));
HLarvaHuntGammaRate_NL <- (tail(drawDurNL$q[,,schain],plotsamples)/(1-tail(drawDurNL$q[,,schain],plotsamples)));      
HLarvaHuntGammaShape_LE <- tail(drawDurLE$r[,,schain],plotsamples);
HLarvaHuntGammaShape_LL <- tail(drawDurLL$r[,,schain],plotsamples)
HLarvaHuntGammaShape_DE <- tail(drawDurDE$r[,,schain],plotsamples);
HLarvaHuntGammaShape_DL <- tail(drawDurDL$r[,,schain],plotsamples)
HLarvaHuntGammaShape_NE <- tail(drawDurNE$r[,,schain],plotsamples);
HLarvaHuntGammaShape_NL <- tail(drawDurNL$r[,,schain],plotsamples)


# These Functions Obtain nS Samples from each of the Posteriors of the Modelled hunt events of separate Larva, 
nS <- 300
muEpiDur_NE <- GammaShapeSamples(drawHD_NE,nS)*GammaScaleSamples(drawHD_NE,nS)/G_APPROXFPS
muEpiDur_DE <- GammaShapeSamples(drawHD_DE,nS)*GammaScaleSamples(drawHD_DE,nS)/G_APPROXFPS
muEpiDur_LE <- GammaShapeSamples(drawHD_LE,nS)*GammaScaleSamples(drawHD_LE,nS)/G_APPROXFPS
muEpiDur_NL <- GammaShapeSamples(drawHD_NL,nS)*GammaScaleSamples(drawHD_NL,nS)/G_APPROXFPS
muEpiDur_DL <- GammaShapeSamples(drawHD_DL,nS)*GammaScaleSamples(drawHD_DL,nS)/G_APPROXFPS
muEpiDur_LL <- GammaShapeSamples(drawHD_LL,nS)*GammaScaleSamples(drawHD_LL,nS)/G_APPROXFPS

pBW <- 0.25
densEpiDur_NE <- density(muEpiDur_NE,bw=pBW)
densEpiDur_LE <- density(muEpiDur_LE,bw=pBW)
densEpiDur_DE <- density(muEpiDur_DE,bw=pBW)
densEpiDur_NL <- density(muEpiDur_NL,bw=pBW)
densEpiDur_LL <- density(muEpiDur_LL,bw=pBW)
densEpiDur_DL <- density(muEpiDur_DL,bw=pBW)
## MAIN PLOT ###
#### HUNT EVENT PER LARVA PLOT #####
## Comprehensive Plot On Number of Hunt Events
pdf(file= paste(strPlotExportPath,"/stat/fig2.B_statComparePoissonHuntDurations",".pdf",sep=""),width = 14,height = 7)
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 3.2
cex = 1.4
adj  = 2.5
padj <- -7.5
las <- 1

layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2,6, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,1,1))

plotHuntDurationDistribution_cdf(datHuntDurVsPreyNE,drawDurNE,colourHE[1],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntDurVsPreyNL,drawDurNL,colourHL[1],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)
legend("bottomright",legend = c(  expression (),
                                  bquote(NF["s"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyNE))  )
                                  ,bquote(NF["s"] ~"Model " ), 
                                  bquote(NF["e"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyNL)) )
                                  ,bquote( NF["e"] ~ "Model " ) ), 
       col=c(colourP[4], colourLegE[1],colourP[4],colourLegL[1]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=cex,bg="white" )
mtext("F",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
mtext(side = 1,cex=cex, line = lineAxis, " Total time spent hunting (sec)")
mtext(side = 2,cex=cex, line = lineAxis, " Cumulative function ")

plotHuntDurationDistribution_cdf(datHuntDurVsPreyLE,drawDurLE,colourHE[2],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntDurVsPreyLL,drawDurLL,colourHL[2],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)
legend("bottomright",legend = c(  expression (),
                                  bquote(LF["s"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyLE))  )
                                  ,bquote(LF["s"] ~"Model " ), 
                                  bquote(LF["e"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyLL)) )
                                  ,bquote( LF["e"] ~ "Model " ) ),  
       col=c(colourP[4], colourLegE[2],colourP[4],colourLegL[2]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=cex,bg="white" )
mtext("G",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
mtext(side = 1,cex=cex, line = lineAxis, " Total time spent hunting (sec)")
mtext(side = 2,cex=cex, line = lineAxis, " Cumulative function ")

plotHuntDurationDistribution_cdf(datHuntDurVsPreyDE,drawDurDE,colourHE[3],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntDurVsPreyDL,drawDurDL,colourHL[3],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)
legend("bottomright",legend = c(  expression (),
                                  bquote(DF["s"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyDE))  )
                                  ,bquote(DF["s"] ~"Model " ), 
                                  bquote(DF["e"] ~ 'Data #' ~ .(NROW(datHuntDurVsPreyDL)) )
                                  ,bquote( DF["e"] ~ "Model " ) ),  
       col=c(colourP[4], colourLegE[3],colourP[4],colourLegL[3]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=cex,bg="white" )
mtext("H",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
mtext(side = 1,cex=cex, line = lineAxis, " Total time spent hunting (sec)")
mtext(side = 2,cex=cex, line = lineAxis, " Cumulative function ")

######
###### Gamma Parameters Comparison###

pchL <- c(1,2,0,16,17,15)
### Plot GAMMA Parameters Space
#Xlim <- range(1/HLarvaHuntGammaRate_DL/G_APPROXFPS)[2]
#plot(1/HLarvaHuntGammaRate_NE/G_APPROXFPS,HLarvaHuntGammaShape_NE,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[1],xlab=NA,ylab=NA)
#points(1/HLarvaHuntGammaRate_LE/G_APPROXFPS,HLarvaHuntGammaShape_LE,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[2])
#points(1/HLarvaHuntGammaRate_DE/G_APPROXFPS,HLarvaHuntGammaShape_DE,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[3])
#points(1/HLarvaHuntGammaRate_NL/G_APPROXFPS,HLarvaHuntGammaShape_NL,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[4])
#points(1/HLarvaHuntGammaRate_LL/G_APPROXFPS,HLarvaHuntGammaShape_LL,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[5])
#points(1/HLarvaHuntGammaRate_DL/G_APPROXFPS,HLarvaHuntGammaShape_DL,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[6])
#strXLab <- (expression(paste(Gamma, " scale (r/FPS)") ) )
#mtext(side = 1,cex=0.8, line = 2.2,strXLab ) 
#mtext(side = 2,cex=0.8, line = 2.2, expression(paste(Gamma, " shape (k)") ) )
#legend("topright", legend = strDataLabels, #c(paste("NE" ), paste("LE"),paste("DE"), paste("NL"),paste("LL"),paste("DL")),
#       col=c(colourHL[1],colourHL[2],colourHL[3],colourHL[1],colourHL[2],colourHL[3]) ,pch=pchL,cex=0.9,bg="white",ncol=2)
#mtext("D",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

###

## BoxPlot of Hunt Event Counts - 
strCondTags <- c("NE","NL","LE","LL","DE","DL")
xbarcenters <- boxplot(log10( (datHuntDurVsPreyNE[,3]+1)/G_APPROXFPS ) ,log10( ( datHuntDurVsPreyNL[,3]+1)/G_APPROXFPS ),log10( (datHuntDurVsPreyLE[,3]+1)/G_APPROXFPS ),
                       log10( ( datHuntDurVsPreyLL[,3]+1 )/G_APPROXFPS ),log10(( datHuntDurVsPreyDE[,3]+1)/G_APPROXFPS ) ,log10( ( datHuntDurVsPreyDL[,3]+1)/G_APPROXFPS ),
                       main=NA,notch=TRUE,col=colourD,names=strCondTags,ylim=c(0,2.5),axes = FALSE,cex=cex,cex.axis=cex,cex.lab=cex  )
mtext(side = 2,cex=cex, line =lineAxis, "Total hunting time  (sec) ") #log( (D+1)/fps
vIDTable    <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]
##Take Frame duration Values for each group and Divide by Framerate
vDat        <- rapply(datHuntStat[,"vHDurationPerLarva"],function(x){return (x/G_APPROXFPS) },how="replace")

axis(1,at<-axis(1,labels=NA),cex=cex,cex.axis=cex,cex.lab=cex, labels=c( strDataLabels[1],strDataLabels[4],strDataLabels[2],strDataLabels[5],strDataLabels[3],strDataLabels[6] ))
yticks <-axis(2,labels=NA)
axis(2, at = yticks, labels =round(10^yticks) , col.axis="black", las=2)
## Connect Larvae From EMpty To LIve Test Condition #
plotConnectedHuntDuration(datHuntStat,vDat,strCondTags)
mtext("I",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj+3,cex.main=cex,cex=cex)


#### Show Density Of Hunt Episode Duration per Hunt Event ####
## PLot Expected Duration - as the Gamma Mean 
## RAW Hunt Duration DATA 
plot(density(datHEvent_NE$DurationFrames/G_APPROXFPS ),type='l',xlim=c(0,6),ylim=c(0,1),lty=lineTypeL[1],col=colourLegL[1],lwd=4,ylab=NA,xlab=NA,main=NA,cex=cex,cex.axis=cex,cex.lab=cex)
lines(density(datHEvent_LE$DurationFrames/G_APPROXFPS ),xlim=c(0,6),col=colourLegL[2],lty=lineTypeL[1],lwd=4,ylab=NA,xlab=NA)
lines(density(datHEvent_DE$DurationFrames/G_APPROXFPS ),xlim=c(0,6),col=colourLegL[3],lty=lineTypeL[1],lwd=4,ylab=NA,xlab=NA)

lines(density(datHEvent_LL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[2],lwd=4,ylab=NA,xlab=NA)
lines(density(datHEvent_DL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[3],lwd=4,ylab=NA,xlab=NA)
lines(density(datHEvent_NL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[1],lwd=4,ylab=NA,xlab=NA)
legend("topright",legend = c(paste("Spontaneous " ),paste("Evoked ")), seg.len=3.5,
       col=c(colourR[4], colourR[4]),lty=c(2,1),lwd=4,cex=1.1,bg="white" )
mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Measured duration of each hunt episode  (sec)") )  )

mtext(side = 2,cex=cex, line = lineAxis, " Density function ")
mtext("J",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)


dev.off()

#### END OF MAIN PLOT ###
#########################


######### Gamma Episode Duration plot ###
pdf(file= paste(strPlotExportPath,"/stat/fig2S_statModelGammaHuntEpisodeDurations.pdf",sep=""),width = 14,height = 7)

layout(matrix(c(1,1,2,2), 1,4, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,1,1))

  plot(densEpiDur_NE,type='l',xlim=c(0,6),ylim=c(0,1),lty=lineTypeL[1],col=colourLegL[1],lwd=4,ylab=NA,xlab=NA,main=NA,cex=cex,cex.axis=cex,cex.lab=cex)
  lines(density(datHEvent_NE$DurationFrames/G_APPROXFPS ),xlim=c(0,6),col=colourLegL[1],lty=lineTypeL[1],lwd=2,ylab=NA,xlab=NA)
  lines(densEpiDur_LE,xlim=c(0,6),col=colourLegL[2],lty=lineTypeL[1],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHEvent_LE$DurationFrames/G_APPROXFPS ),xlim=c(0,6),col=colourLegL[2],lty=lineTypeL[1],lwd=2,ylab=NA,xlab=NA)
  lines(densEpiDur_DE,xlim=c(0,6),col=colourLegL[3],lty=lineTypeL[1],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHEvent_DE$DurationFrames/G_APPROXFPS ),xlim=c(0,6),col=colourLegL[3],lty=lineTypeL[1],lwd=2,ylab=NA,xlab=NA)
  
  
  plot(densEpiDur_LL,xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[2],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHEvent_LL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[2],lwd=2,ylab=NA,xlab=NA)
  lines(densEpiDur_DL,xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[3],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHEvent_DL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[3],lwd=2,ylab=NA,xlab=NA)
  lines(densEpiDur_NL,xlim=c(0,8),lty=lineTypeL[2],col=colourLegL[1],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHEvent_NL$DurationFrames/G_APPROXFPS ),xlim=c(0,6),lty=lineTypeL[2],col=colourLegL[1],lwd=2,ylab=NA,xlab=NA)
  
  legend("topright",legend = c(paste("Spontaneous " ),paste("Evoked ")), seg.len=3.5,
         col=c(colourR[4], colourR[4]),lty=c(2,1),lwd=4,cex=1.1,bg="white" )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Estimated duration of each hunt episode  (sec)") )  )
  mtext(side = 2,cex=cex, line = lineAxis, " Density function ")
  
#  mtext("J",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex)
dev.off()  


###Courtesy of https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
getProbOfInterval  <- function(densEpiDur,xmin,xmax)
{
  #### Calculate Percentiles and Compare Results
  xx <- densEpiDur$x  ## 512 evenly spaced points on [min(x) - 3 * d$bw, max(x) + 3 * d$bw]
  yy <- densEpiDur$y ## 512 density values for `xx`
  dx <- xx[2L] - xx[1L]  ## spacing / bin size
  
  afunct <- approxfun(xx, yy)
  Const <- integrate(afunct, min(xx), max(xx))$value 
  p.unscaled <- integrate(afunct, xmin, min(xmax,max(xx) ) ) $value
  p.scaled <- p.unscaled / Const
  
  return(p.scaled)
}
### DURATION DATA Comparison
##Check the duration that contains 95% of LE Values 
pLEMx <- quantile(datHEvent_LE$DurationFrames/G_APPROXFPS,probs=0.95)
##Compare Percentage to Other Groups
getProbOfInterval(density(datHEvent_DE$DurationFrames/G_APPROXFPS),0,pLEMx) ##73% of huntdurations in DE
getProbOfInterval(density(datHEvent_NE$DurationFrames/G_APPROXFPS),0,pLEMx) ## 69% in NE


##LIVE FED DATA
pLLMx <- quantile(datHEvent_LL$DurationFrames/G_APPROXFPS,probs=0.95)
##Compare Percentage to Other Groups
getProbOfInterval(density(datHEvent_DL$DurationFrames/G_APPROXFPS),0,pLLMx) ##73% of huntdurations in DE
getProbOfInterval(density(datHEvent_NL$DurationFrames/G_APPROXFPS),0,pLLMx) ## 69% in NE

##Modulation By Prey at the 95th percentile
getProbOfInterval(density(datHEvent_LL$DurationFrames/G_APPROXFPS),0,quantile(datHEvent_LE$DurationFrames/G_APPROXFPS,probs=0.80)) ## LL On 95% point of LE
getProbOfInterval(density(datHEvent_NL$DurationFrames/G_APPROXFPS),0,quantile(datHEvent_NE$DurationFrames/G_APPROXFPS,probs=0.80)) ## LL On 95% point of LE
getProbOfInterval(density(datHEvent_DL$DurationFrames/G_APPROXFPS),0,quantile(datHEvent_DE$DurationFrames/G_APPROXFPS,probs=0.80)) ## LL On 95% point of LE

###Means 
errLE<- sd(datHEvent_LE$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_LE$DurationFrames) )
errDE <- sd(datHEvent_DE$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_DE$DurationFrames) )
errNE <- sd(datHEvent_NE$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_NE$DurationFrames) )

message("LE D mean:",mean(datHEvent_LE$DurationFrames/G_APPROXFPS),"+SEM:",errLE,"\n",
        "DE D mean:",mean(datHEvent_DE$DurationFrames/G_APPROXFPS),"+SEM:",errDE,"\n",
        "NE D mean:",mean(datHEvent_NE$DurationFrames/G_APPROXFPS),"+SEM:",errNE)



errLL<- sd(datHEvent_LL$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_LL$DurationFrames) )
errDL <- sd(datHEvent_DL$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_DL$DurationFrames) )
errNL <- sd(datHEvent_NL$DurationFrames/G_APPROXFPS) /sqrt(NROW(datHEvent_NL$DurationFrames) )

message("LL D mean:",mean(datHEvent_LL$DurationFrames/G_APPROXFPS),"+SEM:",errLL,"\n",
      "DL D mean:",mean(datHEvent_DL$DurationFrames/G_APPROXFPS),"+SEM:",errDL,"\n",
      "NL D mean:",mean(datHEvent_NL$DurationFrames/G_APPROXFPS),"+SEM:",errNL)


message("LL D median:",median(datHEvent_LL$DurationFrames/G_APPROXFPS),"+SEM:",errLL,"\n",
        "DL D median:",median(datHEvent_DL$DurationFrames/G_APPROXFPS),"+SEM:",errDL,"\n",
        "NL D median:",median(datHEvent_NL$DurationFrames/G_APPROXFPS),"+SEM:",errNL)


quantile(datHEvent_LE$DurationFrames/G_APPROXFPS,p=0.95)
quantile(datHEvent_LL$DurationFrames/G_APPROXFPS,p=0.95)

quantile(datHEvent_NE$DurationFrames/G_APPROXFPS,p=0.95)
quantile(datHEvent_NL$DurationFrames/G_APPROXFPS,p=0.95)


quantile(datHEvent_DE$DurationFrames/G_APPROXFPS,p=0.95)
quantile(datHEvent_DL$DurationFrames/G_APPROXFPS,p=0.95)


## MODEL ###
##Check the duration that contains 95% of LE Values 
pLEMx <- quantile(muEpiDur_LE,probs=0.95)
##Compare Percentage to Other Groups
getProbOfInterval(densEpiDur_DE,0,pLEMx) ##73% of huntdurations in DE
getProbOfInterval(densEpiDur_NE,0,pLEMx) ## 69% in NE

## Now in Livefed conditions
pLLMx <- quantile(muEpiDur_LL,probs=0.95)
##Compare Percentage to Other Groups
getProbOfInterval(densEpiDur_DL,0,pLLMx) ##73% of huntdurations in DE
getProbOfInterval(densEpiDur_NL,0,pLLMx) ## 69% in NE

quantile(muEpiDur_LE,probs=0.95)
getProbOfInterval(densEpiDur_LE,0,pLEMx) ##73% of huntdurations in DE
quantile(muEpiDur_LL,probs=0.95)

quantile(muEpiDur_NE,probs=0.95)
quantile(muEpiDur_NL,probs=0.95)

quantile(muEpiDur_DE,probs=0.80)
quantile(muEpiDur_DL,probs=0.80)

### Evoked - Spontaneous Overlap
getProbOfInterval(densEpiDur_LL,0,pLEMx) - 0.95 ##73% of huntdurations in DE
getProbOfInterval(densEpiDur_NL,0,quantile(muEpiDur_NE,probs=0.95)) - 0.95 ## Ovelap Between Evoked - Spontaneous NL
getProbOfInterval(densEpiDur_DL,0,quantile(muEpiDur_DE,probs=0.95)) - 0.95## 
























##############Check Successfull Only ###################
##Remove Dublicates - Choose Labels - Duration Needs To be > 5ms
datHuntLabelledEventsSBMerged_Success_filtered <- datHuntLabelledEventsSBMerged_filtered  [
  with(datHuntLabelledEventsSBMerged_filtered , ( ( convertToScoreLabel(huntScore) == "Success-OnApproach" |
                                          convertToScoreLabel(huntScore) == "Success-OnApproach-AfterStrike" |
                                          convertToScoreLabel(huntScore) == "Success" |
                                            convertToScoreLabel(huntScore) == "Success-OnStrike" |
                                            convertToScoreLabel(huntScore) == "Success-SpitBackOut" |
                                            convertToScoreLabel(huntScore) == "Success-OnStrike-SpitBackOut") &
                                          ( (endFrame - startFrame) > 200 ) |  ## limit min event dur to 5ms
                                            eventID == 0)), ] ## Add the 0 Event, In Case Larva Produced No Events

##Remove Dublicates - Choose Labels - Duration Needs To be > 5ms
datHuntLabelledEventsSBMerged_Failed_filtered <- datHuntLabelledEventsSBMerged_filtered  [
  with(datHuntLabelledEventsSBMerged_filtered , ( ( convertToScoreLabel(huntScore) == "Fail-No Strike" |
                                                      convertToScoreLabel(huntScore) == "Fail-With Strike" |
                                                      convertToScoreLabel(huntScore) == "Fail") &
                                                    ( (endFrame - startFrame) > 200 ) |  ## limit min event dur to 5ms
                                                    eventID == 0)), ] ## Add the 0 Event, In Case Larva Produced No Events



datHuntSuccDur_LL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Success_filtered,"LL")
datHuntSuccDur_DL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Success_filtered,"DL")
datHuntSuccDur_NL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Success_filtered,"NL")

datHuntFailDur_LL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Failed_filtered,"LL")
datHuntFailDur_DL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Failed_filtered,"DL")
datHuntFailDur_NL <- getHuntEventDuration(datHuntLabelledEventsSBMerged_Failed_filtered,"NL")


datHuntSuccessStat <- makeHuntStat(datHuntLabelledEventsSBMerged_Success_filtered)

#### Compare How Success Episode duratins shift / Appears differences are beyond the 50% 
quantile(datHuntSuccDur_NL$DurationFrames/G_APPROXFPS)
quantile(datHuntSuccDur_LL$DurationFrames/G_APPROXFPS)
quantile(datHuntSuccDur_DL$DurationFrames/G_APPROXFPS)


#### Compare How Success Episode duratins shift / Appears differences are beyond the 50% 
quantile(datHuntFailDur_NL$DurationFrames/G_APPROXFPS)
quantile(datHuntFailDur_LL$DurationFrames/G_APPROXFPS)
quantile(datHuntFailDur_DL$DurationFrames/G_APPROXFPS)


#     ####PLOT SUCC DENSITY
pdf(file= paste(strPlotExportPath,"/stat/fig2.B_statComparePoissonHuntDurations",".pdf",sep=""),width = 14,height = 7)
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 3.2
cex = 1.4
adj  = 2.5
padj <- -7.5
las <- 1

  plot(density(datHuntSuccDur_NL$DurationFrames/G_APPROXFPS),type='l',xlim=c(0,8),ylim=c(0,1),lty=lineTypeL[3],col=colourL[1],lwd=4,ylab=NA,xlab=NA,main=NA,cex=cex,cex.axis=cex,cex.lab=cex)
  lines(density(datHuntSuccDur_LL$DurationFrames/G_APPROXFPS),xlim=c(0,8),col=colourL[2],lty=lineTypeL[3],lwd=4,ylab=NA,xlab=NA)
  lines(density(datHuntSuccDur_DL$DurationFrames/G_APPROXFPS),xlim=c(0,8),col=colourL[3],lty=lineTypeL[3],lwd=4,ylab=NA,xlab=NA)
  
  legend("topright",legend = c(paste("Evoked Successful ")), seg.len=3.5,
         col=c(colourR[4], colourR[4]),lty=c(3,1),lwd=4,cex=1.1,bg="white" )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Measured duration of each hunt episode  (sec)") )  )
  mtext(side = 2,cex=cex, line = lineAxis, " Density function ")

dev.off()

# 
# 
# ## Continued...
# 
# 
# 
# 
# 
# 
# #### Further Plot Checks ####
# 
# ###########
# ##CHECK fIT
# plot(1:15000,dnbinom(1:15000, size=tail(drawDurLE$r,50),  prob=tail(drawDurLE$q,50)),cex=0.3,col=colourHE[2] )
# lines(density(drawDurLE$d,bw=1000 ) )
# plot(1:15000,dnbinom(1:15000, size=tail(drawDurNE$r,50),  prob=tail(drawDurNE$q,50)),cex=0.3,col=colourHE[1] )
# lines(density(drawDurNE$d,bw=1000 ) )
# plot(1:15000,dnbinom(1:15000, size=tail(drawDurDE$r,50),  prob=tail(drawDurDE$q,50)),cex=0.3,col=colourHE[3] )
# lines(density(drawDurDE$d,bw=1000 ) )
# 
# plot(1:15000,dnbinom(1:15000, size=tail(drawDurNL$r,50),  prob=tail(drawDurNL$q,50)),cex=0.3,col=colourHL[1] )
# lines(density(drawDurNL$d,bw=1000 ) )
# plot(1:15000,dnbinom(1:15000, size=tail(drawDurLL$r,50),  prob=tail(drawDurLL$q,50)),cex=0.3,col=colourHL[2] )
# lines(density(drawDurLL$d,bw=1000 ) )
# 
# #####################################################   ############### ######################################
# 
# ############### Model Each Larva s Event Duration Separatelly ################################################
# ###            Use The Labelled Events Register               #### ##################
# ## Run Baysian Inference on Model for Hunt Event Counts IN A Group/ Test Condition
# ## Return Samples Drawn structure                             ########
# 
# hist(datHEvent_LE$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# hist(datHEvent_LL$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# 
# hist(datHEvent_NE$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# hist(datHEvent_NL$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# 
# hist(datHEvent_DE$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# hist(datHEvent_DL$DurationFrames/G_APPROXFPS,xlim=c(0,10),breaks=30)
# 
# 
# ##
# ### Cut And Examine The data Where There Are Between L and M rotifers Initially
# preyCntRange <- c(0,100)
# 
# ##Retrieve Draws From Each Larva's Episodes
# ## Draw Structure Includes Infered Gamma Distributions to Fit Duration of Hunt Episode of Each Larva.
# plotDurationGammaSamples(drawHD_LE$s,drawHD_LE$r,"green")
# plotDurationGammaSamples(drawHD_DE$s,drawHD_DE$r,"blue")
# plotDurationGammaSamples(drawHD_NE$s,drawHD_NE$r,"red")
# 
# plotDurationGammaSamples(drawHD_LL$s,drawHD_LL$r,"green")
# plotDurationGammaSamples(drawHD_NL$s,drawHD_NL$r,"red")
# plotDurationGammaSamples(drawHD_DL$s,drawHD_DL$r,"blue")
# 
# 
# ###############
# ### Plot GAMMA Parameters Space
# Xlim <- 5 #range(1/drawHD_LL$r)[2]
# Ylim <- range(drawHD_LL$s)[2]
# nS <- 2
# plot(GammaScaleSamples(drawHD_NE,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_NE,nS),col=colourHL[1],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[1],xlab=NA,ylab=NA)
# points(GammaScaleSamples(drawHD_LE,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_LE,nS),col=colourHL[2],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[2])
# points(GammaScaleSamples(drawHD_DE,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_DE,nS),col=colourHL[3],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[3])
# points(GammaScaleSamples(drawHD_NL,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_NL,nS),col=colourHL[1],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[4])
# points(GammaScaleSamples(drawHD_LL,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_LL,nS),col=colourHL[2],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[5])
# points(GammaScaleSamples(drawHD_DL,nS)/G_APPROXFPS,GammaShapeSamples(drawHD_DL,nS),col=colourHL[3],ylim=c(0,Ylim),xlim=c(0,Xlim),pch=pchL[6])
# strXLab <- (expression(paste(Gamma, " scale (r/FPS)") ) )
# mtext(side = 1,cex=0.8, line = 2.2,strXLab ) 
# mtext(side = 2,cex=0.8, line = 2.2, expression(paste(Gamma, " shape (k)") ) )
# legend("topright",legend = c(paste("NE" ),
#                              paste("LE"),paste("DE"), paste("NL"),paste("LL"),paste("DL")),
#        col=c(colourHL[1],colourHL[2],colourHL[3],colourHL[1],colourHL[2],colourHL[3]) ,pch=pchL,cex=0.9,bg="white",ncol=2)

# 
# ######## Plot Comparison Of Duration Data ##########
# ##Now Plot Infered Distributions
# 
# 
# ## Plot Histogram Of Durations in approx SEC
# pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
# 
# layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
# hist(datHDuration_LE$DurationFrames /G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[2],xlab="",ylim=c(0,60),main="LE")
# hist(datHDuration_NE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[3],main="NE",xlab="",ylim=c(0,60))
# hist(datHDuration_DE$DurationFrames/G_APPROXFPS,breaks=seq(0,12,0.5),col=colourR[1],main="DE",xlab="Duration of Spontaneous Hunt Events (sec) ",ylim=c(0,60))
# dev.off()
# 
# ## Raw Histograms for Total Hunt Duration Per Larva
# pdf(file= paste(strPlotExportPath,"/stat/stat_SpontaneousTotalHuntDurationPerLarva",preyCntRange[1],"-",preyCntRange[2], "_hist.pdf",sep=""))
# layout(matrix(c(1,2,3), 3,1, byrow = FALSE))
# hist(datHuntVsPreyL[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[2],main="LE",xlab="",xlim=c(0,40),ylim=c(0,40))
# hist(datHuntVsPreyN[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[3],main="NE",xlab="",xlim=c(0,40),ylim=c(0,40))
# hist(datHuntVsPreyD[,3]/G_APPROXFPS,breaks=seq(0,51,3),col=colourR[1],main="DE",xlab=" Total Duration per Larva spent in Spontaneous Hunt Events ",xlim=c(0,40),ylim=c(0,40))
# dev.off()
# #### Also Plo
# 
# ## Box plot Of Total Duration Per Larva
# boxplot(datHuntVsPreyN[,3]/G_APPROXFPS,datHuntVsPreyL[,3]/G_APPROXFPS,datHuntVsPreyD[,3]/G_APPROXFPS,
#         main=NA,notch=TRUE,names=c("NE","LE","DE"),ylim=c(0,40), ylab="(sec)",col=colourH )
# 
# 
# 
# 
# #### Plot Density ###
# ###Plot Density of Slope
# 
# ### Make Plot Of Histograms and Gamma Fit 
# strPlotName <- paste(strPlotExportPath,"/stat/stat_SpontaneousHuntEventDuration_p",preyCntRange[1],"-",preyCntRange[2], ".pdf",sep="")
# pdf(strPlotName,width=16,height=10,title="Comparing the Duration of spontaneous hunt events " ) 
# 
# 
# Plim <- max( round(range(datHDuration_L[,1])[2] ),round(range(datHDuration_D[,1])[2] ),round(range(datHDuration_N[,1])[2] ))*1.1   
# layout(matrix(c(1,2,3,4,4,5), 3,2, byrow = FALSE))
# plotDurationDensityFitComparison(datHDuration_L,drawDurL,colourH[2],Plim,10) #colourR[2]
# legend("topright",legend=paste(c("Empirical Density ","Baysian Gamma Fit","Histogram ") )
#        ,pch=c(NA,NA,21),lwd=c(2,1,2),lty=c(2,1,NA),col=c("black",colourL[2],colourH[4]) )
# 
# plotDurationDensityFitComparison(datHDuration_D,drawDurD,colourH[1],Plim,10)
# plotDurationDensityFitComparison(datHDuration_N,drawDurN,colourH[3],Plim,10)
# ## Plot Distrib Of Params
# ns <- 200
# plot(drawDurL$r[,(steps/thin-ns):(steps/thin),],drawDurL$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[2],xlim=c(0,0.01),ylim=c(0,6),
#      xlab="Rate Parameter",ylab="Shape Parameter") 
# points(drawDurD$r[,(steps/thin-ns):(steps/thin),],drawDurD$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[1] ) 
# points(drawDurN$r[,(steps/thin-ns):(steps/thin),],drawDurN$s[,(steps/thin-ns):(steps/thin),] ,type="p",pch=16,cex=1.2,col=colourR[3] ) 
# 
# legend("topright",legend=paste(c("DF #","LF #","NF #"),c(nDatDF,nDatLF ,nDatNF ) )
#        ,pch=16,col=colourL)
# 
# boxplot(datHDuration_L$DurationFrames/G_APPROXFPS,datHDuration_D$DurationFrames/G_APPROXFPS,datHDuration_N$DurationFrames/G_APPROXFPS,
#         main="Hunt Duration per Larva ",notch=TRUE,names=c("LE","DE","NE"),ylim=c(0,6), ylab="(sec)" )
# 
# 
# dev.off()
