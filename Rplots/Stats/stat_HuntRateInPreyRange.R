### Makes Bayssian inference on  hunt events rate and duration - within a given rannge of prey counts ##
## Used to plot spontaneous eye vergence events count in the EMPTY Test conditions ##
## Using an Exp distribution for Event Count, as
#  Poisson Distributions would be justified as a model of random occurance of huntevents through the 10 min recording time. Then at the group level the sum of individual poissons would still give
# render a poisson, thus modelling the group rate. 
# ** Yet the  empirical distribution does not look like poisson- but rather EXP like, heavy on the low event counts, and with a long tail. A Poisson For the group would imply an underlying 
# 
### TODO Also Plot Duration Of Eye Vergence Frames ###
#library(Hmisc) ##For Minor Ticks

#library(fitdistrplus) ## Fpr testing Dist. Fit

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")
source("config_lib.R")

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
## Assuming nbinom(r,p) Poisson(L|a,b) Gamma(a,b) then r=a, p=1/(b+1) -> b=(1-p)/p
## Give geometric
modelGEventRateGeom="model { 
q ~ dunif(0.0,1)
r ~ dgamma(1,1)

for(j in 1:NTOT){
  n[j] ~  dnegbin(q,r) ##Model Number Of Hunt Events Per LArvae
  }
}"



library(rjags)
strModelName = "modelGroupEventRate.tmp"
fileConn=file(strModelName)
writeLines(modelGEventRateGeom,fileConn);
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




## Compare Model TO Data Using CDF ##
plotEventCountDistribution_cdf <- function(datHEventCount,drawHEvent,lcolour,lpch,lty,Plim,nplotSamples=100,newPlot = FALSE)
{
  XLim <- 80
  x <- seq(0,XLim,1)

  cdfD_N <- ecdf(datHEventCount[,2])

  plot(cdfD_N,col=lcolour,pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim),ylim=c(0,1),cex=1.5,cex.lab=1.5,add=!newPlot)
  ##Construct CDF of Model by Sampling randomly from Model distribution for exp rate parameter
  for (c in 1:NROW(drawHEvent$q[1,1,])) {
    for (j in (NROW(drawHEvent$q[,,c])-nplotSamples):NROW(drawHEvent$q[,,c]) )
    {
      cdfM <- dnbinom(x,size=drawHEvent$r[,j,c],prob=  drawHEvent$q[,j,c]  )##1-exp(-q*x) ##ecdf(  dexp( x, q  ) )
      lines(x,cumsum(cdfM),col=lcolour,lty=lty) #add=TRUE,
    }
  }
  plot(cdfD_N,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim),ylim=c(0,1),cex=1.5,cex.lab=1.5,add=TRUE)
  #axis(side = 4)
  #mtext(side = 4, line = 2.1, 'Counts')
  ### Draw Distribution oF Hunt Rates - 
  ## For Poisson- gamma mixture we recover Gamma Params from nbinomial as shape r>shape r and scale theta:(1-p)/p
  
  ##For the EXP Mixture with Poisson The Exp Rate was recovered as :
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

## Plots The Spontaneous and Invoked Hunt Rates, Inferred via the gamma given by the -ve binomial
plotGammaHuntRates <- function(x,HEventHuntGammaShape,HEventHuntGammaRate,lcolour,plotsamples,lineType=1,bnewPlot=FALSE)
{
  schain <- 1:3
  if (bnewPlot)
    plot(x,(dgamma(x,shape=HEventHuntGammaShape[1,1],rate=HEventHuntGammaRate[1,1]) ),main="",xlab=NA,ylab=NA,
       xlim=c(0,15),ylim=c(0,0.5),col=lcolour,type="l",lwd=1, lty=lineType)

  for (c in schain)
  {
    for (i in 1:plotsamples)
    {
      lines(x, ( dgamma(x,shape=HEventHuntGammaShape[i,c],rate=HEventHuntGammaRate[i,c]) ) ,col=lcolour[1],lwd=1, lty=lineType)
      #lines(x,dgamma(x,rate=HEventHuntGammaShape[i,c],HEventHuntGammaRate[i,c]),col=lcolour[2],lwd=3,lty=3)
    }
  }
}
  

## Box Plot Assist To Connect Individual Hunt Event Counts Between Empty (Spontaneous) and Live Test Conditions (Evoked) 
## For Each Larva Experiment
plotConnectedEventCounts <- function(datHuntStat,vDat,strCondTags)
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
      
      segments(gIdx,log10(ptSrc+1),gIdx+1,log10(ptDest+1) ,col=colourP[4])
      
      
    }##For Each Experiment
    
  } ## Go Through Pairs Of Conditions ##
}##End Of Function ConnectEventPoints

  
fromchain=1000

#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq
## Load PreCalc Draws From Model ?? ##
load(file =paste(strDataExportDir,"stat_HuntRateInPreyRange_nbinomRJags.RData",sep=""))

############ LOAD EVENTS LIst and Fix ####
## Warning Set Includes Repeated Test For some LF fish - One In Different Food Density
## Merged2 Contains the Fixed, Remerged EventID 0 files, so event Counts appear for all larvae recorded.
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged3"  ##Load the merged Frames
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSBMerged <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))


##Remove Dublicates - Choose Labels - Duration Needs To be > 5ms
datHuntLabelledEventsSBMerged_filtered <- datHuntLabelledEventsSBMerged [
                with(datHuntLabelledEventsSBMerged, ( convertToScoreLabel(huntScore) != "Not_HuntMode/Delete" &
                                                                 convertToScoreLabel(huntScore) != "Duplicate/Overlapping" &
                                                                  (endFrame - startFrame) > 200 ) |  ## limit min event dur to 5ms
                                                                       eventID == 0), ] ## Add the 0 Event, In Case Larva Produced No Events

##These Are Double/2nd Trials on LL, or Simply LL unpaired to any LE (Was checking Rates)
#AutoSet420fps_14-12-17_WTNotFed2RotiR_297_003.mp4
vxCludeExpID <- c(4421,4611,4541,4351,4481,4501,4411)
vWeirdDataSetID <- c(11,17,18,19) ##These Dataset Have a total N  Exp Less than 4*2*3=24
##Check For Missing Exp Less than 24 

datHuntLabelledEventsSBMerged_fixed <- datHuntLabelledEventsSBMerged_filtered[!is.na(datHuntLabelledEventsSBMerged_filtered$groupID) & 
                                                                                !(datHuntLabelledEventsSBMerged_filtered$expID %in% vxCludeExpID),]
for (dID in vWeirdDataSetID )
  print(NROW(unique(datHuntLabelledEventsSBMerged_fixed[datHuntLabelledEventsSBMerged_fixed$dataSetID ==  dID ,]$expID)))

################# # ## # # 



## Get Summarized Hunt Results Per Larva ####
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged_fixed)


## Get Event Counts Within Range  - Along With Total Number of Hunting frames for each Larva##
## Added Larva ID to Check for Correlation Through Time of Day - Surrogate as LarvaID;s increased through the day of the experiment from 1-4
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LL ),datHuntStat[,"vIDLookupTable"]$LL$larvaID )
datHuntVsPreyLL <- datHuntVsPreyLL[!is.na(datHuntVsPreyLL[,1]) ,]
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$LE ),datHuntStat[,"vIDLookupTable"]$LE$larvaID  )
datHuntVsPreyLE <- datHuntVsPreyLE[!is.na(datHuntVsPreyLE[,1]) ,]



datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NL),datHuntStat[,"vIDLookupTable"]$NL$larvaID )
datHuntVsPreyNL <- datHuntVsPreyNL[!is.na(datHuntVsPreyNL[,1]) ,]
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$NE),datHuntStat[,"vIDLookupTable"]$NE$larvaID  )
datHuntVsPreyNE <- datHuntVsPreyNE[!is.na(datHuntVsPreyNE[,1]) ,]

datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DL ),datHuntStat[,"vIDLookupTable"]$DL$larvaID  )
datHuntVsPreyDL <- datHuntVsPreyDL[!is.na(datHuntVsPreyDL[,1]),] ##Remove NA And High Fliers
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),as.numeric(datHuntStat[,"vHDurationPerLarva"]$DE ),datHuntStat[,"vIDLookupTable"]$DE$larvaID  )
datHuntVsPreyDE <- datHuntVsPreyDE[!is.na(datHuntVsPreyDE[,1]),] ##Remove NA And High Fliers


## Check Number of Hunt Events For A LarvaID 
vEventCount <- vector()
datHuntVsPreyC <- datHuntVsPreyDL ##Check Which Group
for (i in 1:4)
{
  vEventCount[i] <- ( sum(datHuntVsPreyC[ datHuntVsPreyC[,4] == i & !is.na(datHuntVsPreyC[,4]) ,2] ) )
}
print(vEventCount)

### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)

### Rum The Sampler ###
drawLL2 <- mcmc_drawEventCountModels(datHuntVsPreyLL,preyCntRange,"modelGroupEventRate.tmp")
drawNL2 <- mcmc_drawEventCountModels(datHuntVsPreyNL,preyCntRange,"modelGroupEventRate.tmp")
drawDL2 <- mcmc_drawEventCountModels(datHuntVsPreyDL,preyCntRange,"modelGroupEventRate.tmp")
drawLE2 <- mcmc_drawEventCountModels(datHuntVsPreyLE,preyCntRange,"modelGroupEventRate.tmp")
drawNE2 <- mcmc_drawEventCountModels(datHuntVsPreyNE,preyCntRange,"modelGroupEventRate.tmp")
drawDE2 <- mcmc_drawEventCountModels(datHuntVsPreyDE,preyCntRange,"modelGroupEventRate.tmp")

save(drawLL2,drawNL2,drawDL2,drawLE2,drawNE2,drawDE2,file =paste(strDataExportDir,"stat_HuntRateInPreyRange_nbinomRJags.RData",sep=""))


### Draw Distribution oF Hunt Rates - 
## for the exp draw (z= p/(1-p)) ## But it is the same for Rate Of Gamma Too / Or inverse for scale
plotsamples <- 500
schain <-1:3

### The Prob Of Success p from NegBinom translates to Gamma Rate p/(1-p), or scale: (1-p)/p
HEventHuntGammaRate_LE <-tail(drawLE2$q[,,schain],plotsamples)/(1-tail(drawLE2$q[,,schain],plotsamples));
HEventHuntGammaRate_LL <-tail(drawLL2$q[,,schain],plotsamples)/(1-tail(drawLL2$q[,,schain],plotsamples));
HEventHuntGammaRate_DE <- tail(drawDE2$q[,,schain],plotsamples)/(1-tail(drawDE2$q[,,schain],plotsamples));
HEventHuntGammaRate_DL <- (tail(drawDL2$q[,,schain],plotsamples)/(1-tail(drawDL2$q[,,schain],plotsamples)));      
HEventHuntGammaRate_NE <- (tail(drawNE2$q[,,schain],plotsamples)/(1-tail(drawNE2$q[,,schain],plotsamples)));
HEventHuntGammaRate_NL <- (tail(drawNL2$q[,,schain],plotsamples)/(1-tail(drawNL2$q[,,schain],plotsamples)));      
HEventHuntGammaShape_LE <- tail(drawLE2$r[,,schain],plotsamples);
HEventHuntGammaShape_LL <- tail(drawLL2$r[,,schain],plotsamples)
HEventHuntGammaShape_DE <- tail(drawDE2$r[,,schain],plotsamples);
HEventHuntGammaShape_DL <- tail(drawDL2$r[,,schain],plotsamples)
HEventHuntGammaShape_NE <- tail(drawNE2$r[,,schain],plotsamples);
HEventHuntGammaShape_NL <- tail(drawNL2$r[,,schain],plotsamples)

pBW = 0.5
densHPoissonRate_LL <- density( HEventHuntGammaShape_LL*1/HEventHuntGammaRate_LL,bw=pBW)
densHPoissonRate_LE <- density( HEventHuntGammaShape_LE*1/HEventHuntGammaRate_LE,bw=pBW)
densHPoissonRate_DE <- density( HEventHuntGammaShape_DE*1/HEventHuntGammaRate_DE,bw=pBW)
densHPoissonRate_DL <- density( HEventHuntGammaShape_DL*1/HEventHuntGammaRate_DL,bw=pBW)
densHPoissonRate_NE <- density( HEventHuntGammaShape_NE*1/HEventHuntGammaRate_NE,bw=pBW)
densHPoissonRate_NL <- density( HEventHuntGammaShape_NL*1/HEventHuntGammaRate_NL,bw=pBW)


Plim <- max(range(datHuntVsPreyLL[,2])[2],range(datHuntVsPreyDL[,2])[2],range(datHuntVsPreyNL[,2])[2])
x <- seq(0.01,Plim,0.05)

#### HUNT EVENT PER LARVA PLOT #####
## Comprehensive Plot On Number of Hunt Events
pdf(file= paste(strPlotExportPath,"/stat/fig2_statComparePoissonHuntRates",".pdf",sep=""),width = 14,height = 7)
##Now Plot Infered Distributions
##Show Alignment with Empirical Distribution of HuntEvent Numbers
## Number of Hunt Events Per Larva
outer = FALSE
line = 1 ## SubFig Label Params
cex = 1.1
adj  = 2.5
padj <- -10.5
las <- 1

layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2,6, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,3.3,1,1))
#lineTypeL[1] <- 1
plotEventCountDistribution_cdf(datHuntVsPreyNE,drawNE2,colourHE[1],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE )
plotEventCountDistribution_cdf(datHuntVsPreyNL,drawNL2,colourHL[1],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE )
legend("bottomright",legend = c(  expression (),
                                  bquote(NF["s"] ~ 'Data #' ~ .(NROW(datHuntVsPreyNE))  )
                                 ,bquote(NF["s"] ~"Model " ), 
                                 bquote(NF["e"] ~ 'Data #' ~ .(NROW(datHuntVsPreyNL)) )
                                 ,bquote( NF["e"] ~ "Model " ) ), 
       col=c(colourP[4], colourLegE[1],colourP[4],colourLegL[1]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext(side = 1,cex=0.8, line = 2.2, "Number of hunt events ")
mtext(side = 2,cex=0.8, line = 2.2, " Cumulative function ")


plotEventCountDistribution_cdf(datHuntVsPreyLE,drawLE2,colourHE[2],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE )
plotEventCountDistribution_cdf(datHuntVsPreyLL,drawLL2,colourHL[2],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE )
legend("bottomright",legend = c(  expression (),
                                  bquote(LF["s"] ~ 'Data #' ~ .(NROW(datHuntVsPreyLE))  )
                                  ,bquote(LF["s"] ~"Model " ), 
                                  bquote(LF["e"] ~ 'Data #' ~ .(NROW(datHuntVsPreyLL)) )
                                  ,bquote( LF["e"] ~ "Model " ) ), 
      col=c(colourP[4], colourLegE[2],colourP[4],colourLegL[2]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext(side = 1,cex=0.8, line = 2.2, "Number of hunt events ")
mtext(side = 2,cex=0.8, line = 2.2, " Cumulative function ")


plotEventCountDistribution_cdf(datHuntVsPreyDE,drawDE2,colourHE[3],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE  )
plotEventCountDistribution_cdf(datHuntVsPreyDL,drawDL2,colourHL[3],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE  )
legend("bottomright",legend = c(  expression (),
                                  bquote(DF["s"] ~ 'Data #' ~ .(NROW(datHuntVsPreyDE))  )
                                  ,bquote(DF["s"] ~"Model " ), 
                                  bquote(DF["e"] ~ 'Data #' ~ .(NROW(datHuntVsPreyDL)) )
                                  ,bquote( DF["e"] ~ "Model " ) ), 
       col=c(colourP[4], colourLegE[3],colourP[4],colourLegL[3]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
mtext(side = 1,cex=0.8, line = 2.2, "Number of hunt events ")
mtext(side = 2,cex=0.8, line = 2.2, " Cumulative function ")
mtext("C",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)



pchL <- c(1,2,0,16,17,15)
### Plot GAMMA Parameters Space ####
#Xlim <- 50
#plot(1/HEventHuntGammaRate_NE,HEventHuntGammaShape_NE,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[1],xlab=NA,ylab=NA)
#points(1/HEventHuntGammaRate_LE,HEventHuntGammaShape_LE,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[2])
#points(1/HEventHuntGammaRate_DE,HEventHuntGammaShape_DE,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[3])
#points(1/HEventHuntGammaRate_NL,HEventHuntGammaShape_NL,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[4])
#points(1/HEventHuntGammaRate_LL,HEventHuntGammaShape_LL,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[5])
#points(1/HEventHuntGammaRate_DL,HEventHuntGammaShape_DL,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[6])
#mtext(side = 1,cex=0.8, line = 2.2, expression(paste(Gamma, " scale (r)") ) )
#mtext(side = 2,cex=0.8, line = 2.2, expression(paste(Gamma, " shape (k)") ) )
#legend("topright",legend = strDataLabels,  #c(paste("NE" ), paste("LE"),paste("DE"), paste("NL"),paste("LL"),paste("DL")),
#       col=c(colourHL[1],colourHL[2],colourHL[3],colourHL[1],colourHL[2],colourHL[3]) ,pch=pchL,cex=0.9,bg="white",ncol=2)
#mtext("D",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

###

## BoxPlot of Hunt Event Counts - 
strCondTags <- c("NE","NL","LE","LL","DE","DL")
xbarcenters <- boxplot(log10(datHuntVsPreyNE[,2]+1),log10(datHuntVsPreyNL[,2]+1),log10(datHuntVsPreyLE[,2]+1),log10(datHuntVsPreyLL[,2]+1),log10(datHuntVsPreyDE[,2]+1),log10(datHuntVsPreyDL[,2]+1),
        main=NA,notch=TRUE,col=colourD,names=strCondTags,ylim=c(0,2),axes = FALSE  )
mtext(side = 2,cex=0.8, line =2.2, "Number of hunt events " ) #log(N+1)
vIDTable    <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]
vDat        <- (datHuntStat[,"vHLarvaEventCount"])

axis(1,at<-axis(1,labels=NA), labels=c( strDataLabels[1],strDataLabels[4],strDataLabels[2],strDataLabels[5],strDataLabels[3],strDataLabels[6] ))
yticks <-axis(2,labels=NA)
axis(2, at = yticks, labels =round(10^yticks)-1 , col.axis="black", las=2)
#minor.tick() ##minor.tick(nx=4, ny=4, tick.ratio=4,x.args=NA) ## Can COnfuse to think scale is linear
## Connect Larvae From EMpty To LIve Test Condition #
plotConnectedEventCounts(datHuntStat,vDat,strCondTags)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)


## Compare Distrib Of Poisson Rate Derived from Model  Params ##
schain <- 1:3
Ylim <- 3 
pBW <- 0.001


Ylim <- 0.7
plot(densHPoissonRate_NL$x, densHPoissonRate_NL$y,type='l',lty=lineTypeL[2],col=colourHLine[1] ,lwd=4,ylab=NA,xlab=NA,xlim=c(0,25),ylim=c(0,Ylim))
lines(densHPoissonRate_LL$x, densHPoissonRate_LL$y,type='l',lty=lineTypeL[2],col=colourHLine[2],lwd=4,ylab=NA,xlab=NA)
lines(densHPoissonRate_DL$x, densHPoissonRate_DL$y,type='l',lty=lineTypeL[2],col=colourHLine[3],lwd=4,ylab=NA,xlab=NA)

lines(densHPoissonRate_NE$x, densHPoissonRate_NE$y,type='l',lty=lineTypeL[1],col=colourHLine[1],lwd=4,ylab=NA,xlab=NA)
lines(densHPoissonRate_LE$x, densHPoissonRate_LE$y,type='l',lty=lineTypeL[1],col=colourHLine[2],lwd=4,ylab=NA,xlab=NA)
lines(densHPoissonRate_DE$x, densHPoissonRate_DE$y,type='l',lty=lineTypeL[1],col=colourHLine[3],lwd=4,ylab=NA,xlab=NA)

legend("topright",legend = c(paste("Spontaneous " ),paste("Evoked ")),seg.len=3.5
       , col=c(colourR[4], colourR[4]),lty=c(2,1),lwd=4,cex=1.1,bg="white" )
mtext(side = 1,cex=0.8, line = 2.5, expression(paste("Estimated hunt rate  (",lambda," )") )  )
mtext(side = 2,cex=0.8, line = 2.2, " Density function ")
mtext("E",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)

dev.off() 


################## ############# ### # # 
# legend("topright",legend = c(paste("NF " ),paste("LF ") ,paste("DF ") ) ,
#        col=colourL,lty=lineTypeL,lwd=3,seg.len=3,cex=1.1)


##Debug
#plot(x,log(dgamma(x,shape=HEventHuntGammaShape_LL[1,1],rate=HEventHuntGammaRate_LL[1,1])) ,main="",xlab=NA,ylab=NA,col=colourHL[2],type="l",lwd=1, lty=1)
#lines(x,log(dgamma(x,shape=HEventHuntGammaShape_LE[1,1],rate=HEventHuntGammaRate_LE[1,1])) ,main="",xlab=NA,ylab=NA,col=colourHL[2],type="l",lwd=1, lty=2)



#### Plot Gamma Distributions Of Hunt Rates From Which Hunt Event Counts where Drawn 
#plotGammaHuntRates(x,HEventHuntGammaShape_NE,HEventHuntGammaRate_NE,colourHE[1],plotsamples,2,TRUE)
#plotGammaHuntRates(x,HEventHuntGammaShape_NL,HEventHuntGammaRate_NL,colourHL[1],plotsamples,1,FALSE)
#legend("topright",legend = c(paste("Spontaneous (NE)" ),paste("Evoked (NL)")), col=c(colourHE[1], colourHL[1]),lty=c(2,1),lwd=2,cex=1.1,bg="white" )
#plotGammaHuntRates(x,HEventHuntGammaShape_LE,HEventHuntGammaRate_LE,c(colourHE[2]),plotsamples,2,TRUE)
#plotGammaHuntRates(x,HEventHuntGammaShape_LL,HEventHuntGammaRate_LL,c(colourHL[2]),plotsamples,1,FALSE)
#legend("topright",legend = c(paste("Spontaneous (LE)" ),paste("Evoked (LL)")), col=c(colourHE[2], colourHL[2]),lty=c(2,1),lwd=2,cex=1.1,bg="white" )
#plotGammaHuntRates(x,HEventHuntGammaShape_DE,HEventHuntGammaRate_DE,c(colourHE[3]),plotsamples,2,TRUE)
#plotGammaHuntRates(x,HEventHuntGammaShape_DL,HEventHuntGammaRate_DL,c(colourHL[3]),plotsamples,1,FALSE)


# ## Compare Distrib Of Model Params ##
# schain <- 1:3
# Ylim <- 50 
# pBW <- 0.001
# densHEventSampled_LL <-density(drawLL2$q[1,,schain],bw=pBW);densHEventSampled_DL <-density(drawDL2$q[1,,schain],bw=pBW);densHEventSampled_NL <-density(drawNL2$q[1,,schain],bw=pBW)
# densHEventSampled_LE <-density(drawLL2$q[1,,schain],bw=pBW);densHEventSampled_DE <-density(drawDL2$q[1,,schain],bw=pBW);densHEventSampled_NE <-density(drawNL2$q[1,,schain],bw=pBW)   
# 
# plot(densHEventSampled_NL$x, densHEventSampled_NL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[1],col=colourL[1],lwd=4,ylab=NA,xlab=NA)
# lines(densHEventSampled_LL$x, densHEventSampled_LL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[2],col=colourL[2],lwd=4,ylab=NA,xlab=NA)
# lines(densHEventSampled_DL$x, densHEventSampled_DL$y,type='l',xlim=c(0,0.4),ylim=c(0,Ylim),lty=lineTypeL[3],col=colourL[3],lwd=4,ylab=NA,xlab=NA)
# mtext(side = 1,cex=0.8, line =2.2, expression(paste("Model  Parameter (p",") ") ) )
# mtext(side = 2,cex=0.8, line =2.2, paste("Density Bw:",densHEventSampled_LL$bw))
# legend("topright",legend = c(paste("NF " ),paste("LF ") ,paste("DF ") ) ,
#        col=colourL,lty=lineTypeL,lwd=3,seg.len=3,cex=1.1)




#points(x/G_APPROXFPS,dgamma(x,rate=drawDurDE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurDE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="DE",xlim=c(0,6),col=colourR[1] ) 
#points(x/G_APPROXFPS,dgamma(x,rate=drawDurNE$r[,(steps/thin-ns):(steps/thin),],shape=drawDurNE$s[,(steps/thin-ns):(steps/thin),]),type="p",pch=16,cex=0.4,main="NE",xlim=c(0,6),col=colourR[3] ) 

## Sand Box ##

labexpID <- table(datHuntLabelledEventsSBMerged[datHuntLabelledEventsSBMerged$groupID == "NL",]$expID)
regexpID <- table(datTrackedEventsRegister[datTrackedEventsRegister$groupID == "NL" ,]$expID) 

datHuntLabelledEventsSBMerged_filtered[datHuntLabelledEventsSBMerged_filtered$groupID == "NL" & datHuntLabelledEventsSBMerged_filtered$expID == 3911,]

datTrackedEventsRegister[datTrackedEventsRegister$groupID == "NL" & datTrackedEventsRegister$expID == 3911 ,]


#### DEBUG CODE - Fitting An EXP Distribution TEST ###
f1L <- fitdist( datJagsLE$d,"exp")
plot(f1L,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f1D <- fitdist( datJagsDE$d,"exp");
plot(f1D,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f1N <- fitdist(  datJagsNE$d,"exp");
plot(f1N,sub="exp") ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
##### ### 
f2D <- fitdist( datJagsDE$d,"exp",lower = c(0, 0),start = list(size = 10,mu=30))
plot(f2D) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f2N <- fitdist( datHuntVsPreyNL[,2],"nbinom",lower = c(0, 0)) ##,start = list(scale = 1, shape = 1)
plot(f2N) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
f2L <- fitdist( datHuntVsPreyLL[,2],"nbinom",lower = c(0, 0)) #,start = list(scale = 1, shape = 1)
plot(f2L) #,start = list(scale = 1, shape = 1) ##Show Diagnostics Of Fitting A EXP using Standard Methods ##
## ########
