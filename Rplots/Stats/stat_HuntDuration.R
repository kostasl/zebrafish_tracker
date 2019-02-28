### KOstasl 2019

### Produces Figure that compares the statistics of total Larva HUnt Duration for each experiment in a group ###
### 



#####
##Models Each Event in the population of larvaer Individually 
modelEventDuration="model { 

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
strModelName = "modelLarvaEventDuration.tmp"
fileConn=file(strModelName)
writeLines(modelEventDuration,fileConn);
close(fileConn)

## Discrete - Geometric Cause Mixture of rates - assuming rates drawn from most informative Prior distribution (EXP)
## Assuming nbinom(r,p) Poisson(L|a,b) Gamma(a,b) then r=a, p=1/(b+1) -> b=(1-p)/p
## Give Neg Binomial
modelLarvaHuntDuration="model { 
q ~ dunif(0.0,1)
r ~ dgamma(1,1)

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
plotHuntDurationDistribution_cdf <- function(datHDuration,drawHEvent,lcolour,lpch,lty,Plim,nplotSamples=100,newPlot = FALSE)
{
  XLim <- G_APPROXFPS*180
  x <- seq(0,XLim,1)
  
  cdfD_N <- ecdf(datHDuration[,3]/G_APPROXFPS)
  
  plot(cdfD_N,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0, XLim/G_APPROXFPS),
       ylim=c(0,1),cex=1.5,cex.lab=1.5,add=!newPlot)
  ##Construct CDF of Model by Sampling randomly from Model distribution for exp rate parameter
  for (c in 1:NROW(drawHEvent$q[1,1,])) {
    for (j in (NROW(drawHEvent$q[,,c])-nplotSamples):NROW(drawHEvent$q[,,c]) )
    {
      cdfM <- dnbinom(x,size=drawHEvent$r[,j,c],prob=  drawHEvent$q[,j,c]  )##1-exp(-q*x) ##ecdf(  dexp( x, q  ) )
      lines(x/G_APPROXFPS,cumsum(cdfM),col=lcolour,lty=lty) #add=TRUE,
    }
  }
  plot(cdfD_N ,col=colourP[4],pch=lpch,xlab=NA,ylab=NA,main="",xlim=c(0,XLim/G_APPROXFPS),ylim=c(0,1),cex=1.1,cex.lab=1.5,add=TRUE)
  
}


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
      
      segments(gIdx,log10(ptSrc+1),gIdx+1,log10(ptDest+1) ,col=colourP[4])
      
      
    }##For Each Experiment
    
  } ## Go Through Pairs Of Conditions ##
}##End Of Function ConnectEventPoints


##Random Init Of Chain 
initLarvaHuntDurfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(r=abs(rnorm(N,3,1)), ##Base Line Vergence Prior to HuntOn
                                     s=abs(rnorm(N,3,1)) ),
                        simplify=FALSE)
  
  return(initlist)
}

############ LOAD DATA #################

############ LOAD EVENTS LIst and Fix ####
## Warning Set Includes Repeated Test For some LF fish - One In Different Food Density
## Merged2 Contains the Fixed, Remerged EventID 0 files, so event Counts appear for all larvae recorded.
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged2" 

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
  
for (dID in vWeirdDataSetID )
  print(NROW(unique(datHuntLabelledEventsSBMerged_fixed[datHuntLabelledEventsSBMerged_fixed$dataSetID ==  dID ,]$expID)))

datHuntLabelledEventsSBMerged_fixed <- datHuntLabelledEventsSBMerged_filtered[!is.na(datHuntLabelledEventsSBMerged_filtered$groupID) & 
                                                                                !(datHuntLabelledEventsSBMerged_filtered$expID %in% vxCludeExpID),]

################# # ## # # 



## Get Summarized Hunt Results Per Larva ####
datHuntStat <- makeHuntStat(datHuntLabelledEventsSBMerged_fixed)

## Baysian Inference Fiting a Gamma distribution to the Hunt Event Duration Data ##
##Setup Data Structure To Pass To RJAgs

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


plotsamples = 20
##Remove Rec Of No Duration ###
### The Hunt Duration Is Discontinuous if we include the Larvae With 0 events, and this creates modelling problems 
## Thus when Looking for Hunt Duration per larva we need to exclude the ones that did not hunt.
datHuntVsPreyLE <- datHuntVsPreyLE[datHuntVsPreyLE[,2] > 0,]
datHuntVsPreyNE <- datHuntVsPreyNE[datHuntVsPreyNE[,2] > 0,]
datHuntVsPreyDE <- datHuntVsPreyDE[datHuntVsPreyDE[,2] > 0,]
datHuntVsPreyLL <- datHuntVsPreyLL[datHuntVsPreyLL[,2] > 0,]
datHuntVsPreyNL <- datHuntVsPreyNL[datHuntVsPreyNL[,2] > 0,]
datHuntVsPreyDL <- datHuntVsPreyDL[datHuntVsPreyDL[,2] > 0,]


load(file =paste(strDataExportDir,"stat_HuntDurationInPreyRange_nbinomRJags.RData",sep=""))

drawDurLE <- mcmc_drawHuntDurationModels(datHuntVsPreyLE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurNE <- mcmc_drawHuntDurationModels(datHuntVsPreyNE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurDE <- mcmc_drawHuntDurationModels(datHuntVsPreyDE,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurLL <- mcmc_drawHuntDurationModels(datHuntVsPreyLL,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurDL <- mcmc_drawHuntDurationModels(datHuntVsPreyDL,preyCntRange,"modelLarvaHuntDuration.tmp" )
drawDurNL <- mcmc_drawHuntDurationModels(datHuntVsPreyNL,preyCntRange,"modelLarvaHuntDuration.tmp" )

save(drawDurLE,drawDurNE,drawDurDE,drawDurLL,drawDurNL,drawDurDL,
     file =paste(strDataExportDir,"stat_HuntDurationInPreyRange_nbinomRJags.RData",sep=""))


plotsamples <- 15
schain <-1:3

### The Prob Of Success p from NegBinom translates to Gamma Rate p/(1-p), or scale: (1-p)/p
HEventHuntGammaRate_LE <-tail(drawDurLE$q[,,schain],plotsamples)/(1-tail(drawDurLE$q[,,schain],plotsamples));
HEventHuntGammaRate_LL <-tail(drawDurLL$q[,,schain],plotsamples)/(1-tail(drawDurLL$q[,,schain],plotsamples));
HEventHuntGammaRate_DE <- tail(drawDurDE$q[,,schain],plotsamples)/(1-tail(drawDurDE$q[,,schain],plotsamples));
HEventHuntGammaRate_DL <- (tail(drawDurDL$q[,,schain],plotsamples)/(1-tail(drawDurDL$q[,,schain],plotsamples)));      
HEventHuntGammaRate_NE <- (tail(drawDurNE$q[,,schain],plotsamples)/(1-tail(drawDurNE$q[,,schain],plotsamples)));
HEventHuntGammaRate_NL <- (tail(drawDurNL$q[,,schain],plotsamples)/(1-tail(drawDurNL$q[,,schain],plotsamples)));      
HEventHuntGammaShape_LE <- tail(drawDurLE$r[,,schain],plotsamples);
HEventHuntGammaShape_LL <- tail(drawDurLL$r[,,schain],plotsamples)
HEventHuntGammaShape_DE <- tail(drawDurDE$r[,,schain],plotsamples);
HEventHuntGammaShape_DL <- tail(drawDurDL$r[,,schain],plotsamples)
HEventHuntGammaShape_NE <- tail(drawDurNE$r[,,schain],plotsamples);
HEventHuntGammaShape_NL <- tail(drawDurNL$r[,,schain],plotsamples)



layout(matrix(c(1,2,3,4,5,6), 3,2, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,3.3,1,1))
plotHuntDurationDistribution_cdf(datHuntVsPreyLE,drawDurLE,colourHE[2],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntVsPreyLL,drawDurLL,colourHL[2],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)
legend("bottomright",legend = c(paste("Data LE #",NROW(datHuntVsPreyLE) ),paste("Model LE "),
                                paste("Data LL #",NROW(datHuntVsPreyLL) ),paste("Model LL ")), 
       col=c(colourP[4], colourLegE[2],colourP[4],colourLegL[2]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )


plotHuntDurationDistribution_cdf(datHuntVsPreyNE,drawDurNE,colourHE[1],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntVsPreyNL,drawDurNL,colourHL[1],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)
legend("bottomright",legend = c(paste("Data NE #",NROW(datHuntVsPreyNE) ),paste("Model NE "),
                                paste("Data NL #",NROW(datHuntVsPreyNL) ),paste("Model NL ")), 
       col=c(colourP[4], colourLegE[1],colourP[4],colourLegL[1]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )

plotHuntDurationDistribution_cdf(datHuntVsPreyDE,drawDurDE,colourHE[3],pchL[1],lineTypeL[2],Plim,plotsamples,newPlot=TRUE)
plotHuntDurationDistribution_cdf(datHuntVsPreyDL,drawDurDL,colourHL[3],pchL[3],lineTypeL[2],Plim,plotsamples,newPlot=FALSE)

legend("bottomright",legend = c(paste("Data DE #",NROW(datHuntVsPreyDE) ),paste("Model DE "),
                                paste("Data DL #",NROW(datHuntVsPreyDL) ),paste("Model DL ")), 
       col=c(colourP[4], colourLegE[3],colourP[4],colourLegL[3]), pch=c(pchL[1],NA,pchL[3],NA),lty=c(NA,1),lwd=2,cex=1.1,bg="white" )
mtext(side = 1,cex=0.8, line = 2.2, " Hunt Duration per Larva (sec)")
mtext(side = 2,cex=0.8, line = 2.2, " F(x < N) ")

######
###### Gamma Parameters Comparison###

pchL <- c(1,2,0,16,17,15)
### Plot GAMMA Parameters Space
Xlim <- range(1/HEventHuntGammaRate_LL/G_APPROXFPS)[2]
plot(1/HEventHuntGammaRate_NE/G_APPROXFPS,HEventHuntGammaShape_NE,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[1],xlab=NA,ylab=NA)
points(1/HEventHuntGammaRate_LE/G_APPROXFPS,HEventHuntGammaShape_LE,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[2])
points(1/HEventHuntGammaRate_DE/G_APPROXFPS,HEventHuntGammaShape_DE,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[3])
points(1/HEventHuntGammaRate_NL/G_APPROXFPS,HEventHuntGammaShape_NL,col=colourHL[1],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[4])
points(1/HEventHuntGammaRate_LL/G_APPROXFPS,HEventHuntGammaShape_LL,col=colourHL[2],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[5])
points(1/HEventHuntGammaRate_DL/G_APPROXFPS,HEventHuntGammaShape_DL,col=colourHL[3],ylim=c(0,3),xlim=c(0,Xlim),pch=pchL[6])
strXLab <- (expression(paste(Gamma, " scale (r/FPS)") ) )
mtext(side = 1,cex=0.8, line = 2.2,strXLab ) 
mtext(side = 2,cex=0.8, line = 2.2, expression(paste(Gamma, " shape (k)") ) )
legend("topright",legend = c(paste("NE" ),
                             paste("LE"),paste("DE"), paste("NL"),paste("LL"),paste("DL")),
       col=c(colourHL[1],colourHL[2],colourHL[3],colourHL[1],colourHL[2],colourHL[3]) ,pch=pchL,cex=0.9,bg="white",ncol=2)
###

## BoxPlot of Hunt Event Counts - 
strCondTags <- c("NE","NL","LE","LL","DE","DL")
xbarcenters <- boxplot(log10( (datHuntVsPreyNE[,3]+1)/G_APPROXFPS ) ,log10( ( datHuntVsPreyNL[,3]+1)/G_APPROXFPS ),log10( (datHuntVsPreyLE[,3]+1)/G_APPROXFPS ),
                       log10( ( datHuntVsPreyLL[,3]+1 )/G_APPROXFPS ),log10(( datHuntVsPreyDE[,3]+1)/G_APPROXFPS ) ,log10( ( datHuntVsPreyDL[,3]+1)/G_APPROXFPS ),
                       main=NA,notch=TRUE,col=colourD,names=strCondTags,ylim=c(0,2),axes = FALSE  )
mtext(side = 2,cex=0.8, line =2.2, "Hunt Duration log(D+1) ")
vIDTable    <- datHuntStat[,"vIDLookupTable"] ##vIDTable$DL <- vIDTable$DL[vIDTable$DL$expID!=3830,]
vDat        <- (datHuntStat[,"vHLarvaEventCount"])

axis(1,at<-axis(1,labels=NA), labels=strCondTags)
yticks <-axis(2,labels=NA)
axis(2, at = yticks, labels =round(10^yticks) , col.axis="black", las=2)
## Connect Larvae From EMpty To LIve Test Condition #
plotConnectedEventCounts(datHuntStat,strCondTags)

###########
##CHECK fIT
plot(1:15000,dnbinom(1:15000, size=tail(drawDurLE$r,50),  prob=tail(drawDurLE$q,50)),cex=0.3,col=colourHE[2] )
lines(density(drawDurLE$d,bw=1000 ) )
plot(1:15000,dnbinom(1:15000, size=tail(drawDurNE$r,50),  prob=tail(drawDurNE$q,50)),cex=0.3,col=colourHE[1] )
lines(density(drawDurNE$d,bw=1000 ) )
plot(1:15000,dnbinom(1:15000, size=tail(drawDurDE$r,50),  prob=tail(drawDurDE$q,50)),cex=0.3,col=colourHE[3] )
lines(density(drawDurDE$d,bw=1000 ) )

plot(1:15000,dnbinom(1:15000, size=tail(drawDurNL$r,50),  prob=tail(drawDurNL$q,50)),cex=0.3,col=colourHL[1] )
lines(density(drawDurNL$d,bw=1000 ) )
plot(1:15000,dnbinom(1:15000, size=tail(drawDurLL$r,50),  prob=tail(drawDurLL$q,50)),cex=0.3,col=colourHL[2] )
lines(density(drawDurLL$d,bw=1000 ) )

#####################################################   ############### ######################################

############### Model Each Larva s Event Duration Separatelly ################################################
###            Use The Labelled Events Register               #### ##################
## Run Baysian Inference on Model for Hunt Event Counts IN A Group/ Test Condition
## Return Samples Drawn structure                             ########
####### Function Returns Hunt Event Durations for Group ID, excluding events 0 (Food Count Event) 
getdatHuntEventDuration <- function(strGroupID)
{
  
  datDurationPerEpisodePerLarva <- (with(datHuntLabelledEventsSBMerged_fixed,
                                         data.frame(DurationFrames=endFrame[groupID == strGroupID & eventID != 0]-startFrame[groupID == strGroupID & eventID != 0],
                                                    expID=expID[groupID == strGroupID & eventID != 0],
                                                    hidx=as.numeric(factor(expID[groupID == strGroupID & eventID != 0]))) )) ##Add hidx to use for correct prior Init in RJags
  
  datDurationPerEpisodePerLarva
  
  return (datDurationPerEpisodePerLarva)
}

##Random Init Of Chain 
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
  plotsamples = 200
  thin=2;
  chains = 3
  
  
  ##Larva Event Counts Slice
  datSliceF <- datHuntVsPrey[datHuntVsPrey[,1] >= preyCntRange[1] & datHuntVsPrey[,1] <= preyCntRange[2], ]
  
  datJags=list(d=datHuntVsPrey[,3],NTOT=NROW(datHuntVsPrey));
  
  nDat = NROW(datHuntVsPrey);
  dataG=list(d=datHuntVsPrey$DurationFrames,NTOT=nDat,hidx=datHuntVsPrey$hidx,food=as.integer(datSliceF[,1]));
  
  model=jags.model(file=strModelFilename,data=dataG,n.chains=chains);
  
  update(model,burn_in)
  
  drawSamples=jags.samples(model,steps,thin=thin,variable.names=varnames1)
  
  return(drawSamples) 
}


datHE_LE <- getdatHuntEventDuration("LE")
##
### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)

drawHD_LE<-mcmc_drawEventDurationModels (datHE_LE,preyCntRange, "modelLarvaEventDuration.tmp")

hist(drawDurLE$q)
hist(drawDurN$s)
hist(drawDurD$s)

plot(tail(drawDurLE$r,1000), tail(drawDurLE$s,1000),ylim=c(0,200))
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
