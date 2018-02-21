##plot Hunting Event Statistics


## Box Plots Used to Compare Conditions On Mean Stats - Saves Output As Pdf
boxplotPerCondition <- function(datStat,datMean,datSe,strtitle,strsubt,stroutFileName,plotTop)
{
  
  par(mar = c(5, 6, 4, 5) + 2.5)

  datN = vector()
  if (is.null(dim(datStat))) ##This is probably A List
  {
    datlbls <-names(datStat) ##Input is List
    for (i in seq(1:NROW(datStat)))
    {
      datN[i] <-NROW( datStat[[i]][!is.na(datStat[[i]]) ]   )
    }
    # datN <- NROW(datStat[[1]])
  }
  else ##A Known Data Frame Struct
  { 
    datN <- unlist(datStat[,"nLarva"],use.names = FALSE)
    datlbls <-row.names(datStat) ##Multidim Use Row Names
  }
  ##Add N Numbers to Labels
  datlbls <- paste(datlbls,"\nn=",datN,sep="")
  
  if(missing(plotTop)) 
  {
    plotTop <- max(datMean) +
      unique(datSe[datMean== max(datMean)]) * 3
  }else
  {
    plotTop <- 1.1*plotTop ##Increase by 10%
    
  }
  plotBottom <- 2*min(min(datMean),0)
  
  barCenters <- barplot(height = datMean,
                        names.arg = datlbls,
                        beside = true, las = 2,
                        ylim = c(plotBottom, plotTop),
                        cex.names = 0.75, xaxt = "n",
                        main = strtitle,
                        sub = strsubt,
                        ylab = "#",
                        border = "black", axes = TRUE)
  
  # Specify the groupings. We use srt = 45 for a
  # 45 degree string rotation
  #text(x = barCenters, y = par("usr")[3] - 0.01, srt = 45,    adj = 1, labels = datlbls, xpd = TRUE)
  text(x = barCenters+0.2, y = 0, srt = 45,
       adj = 1.8, labels = datlbls, xpd = TRUE)
  
  segments(barCenters, datMean - datSe * 2, barCenters,
           datMean + datSe * 2, lwd = 1.5)
  #title(sub=strsubt)
  #  dev.off()
  
  message(strsubt)
  
  return(barCenters)
}




plotConnectedPointsPairs <- function(vIDTable,vDat,strCondTags,xbarcenters)
{
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    
    datSetID <- levels(vIDTable[[gE]]$dataSetID)[ vIDTable[[gE]]$dataSetID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID] ]

    idsE <- vIDTable[[gE]]$expID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID ==vIDTable[[gL]]$dataSetID]
    idsL <- vIDTable[[gL]]$expID[vIDTable[[gL]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID]
    ptSrc  <- vDat[[gE]][levels(idsE)[idsE]]
    ptDest <- vDat[[gL]][levels(idsL)[idsL]]
    
    ##OPTIONAL: Replace NAs With 0 when Plotting So as to To Show Direction For All points
    #ptSrc <- replace(ptSrc,is.na(ptSrc),0)
    #ptDest <- replace(ptDest,is.na(ptDest),0)
    ##Otherwise Do not Plot NA connecting line
    
    ##Plot The Lines Connect Each Empty Tested Larva With Itself In THe Live Fed Conditions 
    idxSrc  <- match(gE,strCondTags) ##Bar Center Idx for Each Condition E. Fed
    idxDest <- match(gL ,strCondTags)           
    
    for (p in 1:NROW(ptSrc))
    {
      pcolour <- rDataset[as.numeric( replace( datSetID[p],is.na(datSetID[p]) ,1 ) ) ]
      points(xbarcenters[idxSrc],ptSrc[p],pch=as.numeric(datSetID[p]),
             col=pcolour )
      points(xbarcenters[idxDest],ptDest[p],pch=as.numeric(datSetID[p]),
             col=pcolour )
      
      ccLine <- rbind(c(xbarcenters[idxSrc],ptSrc[p]),c(xbarcenters[idxDest],ptDest[p] ) ) 
      
      lines(ccLine,col=pcolour)
    }

  } ## Go Through Pairs Of Conditions ##
}


plotPairedChangeHistogram <- function(vIDTable,vDat,strCondTags,uLim,lLim)
{
  lChangePairs <- list()
  idx = 1
  for (gIdx in seq(1,NROW(strCondTags),2)  ) ##Iterated Through LF DF And NF Groups
  {
    gE <- strCondTags[gIdx] ##Empty Condution
    gL <- strCondTags[gIdx+1] ##With ROtifers Test Condition 
    
    datSetID <- levels(vIDTable[[gE]]$dataSetID)[ vIDTable[[gE]]$dataSetID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID] ]
    idsE <- vIDTable[[gE]]$expID[vIDTable[[gE]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID ==vIDTable[[gL]]$dataSetID]
    idsL <- vIDTable[[gL]]$expID[vIDTable[[gL]]$larvaID == vIDTable[[gL]]$larvaID && vIDTable[[gE]]$dataSetID == vIDTable[[gL]]$dataSetID]
    ptSrc <-  vDat[[gE]][levels(idsE)[idsE]]
    ptDest <- vDat[[gL]][levels(idsL)[idsL]]
    
    ##OPTIONAL: Replace NAs With 0 when Plotting So as to To Show Direction For All points
    #ptSrc <- replace(ptSrc,is.na(ptSrc),0)
    #ptDest <- replace(ptDest,is.na(ptDest),0)
    ##Otherwise Do not Plot NA connecting line
    
    ##Plot The Lines Connect Each Empty Tested Larva With Itself In THe Live Fed Conditions 
    idxSrc  <- match(gE,strCondTags) ##Bar Center Idx for Each Condition E. Fed
    idxDest <- match(gL ,strCondTags)           
    
    ##PLot Histogram Of THis Pair
    ptDelta <- ptDest-ptSrc
    ##Saturate Limits 
    ptDelta[ptDelta > uLim] = uLim
    ptDelta[ptDelta < lLim] = lLim
    
    lChangePairs[[idx]] <- ptDelta
    idx = idx+1
    hist(ptDelta,ylim = c(0,NROW(ptDelta)/3),breaks = seq(lLim,uLim,(uLim-lLim)/30),
                main =paste("                                                \t \t Change ",gE,"-",gL,sep="") )
    
  } ## Go Through Pairs Of Conditions ##
  
  names(lChangePairs) <- c("LE-LL","NE-NL","DE-DL")
  
  return (lChangePairs)
}


## Show In Which Prey Count Density Most Hunt Events Occurred - 
## Note: NOT INITIAL #Prey count but rather prey count sample
plotHuntEventPreyCountHist <- function(strCondTags,dataSetsToProcess)
{
  
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  
  for (i in strCondTags)
  {
    strDataFileName <- paste("out/setn",NROW(dataSetsToProcess),"HuntEvents",i,sep="-") ##To Which To Save After Loading
    message(paste(" Loading Hunt Events: ",strDataFileName))
    ##ExPORT 
    load(file=paste(strDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
    hist(datHuntEvent[,"PreyCount"],main=paste(i," Events #",NROW(datHuntEvent[,"PreyCount"]),sub = ""),
         breaks=seq(0,85,5),xlab="# Prey",xlim=c(0,70),ylim=c(0,300) )
  }
    
}
##par(bg="black")

colourH <- c(rgb(0.01,0.7,0.01,0.2),rgb(0.9,0.01,0.01,0.2),rgb(0.01,0.01,0.9,0.2))

################-HUNTING STAT PLOTS -###################
## Bar Plot Mean Hunting Events Per Animal #
## Common Subtitle info            ########
sampleSize = sum(unlist(datHuntStat[,"nLarva"],use.names = FALSE))
totalFrames = sum(unlist(datHuntStat[,"totalFrames"],use.names = FALSE))
strPlotName = paste("plots/HuntStat_N",sampleSize,".pdf",sep="")

## PDF OUTPUT ##
bonefile = FALSE

if (bonefile)
{
  pdf(strPlotName,paper = "a4",title="zfish Hunt Stat",width=16.27,height=22,onefile=TRUE) #col=(as.integer(filtereddatAllFrames$expID)) width=8,height=8
  par(mfrow=c(2,2)) ##MultiPlot Page
}



#huntFrames = sum(unlist(datMotionStat[,"huntframes"],use.names = FALSE))

FPS = 420;
strsub = paste("#n=", sampleSize, " #F:",totalFrames,
               "(",format(totalFrames/FPS/60,digits =3),"min)",
               " F_H/F:",format(sum(unlist(datHuntStat[,"totalHuntFrames"] ,use.names = FALSE) )/totalFrames,digits =3) ,
               " #Hunts:",sum(unlist(datHuntStat[,"groupHuntEvents"] ,use.names = FALSE) ),
               collapse=NULL)
##### Done Subtitle ##


##PLot Distribution Of Hunt Events Among Their Sampled Prey Counts - Not the initial Prey Counts##
strPlotName = "plots/HuntEventsVsPreyCount_Hist.pdf"

pdf(strPlotName,width=8,height=10,title="Hunt Vs the Prey Count they Occured under, for each Condition") #col=(as.integer(filtereddatAllFrames$expID))
#X11()
plotHuntEventPreyCountHist(strCondTags, dataSetsToProcess)
dev.off()

########### MEAN and Distribution of Prey Count At Start of Hunt EVENTS ##### 
strPlotName = "plots/meanInitialPreyCountPerCond.pdf"
vDat <- datHuntStat[,"vHInitialPreyCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]


datmean <- unlist(datHuntStat[,"initPreyCount"],use.names = FALSE)
datse <- unlist(datHuntStat[,"initsePreyCount"],use.names = FALSE)
strtitle <- "Mean Initial Prey Count per Experiment"


yl <- c(0,sampleSize/5)
xl <- c(0,max(vDat$LL,vDat$NL,vDat$DL)+10)


if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE),na.rm=TRUE)
xbarcenters <- boxplotPerCondition(vDat,datmean,datse,strtitle,strsub,strPlotName,ylim)
plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
dev.off()

##Do Initial Prey  Histogram for All Groups
strPlotName = "plots/meanInitialPreyCountPerCond_Hist.pdf"
pdf(strPlotName,width=8,height=8,title="Init Prey Per Experiment Cond Histograms") #col=(as.integer(filtereddatAllFrames$expID))

par(mfrow=c(3,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))

hist(vDat$LE,col=colourH[1],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="LE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")
hist(vDat$LL,col=colourH[1],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="LL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

hist(vDat$NE,col=colourH[2],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="NE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")
hist(vDat$NL,col=colourH[2],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="NL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

hist(vDat$DE,col=colourH[3],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="DE",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=seq(0,70,10),
     main="DL",
     xlab = "Number of Prey at start of experiment",
     ylab = "# Experiments")

dev.off()
###### # # # # # ## ## # ##  # # ## #  ## # # # 




#####  Scattter OF Hunt Events Vs Initial  Prey Counts ### ##
strPlotName = "plots/HuntEventsInFedTestVsInitPreyCount_scatter.pdf"
pdf(strPlotName,width=8,height=8,title="Number of Prey In Live Test Conditions") 


yl <- c(0,ceiling(max(datHuntStat[,"vHLarvaEventCount"]$LL,datHuntStat[,"vHLarvaEventCount"]$DL,datHuntStat[,"vHLarvaEventCount"]$NL)/10)*10 )
xl <- c(0,60)


hist(vDat$NL,col=colourH[1],ylim=yl,xlim=xl,breaks=10,
     main="Number of Prey Vs Hunting Events",
     xlab = "Number of Prey at start of experiment",
     ylab = "Number of Hunt Events")
hist(vDat$LL,col=colourH[2],ylim=yl,xlim=xl,breaks=10,add=T)
hist(vDat$DL,col=colourH[3],ylim=yl,xlim=xl,breaks=10,add=T)
box()
par(new=TRUE)
##Add Scatter Of Hunt Events
plot(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
     pch=4,col=1,
     xlab="",ylab="",ylim=yl,xlim=xl,main="")

points(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
       pch=2 ,col=2,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)

points(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
       pch=19 ,col=6,type ="p",
       xlab="",ylab="",ylim=yl,xlim=xl)

legend(50,55,legend=c("NL","LL","DL"),
       fill=colourH,
       col = c(1, 2,6),pch = c(4,2,19),
       bg = "gray90",lty = c(2, -1, 1),
       merge=TRUE)

dev.off()

### Compare FED to Non Fed Group ##
# strPlotName = "plots/NumberOfHuntEventsVsPreyCount-FedGroups.pdf"
# pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
# plot(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
#      pch=4,col=1,
#      xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,main=" Modulation of hunt rate by stimulus")
# 
# points(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
#        pch=2 ,col=2,type ="p",
#        xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl)
# 
# points(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
#        pch=19 ,col=6,type ="p",
#        xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl)
# 
# 
# legend(0.7, 55, c("NL","LL","DL"), col = c(1, 2,6),
#        text.col = "green4", lty = c(2, -1, 1), pch = c(4,2,19),
#        merge = TRUE, bg = "gray90")  
# dev.off()  
# 
########### END oF Prey Vs Hunt Events Plot ###########











#### SLICES PLot Mean Hunt Events Within Prey Cound Ranges / Slices ########
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )

checkRanges <- list(range(0,20),range(10,30),range(20,40),range(30,100))


strPlotName = paste("plots/HunteventsVsPreyCount_RangeSliced.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Event Counts Within Initial Prey Range ") 
par(mfrow=c(2,2)) ##MultiPlot Page
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
for (rng in checkRanges)
{
  print(rng)
  preyCntRange <- rng


  datSliceLL <- datHuntVsPreyLL[datHuntVsPreyLL[,1] > preyCntRange[1] & datHuntVsPreyLL[,1] < preyCntRange[2],2 ]
  datSliceNL <- datHuntVsPreyNL[datHuntVsPreyNL[,1] > preyCntRange[1] & datHuntVsPreyNL[,1] < preyCntRange[2],2 ]
  datSliceDL <- datHuntVsPreyDL[datHuntVsPreyDL[,1] > preyCntRange[1] & datHuntVsPreyDL[,1] < preyCntRange[2],2 ]

  datmean <- c(mean(datSliceLL,na.rm = TRUE),mean(datSliceNL,na.rm = TRUE),mean(datSliceDL,na.rm = TRUE))
  datse   <-  c(sd(datSliceLL,na.rm = TRUE)/sqrt(length(datSliceLL)),sd(datSliceNL,na.rm = TRUE)/sqrt(length(datSliceNL)),sd(datSliceDL,na.rm = TRUE)/sqrt(length(datSliceDL)) )

  barCenters <- barplot(height = datmean,
        names.arg = c("LL","NL","DL"),
        beside = true, las = 2,
        ylim = c(0, 35),
        cex.names = 0.75, xaxt = "n",
        main = paste("#Hunts in PreyCount Range",preyCntRange[1],preyCntRange[2]),
        sub = strsub,
        ylab = "#",
        border = "black", axes = TRUE)
  
  datlbls <- c( paste("LL \nn=",length(datSliceLL) ),
                paste("NL \nn=",length(datSliceNL) ),
                paste("DL \nn=",length(datSliceDL) ) )
  text(x = barCenters+0.2, y = 0, srt = 45,
       adj = 1.8, labels = datlbls, xpd = TRUE)
  
  
  segments(barCenters, datmean - datse * 2, barCenters,
         datmean + datse * 2, lwd = 1.5)
}
dev.off()
#################################### # # # # # # # ######################




########### MEAN Prey Count REDUCTION from Initial During Hunt EVENTS ##### 

strPlotName = "plots/meanPreyCountReductionPerLarvaHuntEpisode.pdf"
vDat <- datHuntStat[,"vHNablaPreyCount"] #INverse Sign And Denote Reduction
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]


datmean <- unlist(datHuntStat[,"meanNablaPreyCount"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seNablaPreyCount"],use.names = FALSE)
strtitle <- "Mean Prey Count Reduction from Initial at each Hunt Event Per Larva"

if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))


ylim <- max(unlist(vDat,use.names=FALSE),na.rm=TRUE)
xbarcenters <- boxplotPerCondition(vDat,datmean,datse,strtitle,strsub,strPlotName,ylim)

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
strPlotName = "plots/meanDeltaPreyCountPerLarvaHuntEpisode_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 50,lLim = -10)
title("Prey Counts Events", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()
## Plot Mean Change ######

strPlotName = "plots/meanPreyCountPerLarvaHuntEpisode_PairedChange.pdf"
meanDelta <- sapply(lDeltas,mean,na.rm =TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm =TRUE)/sqrt(sapply(lDeltas,NROW))
ylim <- 2*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
boxplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Hunt Prey Count ",strsub,strPlotName,ylim)
dev.off()

####



########### MEAN HUNTING EVENTS ##### 
#X11()
strPlotName = "plots/meanHuntingEventsPerCond.pdf"
vDat <- datHuntStat[,"vHLarvaEventCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]


datmean <- unlist(datHuntStat[,"meanHuntingEventsPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntingEventsPerLarva"],use.names = FALSE)
strtitle <- "Mean Hunting Events Per Condition"

if (!bonefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))


ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()

###### Histogram ##########
  strPlotName = "plots/meanHuntingEventsPerLarva_ChangeHist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = max(vDat$LL,vDat$DL,vDat$NL),lLim = -10)
  title("Hunt Events", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
  dev.off()
  
  ## Plot Mean Change ######
  
  strPlotName = "plots/meanHuntingEventsPairedChange.pdf"
  meanDelta <- sapply(lDeltas,mean)
  seDelta   <-  sapply(lDeltas,sd)/sqrt(sapply(lDeltas,NROW))
  ylim <- 2*max(unlist(meanDelta,use.names=FALSE))
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  #plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
  boxplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Hunt Event Count ",strsub,strPlotName,ylim)
  dev.off()
  
#### HUNT EVENTS HISTOGRAM ###
  uLim = max(vDat$LL,vDat$DL,vDat$NL)+10
  lLim = min(vDat$LL,vDat$DL,vDat$NL)
  res = 20
  strPlotName = "plots/meanHuntingEvents_Hist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE))
  
  hist(vDat$LE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="LE ",xlab="# Hunt Events")
  hist(vDat$LL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="LL ",xlab="# Hunt Events")
  hist(vDat$NE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="NE ",xlab="# Hunt Events")
  hist(vDat$NL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="NL ",xlab="# Hunt Events")
  hist(vDat$DE,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="DE ",xlab="# Hunt Events")
  hist(vDat$DL,ylim = c(0,NROW(vDat$LL)/2),breaks = seq(lLim,uLim,uLim/res), main ="DL ",xlab="# Hunt Events")
  
  dev.off()
  
  #######################################
  
#########HUnT EVENTS VS Initial Prey Count Scatter Plot 

  yl <- c(0, max(vDat$LL,vDat$DL,vDat$NL)+10)
  xl <- c(0,50)
strPlotName = "plots/NumberOfHuntEventsVsInitPreyCount.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle,onefile=TRUE) #col=(as.integer(filtereddatAllFrames$expID))
par(bg="white",fg="black")
plot(datHuntStat[,"vHInitialPreyCount"]$LE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE),
     pch=1,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$LE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="LE - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$LL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL),
     pch=2 ,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$LL ) ],type ="p",
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="LL - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$NE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE),
     pch=3,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$NE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="NE - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$NL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL),
     pch=4,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$NL ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="NF Live - Modulation of hunt rate by stimulus")

plot(datHuntStat[,"vHInitialPreyCount"]$DE,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE),
     pch=5,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$DE ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="DE  - Modulation of hunt rate by stimulus")  

plot(datHuntStat[,"vHInitialPreyCount"]$DL,as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL),
     pch=6,col=rDataset[as.numeric( datHuntStat[,"vDataSetID"]$DL ) ],
     xlab="Initial Prey Count",ylab="# hunt events",ylim=yl,xlim=xl,
     main="DL - Modulation of hunt rate by stimulus")  
dev.off()


######## EPISODE DURATION ############
#X11()

strPlotName = "plots/meanEpisodeDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanEpisodeDuration"],use.names = FALSE) #Of the Group
#vDat    <- datHuntStat[,"vmeanHLarvaDuration"] # Mean Episode Duration of Each LArva
vDat    <- datHuntStat[,"vmeanHEpisodeDurationPerLarva"]

vDatSetID <- datHuntStat[,"vDataSetID"]
datse   <- unlist(datHuntStat[,"seEpisodeDuration"],use.names = FALSE)
strtitle <- "Mean Duration of each Hunting Episode"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

lDat <- unlist(vDat,use.names=FALSE)
ylim <- max(replace(lDat,is.na(lDat),0) )
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt <- unlist(vDat[g],use.names=TRUE)
  ##OPTIONAL - PLot NA points As Zero - 
  ##vpt <-replace(vpt,is.na(vpt),0) ##NA Appear Where for exp Where no Hunting Episodes Occured
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
  strPlotName = "plots/meanEpisodeDurationOfGroup_Hist.pdf"
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 500,lLim = -610)
  title("Episode Duration", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
  dev.off()
  
  ## Plot Mean Change ######
  strPlotName = "plots/meanEpisodeDurationChange.pdf"
  meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
  seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
  ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
  #plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
  boxplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
  dev.off()
  

#######################################




######## HUNTING DURATION PER LARVA ############
#X11()
strPlotName = "plots/meanHuntDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanDuration"],use.names = FALSE)
datse   <- unlist(datHuntStat[,"seDuration"],use.names = FALSE)
vDat    <- datHuntStat[,"vHDurationPerLarva"] #Total H Duration Per Larva
vDatSetID <- datHuntStat[,"vDataSetID"]
strtitle <- "Duration of Hunting per Larva"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
# for (g in strCondTags)
# {
#   idx <- match(g,strCondTags)
#   vpt = unlist(vDat[g],use.names=FALSE)
#   #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
#   points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
# }

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)

if (!bonefile)
  dev.off()


###### Histogram ##########
strPlotName = "plots/HuntDuration_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle,onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 25000,lLim = -10000)
title("Hunt Duration", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()

strPlotName = "plots/HuntDurationChange.pdf"
meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
#plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
boxplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
dev.off()


#######################################


#######################################



############### HUNT RATIOn ###############
#X11()

strPlotName <- "plots/meanHuntRatioOfGroup.pdf"
datmean     <- unlist(datHuntStat[,"meanHuntRatioOfGroup"],use.names = FALSE)
datse       <- unlist(datHuntStat[,"seHuntRatioOfGroup"],use.names = FALSE)
vDat        <- datHuntStat[,"vLarvaHRatio"]
vDatSetID   <- datHuntStat[,"vDataSetID"]
strtitle    <- "Ratio of Time spent Hunting Over all Frames"

if (!bonefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
# for (g in strCondTags)
# {
#   idx <- match(g,strCondTags)
#   vpt = unlist(vDat[g],use.names=FALSE)
#   #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
#   points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
# }

plotConnectedPointsPairs(vIDTable,vDat,strCondTags,xbarcenters)


dev.off()


###### Histogram ##########
strPlotName = "plots/HuntRatio_Hist.pdf"
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
lDeltas <- plotPairedChangeHistogram(vIDTable,vDat,strCondTags,uLim = 0.3,lLim = -0.3)
title("Hunt Ratio", sub = strsub, cex.main = 1.2,   font.main= 1.5, col.main= "black", cex.sub = 1.0, font.sub = 2, col.sub = "black")
dev.off()

## Plot Mean Change ######
strPlotName = "plots/HuntRatioChange.pdf"
meanDelta <- sapply(lDeltas,mean, na.rm = TRUE)
seDelta   <-  sapply(lDeltas,sd,na.rm = TRUE)/sqrt(sapply(lDeltas[is.na(lDeltas)==FALSE],NROW))
ylim <- 4*max(unlist(meanDelta,use.names=FALSE))
pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))
#plot.window(xlim=c(0,100), ylim=c(-10, ylim) )
boxplotPerCondition(lDeltas,meanDelta,seDelta,"Mean Change in Episode Duration ",strsub,strPlotName,ylim)
dev.off()

#######################################


###############################################
