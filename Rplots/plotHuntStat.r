##plot Hunting Event Statistics


## Box Plots Used to Compare Conditions On Mean Stats - Saves Output As Pdf
boxplotPerCondition <- function(datStat,datMean,datse,strtitle,strsubt,stroutFileName,plotTop)
{
  
  par(mar = c(5, 6, 4, 5) + 2.5)
  datN <- unlist(datStat[,"nLarva"],use.names = FALSE)
  
  datlbls <-row.names(datStat)
  ##Add N Numbers to Labels
  datlbls <- paste(datlbls,"\nn=",datN,sep="")
  
  if(missing(plotTop)) 
  {
    plotTop <- max(datmean) +
      unique(datse[datmean== max(datmean)]) * 3
  }else
  {
    plotTop <- 1.1*plotTop ##Increase by 10%
    
  }
  
  
  barCenters <- barplot(height = datmean,
                        names.arg = datlbls,
                        beside = true, las = 2,
                        ylim = c(0, plotTop),
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
  
  segments(barCenters, datmean - datse * 2, barCenters,
           datmean + datse * 2, lwd = 1.5)
  #title(sub=strsubt)
  #  dev.off()
  
  message(strsubt)
  
  return(barCenters)
}




plotConnectedPointsPairs <- function(vIDTable,vDat,strCondTags)
{
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
    for (p in 1:NROW(ptSrc))
    {
      ccLine <- rbind(c(xbarcenters[idxSrc],ptSrc[p]),c(xbarcenters[idxDest],ptDest[p] ) ) 
      lines(ccLine,col=rDataset[as.numeric( replace( datSetID[p],is.na(datSetID[p]) ,1 ) ) ])
    }
  } ## Go Through Pairs Of Conditions ##
}
## DF GROUP ##


##par(bg="black")


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


#X11()

strPlotName = "plots/meanHuntingEventsPerLarva.pdf"
vDat <- datHuntStat[,"vHLarvaEventCount"]
vDatSetID <- datHuntStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]


datmean <- unlist(datHuntStat[,"meanHuntingEventsPerLarva"],use.names = FALSE)
datse <- unlist(datHuntStat[,"seHuntingEventsPerLarva"],use.names = FALSE)
strtitle <- "Hunting Events Per Larva"

if (!onefile)
  pdf(strPlotName,width=8,height=8,title=strtitle) #col=(as.integer(filtereddatAllFrames$expID))


ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #vDatSetID <- ,-1,vDatSetID[[g]] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags)

#dev.off() ##Close PDF





#######################################
######## EPISODE DURATION ############
#X11()

strPlotName = "plots/meanEpisodeDurationOfGroup.pdf"
datmean <- unlist(datHuntStat[,"meanEpisodeDuration"],use.names = FALSE) #Of the Group
#vDat    <- datHuntStat[,"vmeanHLarvaDuration"] # Mean Episode Duration of Each LArva
vDat    <- datHuntStat[,"vmeanHEpisodeDurationPerLarva"]

vDatSetID <- datHuntStat[,"vDataSetID"]
datse   <- unlist(datHuntStat[,"seEpisodeDuration"],use.names = FALSE)
strtitle <- "Mean Duration of each Hunting Episode"

if (!onefile)
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
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags)

if (!onefile)
  dev.off()
#######################################

######## HUNTING DURATION PER LARVA ############
#X11()
strPlotName = "plots/meanHuntDurationPerLarva.pdf"
datmean <- unlist(datHuntStat[,"meanDuration"],use.names = FALSE)
datse   <- unlist(datHuntStat[,"seDuration"],use.names = FALSE)
vDat    <- datHuntStat[,"vHDurationPerLarva"] #Total H Duration Per Larva
vDatSetID <- datHuntStat[,"vDataSetID"]
strtitle <- "Duration of Hunting per Larva"

if (!onefile)
  pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datHuntStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

plotConnectedPointsPairs(vIDTable,vDat,strCondTags)

if (!onefile)
  dev.off()
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
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt = unlist(vDat[g],use.names=FALSE)
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}


plotConnectedPointsPairs(vIDTable,vDat,strCondTags)


dev.off()
###############################################
