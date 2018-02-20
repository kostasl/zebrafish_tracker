
#################   MOTION ######################
##Motion Plots  - Path Length ##
par(bg="white")
sampleSize = sum(unlist(datMotionStat[,"nLarva"],use.names = FALSE))
totalFrames = sum(unlist(datMotionStat[,"totalFrames"],use.names = FALSE))
moveFrames = sum(unlist(datMotionStat[,"totalMotionFrames"],use.names = FALSE))

FPS = 420;
strsub = paste("#n=", sampleSize, " #F:",totalFrames,
               "(",format(totalFrames/FPS/60,digits =3),"min)",
               " F_M/F:",format( moveFrames/totalFrames,use.names = FALSE ,digits =3) ,
               " #Events:",sum(unlist(datMotionStat[,"totalEventCount"])),
               collapse=NULL)


strPlotName = "plots/meanPathLengthLarva.pdf"
vDat <- datMotionStat[,"vPathLengths"]
vDatSetID <- datMotionStat[,"vDataSetID"]
vIDTable <- datHuntStat[,"vIDLookupTable"]

datmean <- unlist(datMotionStat[,"meanPathLength"],use.names = FALSE)
datse <- unlist(datMotionStat[,"sePathLength"],use.names = FALSE)
strtitle <- "Mean Path Length Per Larva"

##Fish With No Hunting Events #
#datHuntStat[,"vexpID"]$NE[datHuntStat[,"vHLarvaEventCount"]$NE == 0]
##*All:unlist(datHuntStat[,"vexpID"],,use.names=FALSE)[unlist(datHuntStat[,"vHLarvaEventCount"],use.names=FALSE) == 0]
pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))

xbarcenters <- boxplotPerCondition(datMotionStat, datmean,datse,strtitle,strsub,strPlotName,ylim)
#vpch = c(0:25,32:127)
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
dev.off();


#plot(rep(xbarcenters[1],NROW(datMotionStat[[1,"vPathLengths"]]) ),datMotionStat[[1,"vPathLengths"]] )

##Motion Plots - Speed ##
strPlotName = "plots/meanSpeedLarva.pdf"
vDat <- datMotionStat[,"vSpeed"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanSpeed"],use.names = FALSE)
datse <- unlist(datMotionStat[,"seSpeed"],use.names = FALSE)
strtitle <- "Movement Speed Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt <- unlist(vDat[g],use.names=TRUE)
  ##OPTIONAL - PLot NA points As Zero - 
  ##vpt <-replace(vpt,is.na(vpt),0) ##NA Appear Where for exp Where no Hunting Episodes Occured
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

dev.off();

##Motion Ratio -  ##
strPlotName = "plots/moveRatioGroup.pdf"
vDat <- datMotionStat[,"vMovementRatio"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanMotionRatio"])
datse <- unlist(datMotionStat[,"seSpeed"],use.names = FALSE)
strtitle <- "Movement Ratio Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt <- unlist(vDat[g],use.names=TRUE)
  ##OPTIONAL - PLot NA points As Zero - 
  ##vpt <-replace(vpt,is.na(vpt),0) ##NA Appear Where for exp Where no Hunting Episodes Occured
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}


dev.off();

##Motion Sinuosity -  ##
strPlotName = "plots/pathSinuosityGroup.pdf"
vDat <- datMotionStat[,"vSinuosity"]
vDatSetID <- datMotionStat[,"vDataSetID"]
datmean <- unlist(datMotionStat[,"meanSinuosity"])
datse <- unlist(datMotionStat[,"seSinuosity"],use.names = FALSE)
strtitle <- "Path Sinuosity Ratio Per Larva"

pdf(strPlotName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))

ylim <- max(unlist(vDat,use.names=FALSE))
xbarcenters <- boxplotPerCondition(datMotionStat,datmean,datse,strtitle,strsub,strPlotName,ylim)
for (g in strCondTags)
{
  idx <- match(g,strCondTags)
  vpt <- unlist(vDat[g],use.names=TRUE)
  ##OPTIONAL - PLot NA points As Zero - 
  ##vpt <-replace(vpt,is.na(vpt),0) ##NA Appear Where for exp Where no Hunting Episodes Occured
  
  #points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=idx,col=r[idx] )
  points(rep(xbarcenters[idx],NROW(vpt) ),vpt,pch=as.numeric(vDatSetID[[g]]),col=rDataset[as.numeric(vDatSetID[[g]])] )
}

dev.off()
