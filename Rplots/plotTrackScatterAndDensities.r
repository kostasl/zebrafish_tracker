plotGroupMotion <- function(filtereddatAllFrames,groupStat,vlarvaID)
{
  yTop <- 500
  ##Note Y plotting is inverted to match video orientation ie. yTop - posY
  
  message("PLOT Motion Tracks of each Larva noting Hunting Episodes")
  ### INDIVIDUAL TRAJECTORIES - With distinct colour for each larva ####
  for (i in vlarvaID)
  {
    strTrajectoryplotFileName <- paste("plots/scatter/Motion/larva/MotionTrajectories-Set-",strCond,"-lID_",i,".pdf",collapse=NULL);
    message(strTrajectoryplotFileName)
    pdf(strTrajectoryplotFileName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))
    par(bg="black")
    par(fg="yellow")
    
    datLarvalAllFramesHunt <- filtereddatAllFrames[filtereddatAllFrames$larvaID == i & filtereddatAllFrames$LEyeAngle >=G_THRESHUNTANGLE & filtereddatAllFrames$REyeAngle <= -G_THRESHUNTANGLE,]
    datLarvalAllFramesAll <- filtereddatAllFrames[filtereddatAllFrames$larvaID == i,]
    vEvent <- unique(datLarvalAllFramesAll$fileIdx)
    #points(datLarvalAllFramesAll$posX,datLarvalAllFramesAll$posY,pch='.',col="white",xlim=c(80,565),ylim=c(0,500),col.axis="red")
    plot(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,type='p',pch='.',col="white",xlim=c(80,600),ylim=c(0,yTop),col.axis="red")
    points(datLarvalAllFramesHunt$posX,yTop-datLarvalAllFramesHunt$posY,pch='.',col="red",xlim=c(80,600),ylim=c(0,yTop),col.axis="red")
     
    for (j in vEvent)
    {
      datEvent <- datLarvalAllFramesAll[datLarvalAllFramesAll$fileIdx==j,]
      points(datEvent[1]$posX,yTop-datEvent[1]$posY,pch=12,col="red",xlim=c(80,600),ylim=c(0,500),col.axis="red")
      
    }
    
    sampleSize  <- length(unique(datLarvalAllFramesAll$fileIdx)) #Number of Larvae Used 
    strtitle = paste(strCond,"Motion",collapse=NULL)
    strsub = paste("#e=", sampleSize, " #F:",groupStat$totalFrames,collapse=NULL)
    title(strtitle, sub = strsub, cex.main = 1.5,   font.main= 1.5, col.main= "yellow", cex.sub = 1.0, font.sub = 2, col.sub = "red")
    dev.off()
  }
  
  ##### Plot ALL Larvae In The group TOgether ##
  message("PLOT - Overlay All Larva of Group / noting Hunting Episodes")
  
  strTrajectoryplotFileName <- paste("plots/scatter/Motion/group/MotionTrajectories-Set-",strCond,"-All.pdf",collapse=NULL);
  message(strTrajectoryplotFileName)
  pdf(strTrajectoryplotFileName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$larvaID))
  par(bg="black")
  par(fg="yellow")
  
  colMap = colTraj(filtereddatAllFrames$larvaID);
  
  if (length(filtereddatAllFrames$larvaID) == 0)
  {
    #plot(filtereddatAllFrames$posX,yTop-filtereddatAllFrames$posY,type='p',pch='.',lwd=1,col="grey",xlim=c(80,600),ylim=c(0,500),col.axis="red")
    #plot.new()
    warning(paste("No Data To plot trajectories for :",strCond) )
    message(paste("No Data To plot trajectories for :",strCond) )
  }
  
  bFreshPlot = TRUE
  ##Now PLot All Larval Tracks from the Group On the SAME PLOT ##
  for (i in vlarvaID)
  {
    #message(i)
    datLarvalAllFramesHunt <- filtereddatAllFrames[filtereddatAllFrames$larvaID == i &
                                                     filtereddatAllFrames$REyeAngle <= -G_THRESHUNTANGLE &
                                                     filtereddatAllFrames$LEyeAngle >=G_THRESHUNTANGLE &
                                                     abs(filtereddatAllFrames$LEyeAngle-filtereddatAllFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE,]


    datLarvalAllFramesAll <- filtereddatAllFrames[filtereddatAllFrames$larvaID == i,]
    #points(datLarvalAllFramesAll$posX,datLarvalAllFramesAll$posY,pch='.',col="white",xlim=c(80,565),ylim=c(0,500),col.axis="red")
    if (bFreshPlot)
    {
      plot(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,type='p',pch='.',col=colMap[which(vlarvaID == i)],xlim=c(80,600),ylim=c(0,500),col.axis="red")
      bFreshPlot = FALSE
    }else
    {
      points(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,pch='.',col=colMap[which(vlarvaID == i)],xlim=c(80,600),ylim=c(0,500),col.axis="red")
    }
    
    points(datLarvalAllFramesHunt$posX,yTop-datLarvalAllFramesHunt$posY,pch=1,lwd=2,col="red",xlim=c(80,600),ylim=c(0,500),col.axis="red")
  }##For Each Larva
  sampleSize  <- length(vlarvaID) #Number of Larvae Used 
  strtitle = paste(strCond,"Motion",collapse=NULL)
  strsub = paste("#n=", sampleSize, " #F:",groupStat$totalFrames,
                 "\n #Hunts:",groupStat$groupHuntEvents,
                 " (mu:", format(groupStat$meanHuntingEventsPerLarva,digits =3),
                 " sig:",format(groupStat$stdHuntingEventsPerLarva,digits=3),") #F_h:",groupStat$huntFrames,
                 "R_h:", format(groupStat$groupHuntRatio,digits=2),
                 "(mu:",format(groupStat$meanHuntRatioPerLarva,digits=3),"sd:",format(groupStat$stdHuntRatioPerLarva,digits=3),")" ,collapse=NULL)
  
  title(strtitle, sub = strsub, cex.main = 1.5,   font.main= 1.5, col.main= "yellow", cex.sub = 1.0, font.sub = 2, col.sub = "red")
  dev.off()
}

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


############# PLot Heat Map of Movement Trajectories Across COnditions #####
# strTrajectoryDensityFileName <- paste("plots/densities/MotionDensity-Set-",strCond,".pdf",collapse=NULL);
# pdf(strTrajectoryDensityFileName,width=8,height=8)
# eGroupDens <- kde2d(filtereddatAllFrames$posX,filtereddatAllFrames$posY, n=60, lims=c(range(0,565),range(0,565)) )
# image(eGroupDens,col=r)
# #title(paste(strCond,"Group Motion Densities #n=", sampleSize, " #F:",procDatFrames),collapse=NULL);
# title(strtitle, sub = strsub, cex.main = 1.5,   font.main= 1.5, cex.sub = 1.0, font.sub = 2)
# dev.off()
#############################
########