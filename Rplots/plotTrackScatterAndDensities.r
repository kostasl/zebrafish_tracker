rfHot <- colorRampPalette(rev(brewer.pal(11,'Spectral')));

histj<- function(x,y,x.breaks,y.breaks){
  c1 = as.numeric(cut(x,breaks=x.breaks));
  c2 = as.numeric(cut(y,breaks=y.breaks));
  mat<-matrix(0,ncol=length(y.breaks)-1,nrow=length(x.breaks)-1);
  mat[cbind(c1,c2)] = 1;
  return(mat)
}  


plotGroupMotion <- function(filtereddatAllFrames,groupStat,vexpID)
{
  yTop <- 500
  ##Note Y plotting is inverted to match video orientation ie. yTop - posY
  strDatGroup <- toString(unique(filtereddatAllFrames$group))
  
  message("PLOT Motion Tracks of each Larva noting Hunting Episodes")
  ### INDIVIDUAL TRAJECTORIES - With distinct colour for each larva ####
  for (i in vexpID)
  {
    strTrajectoryplotFileName <- paste("plots/scatter/Motion/larva/MotionTrajectories-Set-",strDatGroup,"-lID_",i,".pdf",sep="",collapse=NULL);
    message(strTrajectoryplotFileName)
    pdf(strTrajectoryplotFileName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))
    par(bg="black")
    par(fg="yellow")
    
    datLarvalAllFramesHunt <- filtereddatAllFrames[filtereddatAllFrames$expID == i & filtereddatAllFrames$LEyeAngle >=G_THRESHUNTANGLE & filtereddatAllFrames$REyeAngle <= -G_THRESHUNTANGLE,]
    datLarvalAllFramesAll <- filtereddatAllFrames[filtereddatAllFrames$expID == i,]
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
    #dev.copy(jpeg,filename=paste(strTrajectoryplotFileName,"-plot.jpg",sep=""));
    dev.off()
  }
  
  ##### Plot ALL Larvae In The group TOgether ##
  message("PLOT - Overlay All Larva of Group / noting Hunting Episodes")
  
  strTrajectoryplotFileName <- paste("plots/scatter/Motion/group/MotionTrajectories-Set-",strDatGroup,"-All.pdf",sep="",collapse=NULL);
  message(strTrajectoryplotFileName)
  pdf(strTrajectoryplotFileName,width=8,height=8) #col=(as.integer(filtereddatAllFrames$expID))
  par(bg="black")
  par(fg="yellow")
  
  colMap = colTraj(filtereddatAllFrames$expID);
  
  if (length(filtereddatAllFrames$expID) == 0)
  {
    #plot(filtereddatAllFrames$posX,yTop-filtereddatAllFrames$posY,type='p',pch='.',lwd=1,col="grey",xlim=c(80,600),ylim=c(0,500),col.axis="red")
    #plot.new()
    warning(paste("No Data To plot trajectories for :",strCond) )
    message(paste("No Data To plot trajectories for :",strCond) )
  }
  
  bFreshPlot = TRUE
  
  procMotFrames = 0;
  procHuntFrames = 0;
  hbinXY = list(); ##List Of Binarized Trajectories
  hbinHXY = list(); ##List Of Binarized Hunting Episode Trajectories
  
  ##Now PLot All Larval Tracks from the Group On the SAME PLOT ##
  idx = 0
  for (i in vexpID)
  {
    idx = idx + 1
    #message(i)
    datLarvalAllFramesHunt <- filtereddatAllFrames[filtereddatAllFrames$expID == i &
                                                     filtereddatAllFrames$REyeAngle <= -G_THRESHUNTANGLE &
                                                     filtereddatAllFrames$LEyeAngle >=G_THRESHUNTANGLE &
                                                     abs(filtereddatAllFrames$LEyeAngle-filtereddatAllFrames$REyeAngle) >= G_THRESHUNTVERGENCEANGLE,]

    procHuntFrames = procHuntFrames + NROW(datLarvalAllFramesHunt)
    
    datLarvalAllFramesAll <- filtereddatAllFrames[filtereddatAllFrames$expID == i,]
    
    procMotFrames = procMotFrames + NROW(datLarvalAllFramesAll)
    
    hbinHXY[[idx]] <- histj(datLarvalAllFramesHunt$posX,yTop-datLarvalAllFramesHunt$posY,seq(0,600,600),seq(0,yTop,20))
    hbinXY[[idx]] <- histj(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,seq(0,640,10),seq(50,yTop,10))
    
    #points(datLarvalAllFramesAll$posX,datLarvalAllFramesAll$posY,pch='.',col="white",xlim=c(80,565),ylim=c(0,500),col.axis="red")
    if (bFreshPlot)
    {
      plot(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,type='p',pch='.',col=colMap[which(vexpID == i)],xlim=c(80,600),ylim=c(0,500),col.axis="red")
      bFreshPlot = FALSE
    }else
    {
      points(datLarvalAllFramesAll$posX,yTop-datLarvalAllFramesAll$posY,pch='.',col=colMap[which(vexpID == i)],xlim=c(80,600),ylim=c(0,500),col.axis="red")
    }
    
    points(datLarvalAllFramesHunt$posX,yTop-datLarvalAllFramesHunt$posY,pch=1,lwd=2,col="red",xlim=c(80,600),ylim=c(0,500),col.axis="red")
  }##For Each Larva
  
  sampleSize  <- length(vexpID) #Number of Larvae Used 
  strtitle = paste(strCond,"Motion",collapse=NULL)
  strsub = paste("#n=", sampleSize, " #F:",groupStat$totalFrames,
                 "\n #Hunts:",groupStat$groupHuntEvents,
                 " (mu:", format(groupStat$meanHuntingEventsPerLarva,digits =3),
                 " sig:",format(groupStat$stdHuntingEventsPerLarva,digits=3),") #F_h:",groupStat$huntFrames,
                 "R_h:", format(groupStat$groupHuntRatio,digits=2),
                 "(mu:",format(groupStat$meanHuntRatioPerLarva,digits=3),"sd:",format(groupStat$stdHuntRatioPerLarva,digits=3),")" ,collapse=NULL)
  
  title(strtitle, sub = strsub, cex.main = 1.5,   font.main= 1.5, col.main= "yellow", cex.sub = 1.0, font.sub = 2, col.sub = "red")
  #dev.copy(device=jpeg,filename=paste(strTrajectoryplotFileName,"-plot.jpg"));
  dev.off()
  
  
  
  ###### BINARIZED HISTOGRAM PER GROUP ###
  ## Now Sum All LArva Binarized Trajectories and Display Heat Map
  hGroupbinDensity <- Reduce('+', hbinXY)
  strDensityplotFileName <- paste("plots/binDensity/MotionDensity-BINSet-",strCond,".pdf",collapse=NULL,sep="");
  pdf(strDensityplotFileName,width=8,height=8)
  sampleSize  <- length(vexpID) #Number of Larvae Used 
  hotMap <- c(rfHot(sampleSize),"#FF0000");
  image(seq(0,640,10),seq(50,yTop,10),hGroupbinDensity,axes=TRUE,col=hotMap,xlab="Pos X",ylab="Pos Y")
  title(paste(strCond,"Motion Trajectory Heatmap  #n=", sampleSize, " #F:",procMotFrames),collapse=NULL);
  #dev.copy(jpeg,filename=paste(strDensityplotFileName,"-plot.jpg"));
  dev.off()
  ###
  
  

  ## Now Sum All LArva Binarized Hunting Episode Trajectories and Display Heat Map
  hGroupbinDensity <- Reduce('+', hbinHXY)
  strDensityplotFileName <- paste("plots/binDensity/MotionHuntingDensity-BINSet-",strCond,".pdf",collapse=NULL,sep="");
  pdf(strDensityplotFileName,width=8,height=8)
  sampleSize  <- length(vexpID) #Number of Larvae Used 
  hotMap <- c(rfHot(sampleSize),"#FF0000");
  image(seq(0,600,600),seq(0,yTop,20),hGroupbinDensity,axes=TRUE,col=hotMap,xlab="Pos X",ylab="Pos Y")
  title(paste(strCond,"Motion Hunting Episode  Heatmap  #n=", sampleSize, " #F:",procHuntFrames),collapse=NULL);
  #dev.copy(jpeg,filename=paste(strDensityplotFileName,"-plot.jpg"));
  dev.off()
  ###
  
} ##End of Function



##Test  PlayBack Plot Hunt Event###
renderHuntEventPlayback <- function(datHuntEventMergedFrames,speed=1)
{
  X11()
  for (i in seq(1,NROW(datHuntEventMergedFrames),speed) )
  {
    tR = (1: min( c(i,NROW(datHuntEventMergedFrames) ) ) )
    posX = datHuntEventMergedFrames[max(tR),]$posX
    posY = 640-datHuntEventMergedFrames[max(tR),]$posY
    bearingRad = pi/180*(datHuntEventMergedFrames[max(tR),]$BodyAngle+90+180)
    posVX = posX+cos(bearingRad)*15
    posVY = posY-sin(bearingRad)*15
    
    plot(datHuntEventMergedFrames[tR,]$posX,640-datHuntEventMergedFrames[tR,]$posY,xlim=c(20,480),ylim=c(0,600),col="black",cex = .5,type='l',xlab="X",ylab="Y")
    points(posX,posY,col="black",pch=16)
    lines(datHuntEventMergedFrames[tR,]$Prey_X,640-datHuntEventMergedFrames[tR,]$Prey_Y,col="red")
    points(datHuntEventMergedFrames[max(tR),]$Prey_X,640-datHuntEventMergedFrames[max(tR),]$Prey_Y,col="red",pch=16)
    arrows(posX,posY,posVX,posVY)
    
    
    
   }
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