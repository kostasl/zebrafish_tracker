

rfHot <- colorRampPalette(rev(brewer.pal(11,'Spectral')));

histj<- function(x,y,x.breaks,y.breaks){
  c1 = as.numeric(cut(x,breaks=x.breaks));
  c2 = as.numeric(cut(y,breaks=y.breaks));
  mat<-matrix(0,ncol=length(y.breaks)-1,nrow=length(x.breaks)-1);
  mat[cbind(c1,c2)] = 1;
  return(mat)
}  

## Plot Mean And Std Error Around Mean ##
plotMeanEyeV <- function(lEyeVDistMatrix,lcolour,addNewPlot=TRUE)
{
  
  lEyeVMatrix <-lEyeVDistMatrix
  nSPerX    <- apply(lEyeVMatrix,2,function(x){return (NROW(x[!is.na(x)])) })
  
  bandUpper <- apply(lEyeVMatrix,2,mean,na.rm=TRUE ) + apply(lEyeVMatrix,2,sd,na.rm=TRUE)/sqrt(nSPerX)
  bandLower <- apply(lEyeVMatrix,2,mean,na.rm=TRUE ) - apply(lEyeVMatrix,2,sd,na.rm=TRUE)/sqrt(nSPerX)
  
  ##plot Only Where we Have more than 1 sample
  x<- seq(0,maxDist,stepDist)[nSPerX > 1]
  if (addNewPlot)
    
    plot(x,smooth( apply(lEyeVMatrix,2,mean,na.rm=TRUE)[nSPerX > 1],kind="3R" ),
        type="l",col=lcolour,lwd=3,ylim=c(0,80),xlim=c(0,5),
        xlab=NA,ylab=NA,cex.lab = FONTSZ_AXISLAB,cex.axis=FONTSZ_AXIS)
  else
    lines(x,smooth(apply(lEyeVMatrix,2,mean,na.rm=TRUE)[nSPerX > 1],kind="3R" ),
         type="l",col=lcolour,lwd=3,ylim=c(40,80),xlim=c(0,5),
         xlab=NA,ylab=NA)
  
  polygon(c(x, rev(x )),
          c(smooth(bandUpper[nSPerX > 1]) ,
            smooth(rev(bandLower[nSPerX > 1])) ),
          col=lcolour,
          lwd=1,ylim=c(40,80),xlab=NA,ylab=NA)
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

bPause <- FALSE
keydown <- function(key) {
  if (key == "p") 
  {
    bPause <- TRUE
    message("KEYPRESS ")
    return( bPause)
  }
    
  NULL
}




##Test  PlayBack Plot Hunt Event###
renderHuntEventPlayback <- function(datHuntEventMergedFrames,preyTargetID,speed=1,saveToFolder=NA)
{

  #datHuntEventMergedFrames$LEyeAngle <- meanf(datHuntEventMergedFrames$LEyeAngle,20)
  #datHuntEventMergedFrames$REyeAngle <- meanf(datHuntEventMergedFrames$REyeAngle,20)
  
  
  frameWidth = 610
  frameHeight = 610
  
  X_FRAME <- c(0,frameWidth)
  Y_FRAME <- c(0,frameHeight)
  
  iConeLength = 100
  ## (see Bianco et al. 2011) : "the functional retinal field as 163˚ after Easter and Nicola (1996)."
  iConeArc = 163/2 ##Degrees Of Assumed Half FOV of Each Eye
  ##Eye Distance taken By Bianco As 453mum, ie 0.5mm , take tracker
  EyeDist         = 0.45/DIM_MMPERPX ##From Head Centre
  BodyArrowLength = DIM_DISTTOMOUTH_PX
  LEyecolour      = "#0C0CFF2A"
  REyecolour      = "#FF00002A"

  #display.brewer.all() to see avaulable options
  Polarrfc <- colorRampPalette(rev(brewer.pal(8,'Dark2')));
  
  
  datHuntEventMergedFrames <- datHuntEventMergedFrames[datHuntEventMergedFrames$posX < frameWidth &
                                                         datHuntEventMergedFrames$posY < frameHeight & 
                                                         !is.na(datHuntEventMergedFrames$frameN) ,]
 
  X11() ##Show On Screen
  setGraphicsEventHandlers(prompt = "Click p to pause",
                           onMouseDown = NULL,
                           onMouseUp = NULL,
                           onIdle=NULL,
                           onKeybd = keydown,
                           consolePrompt="Press Key To Pause")
  
  eventEnv <- getGraphicsEventEnv()
  
  
  startFrame <- min(datHuntEventMergedFrames$frameN,na.rm =TRUE)
  endFrame   <- max(datHuntEventMergedFrames$frameN,na.rm =TRUE)
  vPreyFrameN   <- datHuntEventMergedFrames[datHuntEventMergedFrames$PreyID == preyTargetID ,]$frameN
  lastPreyFrame <- max(vPreyFrameN[!is.na(vPreyFrameN)])
  
  for (i in seq(startFrame,endFrame,speed) )
  {
    while (bPause) 
    {
      key<- readline(prompt="- Press r to continue -")
      if (key == 'r')
        bPause <- FALSE
    }
    
    
    
    tR = (startFrame: min( c(i,endFrame ) ) )
    ##Multiple Copies Of Fish Can Exist As its Joined the Food Records, when tracking more than one Food Item.
    ## Thus When Rendering the fish Choose one of the food items that appears in the current frame range
    
    
    datFishFrames <- datHuntEventMergedFrames[datHuntEventMergedFrames$frameN %in% tR,] ##in Range
    vTrackedPreyIDs <- unique(datFishFrames$PreyID)
    
    lastFrame <- i
    ##If this specific Frame Does not Exist In the Dat, Then Take The Last One within Range
    if (NROW(datFishFrames[datFishFrames$frameN == i,]) < 1)
      lastFrame <- max(datFishFrames$frameN)
     
    
    ##Filter The Fish Motion In the Subset PreyID Selection
    #datFishFrames <- filterEyeTailNoise(datFishFrames)
    recLastFishFrame <- datFishFrames[datFishFrames$frameN == lastFrame,]
    
    ##There Could Be Multiple With Thaty Frame N - Isolate Single Record ##
    if (NROW(recLastFishFrame) > 1)
    {
      if (is.na(preyTargetID) )
        #if (NROW(datFishFrames[datFishFrames$frameN == lastFrame & !is.na(datFishFrames$PreyID),]))
        preyTargetID <- min(c(datFishFrames[datFishFrames$frameN == lastFrame,]$PreyID ) ) ##Choose A Prey ID found on the Last Frame The max Id F
      
      
        recLastFishFrame <- datFishFrames[datFishFrames$PreyID == preyTargetID ,]
    }
    ##Now Isolate Fish Rec, Focus on Single Prey Item
    


    
    
    posX = recLastFishFrame$posX
    posY = frameWidth-recLastFishFrame$posY
    bearingRad = pi/180*(recLastFishFrame$BodyAngle-90)##+90+180 - Body Heading
    TailRad <- vector()
    TailRad[1] =   pi/180*(recLastFishFrame$DThetaSpine_1 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[2] =   pi/180*(recLastFishFrame$DThetaSpine_2 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[3] =   pi/180*(recLastFishFrame$DThetaSpine_3 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[4] =   pi/180*(recLastFishFrame$DThetaSpine_4 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[5] =   pi/180*(recLastFishFrame$DThetaSpine_5 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[6] =   pi/180*(recLastFishFrame$DThetaSpine_6 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    TailRad[7] =   pi/180*(recLastFishFrame$DThetaSpine_7 + recLastFishFrame$ThetaSpine_0 - 90) #Tail - bearingRad+pi
    posVX = posX+cos(bearingRad)*BodyArrowLength
    posVY = posY-sin(bearingRad)*BodyArrowLength
    
    dev.hold()
    ##Plot Track
    par(bg="white") 
    plot(datFishFrames$posX,frameWidth-datFishFrames$posY,xlim=X_FRAME,ylim=Y_FRAME,col="black",cex = .5,type='l',xlab="X",ylab="Y")
    
    
    ##Plot Current Frame Position
    points(posX,posY,col="black",pch=16)
    arrows(posX,posY,posVX,posVY)
    
    ##Draw Heading Line Of Sight In Blue
    posVX2 = posX+cos(bearingRad)*BodyArrowLength*10
    posVY2 = posY-sin(bearingRad)*BodyArrowLength*10
    arrows(posX,posY,posVX2,posVY2,length=0.01,col="blue") ##Draw Heading Forward Arrow
    
    ##Draw Tail Segment Motion
    posTX2 = posX
    posTY2 = posY
    for (s in 1:NROW(TailRad))
    {
      posTX2n <- posTX2+cos(TailRad[s])*BodyArrowLength*1
      posTY2n <- posTY2-sin(TailRad[s])*BodyArrowLength*1
      arrows(posTX2,posTY2,posTX2n,posTY2n,length=0.03,col="magenta") ##Draw Heading Forward Arrow
      posTX2 <- posTX2n ##Next Tail Segment Is Drawn from the end of previous one
      posTY2 <- posTY2n
    }
    
        
    
    ##Draw Eyes 
    ##Left Eye - Requires Inversions due to differences in How Angles Are Calculated in Tracker and In R Plots
    LEyePosX <- posX-cos(bearingRad+pi/180*(45+90))*EyeDist
    LEyePosY <- posY+sin(bearingRad+pi/180*(45+90))*EyeDist
    
    #LEyeConeX <- c(LEyePosX,LEyePosX-cos(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90-iConeArc))*iConeLength,                   LEyePosX-cos(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90+iConeArc))*iConeLength )
    #LEyeConeY <- c(LEyePosY,LEyePosY+sin(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90-iConeArc))*iConeLength,                   LEyePosY+sin(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90+iConeArc))*iConeLength )
    nsteps = 10
    rs <- seq(-iConeArc,+iConeArc,len=nsteps) 
    LEyeConeX <- c(LEyePosX,LEyePosX-cos(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90-rs))*iConeLength)
    LEyeConeY <- c(LEyePosY, LEyePosY+sin(bearingRad+pi/180*(recLastFishFrame$LEyeAngle+90-rs))*iConeLength)
    
    polygon(LEyeConeX,LEyeConeY,col=LEyecolour) #density=20,angle=45

    ##Right Eye
    REyePosX <- posX-cos(bearingRad+pi/180*(-45-90))*EyeDist
    REyePosY <- posY+sin(bearingRad+pi/180*(-45-90))*EyeDist
    
    REyeConeX <- c(REyePosX,
                   REyePosX-cos(bearingRad+pi/180*(recLastFishFrame$REyeAngle-90-rs))*iConeLength )
    
    REyeConeY <- c(REyePosY,
                   REyePosY+sin(bearingRad+pi/180*(recLastFishFrame$REyeAngle-90-rs))*iConeLength)
    
    polygon(REyeConeX,REyeConeY,col=REyecolour) ##,density=25,angle=-45
    
    
    ##Draw Frame Number
    #text(X_FRAME[1] + 60,frameHeight+20,labels=paste(i,"# (",i-startFrame,")",round( (i-startFrame)/(Fs/1000)),"msec" ) ,col="darkblue",cex=0.7)
    mtext(side = 3,cex=1.0, line = 2.2, outer=FALSE, 
          paste(i,"# (",i-startFrame,")",round( (i-startFrame)/(Fs/1000)),"msec" ) ,col="darkblue" )
    
    colR <- c(Polarrfc(NROW(vTrackedPreyIDs) ) ,"#FF0000");
    ###Draw Prey
    nn <- 0
    for (f in vTrackedPreyIDs)
    {
      nn <- nn + 1
      lastPreyFrame <- datHuntEventMergedFrames[datHuntEventMergedFrames$frameN == lastFrame & datHuntEventMergedFrames$PreyID == f,]
      rangePreyFrame <- datHuntEventMergedFrames[datHuntEventMergedFrames$frameN >= startFrame & datHuntEventMergedFrames$frameN <= lastFrame & datHuntEventMergedFrames$PreyID == f,]
      
      if (NROW(lastPreyFrame$Prey_X) > 0 )
      {
        
        points(lastPreyFrame$Prey_X,frameWidth-lastPreyFrame$Prey_Y,col=colR[[nn]],pch=16,cex=lastPreyFrame$Prey_Radius/2)
        if (lastPreyFrame$Prey_Radius < 2) ##Draw X over Prey, If it has likely Dissappeared By Now   
          points(lastPreyFrame$Prey_X,frameWidth-lastPreyFrame$Prey_Y,col=colR[[nn]],pch=4,cex=1.2)
        
        lines(rangePreyFrame$Prey_X,frameWidth-rangePreyFrame$Prey_Y,col="red")
        text(lastPreyFrame$Prey_X+5,frameWidth-lastPreyFrame$Prey_Y+10,labels=f,col="darkred",cex=0.8)
      }
    }
    
    dev.flush()
    if (!is.na(saveToFolder) )
    {
      dev.copy(jpeg,filename=paste(saveToFolder,"/",sprintf("%05d", i) ,".jpg",sep=""), bg="white" ,quality=80);
      dev.off ();
    }
    
   } ##For Each Frame
  
  
}##RenderHunt Event








## PLot The Relative Angle Of Fish Bearing to Prey Over Time On  a Polar Plot - For Each Prey Of this Hunt Event
## \returns the relative Angle Of Each Prey To The Fish;s Heading
polarPlotAngleToPrey <- function(datRenderHuntEvent)
{
### Plot Relative Angle To Each Prey ###
vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)

Range <- ((max(datRenderHuntEvent[!is.na(datRenderHuntEvent$PreyID),]$frameN) 
           - min(datRenderHuntEvent[!is.na(datRenderHuntEvent$PreyID),]$frameN) ) / G_APPROXFPS)+1

relAngle <- list()


txtW <- -0.2# strwidth(parse(text=paste("270", "^o ", sep="")))
plot(1,type='n',xlim=c(-(Range+txtW),(Range+txtW)) ,ylim=c(-(Range+txtW),(Range+txtW) ),main="Angle to Prey Over Time ")


#display.brewer.all() to see avaulable options
Polarrfc <- colorRampPalette(rev(brewer.pal(8,'Dark2')));
  colR <- c(Polarrfc(NROW(vTrackedPreyIDs) ) ,"#FF0000");

  n <- 0
  for (f in vTrackedPreyIDs)
  {
    n<-n+1
    message(f)
    message(colR[n])
    if (is.na(f))
      next
    datRenderPrey <- datRenderHuntEvent[datRenderHuntEvent$PreyID == f,]
    ##Atan2 returns -180 to 180, so 1st add 180 to convert to 360, then sub the fishBody Angle, then Mod 360 to wrap in 360deg circle, then sub 180 to convert to -180 to 180 relative to fish heading angles
    #relAngle[[as.character(f)]] <- ( ((360+180/pi * atan2( datRenderHuntEvent$Prey_X-datRenderHuntEvent$posX,datRenderHuntEvent$posY - datRenderHuntEvent$Prey_Y)) - datRenderHuntEvent$BodyAngle) %% 360) -180
    
    relAngle[[as.character(f)]]  <- (  ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX,datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
    
    #points(relAngle[[as.character(f)]],datRenderPrey$frameN,type='b',cex=0.2,xlim=c(-180,180))
    
    ##Convert Frames To Seconds
    
    d <- (datRenderPrey$frameN-min(datRenderHuntEvent[!is.na(datRenderHuntEvent$PreyID),]$frameN)) / G_APPROXFPS
    x <- (d)*cos(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    y <- (d)*sin(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    points(x,y,type='p',cex=0.2,xlim=c(-(Range),(Range) ) ,ylim=c(-(Range),(Range) ), main="",col=colR[n])
    
    points(0,0,cex=0.8,col="blue")
    for (i in seq(0,Range,0.5 )  )
    {
      lines(i*cos(pi/180 * seq(0,360,1) ),i*sin(pi/180 * seq(0,360,1) ),col="blue")
      txtW <- strwidth(paste(as.character(i),"s",sep="") )/2
      text(i*cos(pi/180 * 0 )+txtW,i*sin(pi/180 * 0 ),labels = paste(as.character(i),"s",sep="")   ,col="blue",cex=0.7)
    }
    
    lines(c(0,0),c(0,Range+Range/30) ,col="blue") 
    txtW <- strwidth(parse(text=paste("270", "^o ", sep="")))/2
    text((Range+txtW)*cos(pi/180 * seq(0,-270,-90) + pi/2)+Range/40,(Range+txtW)*sin(pi/180 *seq(0,-270,-90) + pi/2) ,labels = parse(text=paste(seq(0,270,90), "^o ", sep="")) ,col="blue",cex=0.8)
  }
  
return (relAngle)
}


## Returns A list of vectors showing bearing Angle To Each Prey 
calcRelativeAngleToPrey <- function(datRenderHuntEvent)
{
  ### Plot Relative Angle To Each Prey ###
  vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)
  
  Range <- ((max(datRenderHuntEvent[!is.na(datRenderHuntEvent$PreyID),]$frameN) - min(datRenderHuntEvent$frameN) ) / G_APPROXFPS)+1
  relAngle <- list()
  
  n <- 0
  for (f in vTrackedPreyIDs)
  {
    n<-n+1
    #message(f)
  
    if (is.na(f))
      next
    
    datRenderPrey <- datRenderHuntEvent[datRenderHuntEvent$PreyID == f,]
    ##Atan2 returns -180 to 180, so 1st add 180 to convert to 360, then sub the fishBody Angle, then Mod 360 to wrap in 360deg circle, then sub 180 to convert to -180 to 180 relative to fish heading angles
    ##dd Time Base As frame Number on First Column
    relAngle[[as.character(f)]]  <- cbind(datRenderPrey$frameN, 
                                      ( ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX,datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
                                    )
  }
    #points(relAngle[[as.character(f)]],datRenderPrey$frameN,type='b',cex=0.2,xlim=c(-180,180))
    
    ##Convert Frames To Seconds
  return (relAngle)
}
# ### DUBlicate



### Calc Relative Angle  To Prey / Azimuth - Returns Vector
calcPreyAzimuth <- function(datRenderHuntEvent)
{
  
  ### Plot Relative Angle To Each Prey ###
  vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)
  relAngle <- list()
  
  n <- 0
  for (f in vTrackedPreyIDs)
  {
    n<-n+1
    message(f)
    if (is.na(f))
      next
    
    datRenderPrey <- datRenderHuntEvent[datRenderHuntEvent$PreyID == f,]
    ##Atan2 returns -180 to 180, so 1st add 180 to convert to 360, then sub the fishBody Angle, then Mod 360 to wrap in 360deg circle, then sub 180 to convert to -180 to 180 relative to fish heading angles
    #relAngle[[as.character(f)]] <- ( ((360+180/pi * atan2( datRenderHuntEvent$Prey_X-datRenderHuntEvent$posX,datRenderHuntEvent$posY - datRenderHuntEvent$Prey_Y)) - datRenderHuntEvent$BodyAngle) %% 360) -180
    #points(relAngle[[as.character(f)]],datRenderPrey$frameN,type='b',cex=0.2,xlim=c(-180,180))
    
    ##Convert Frames To Seconds
    bearingRad = pi/180*(datRenderPrey$BodyAngle-90)##+90+180 - Body Heading
    posVX = datRenderPrey$posX -cos(bearingRad)*DIM_DISTTOMOUTH_PX
    posVY = datRenderPrey$posY+sin(bearingRad)*DIM_DISTTOMOUTH_PX
    ##For Rel Angle Use Bladder Centroid So As to minimize angle error
    
    ##For Distance Use Estimated MouthPOint
    d <- sqrt(  (datRenderPrey$Prey_X -posVX )^2 + (datRenderPrey$Prey_Y - posVY)^2   ) 
    relAngle[[n]]  <- cbind(preyID=f,distPX=d,
                                          azimuth=( ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX, datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
    )
    
        x <- (d)*cos(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    y <- (d)*sin(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    
  }
  
  return(relAngle)
}



############# A Linear 2 Axis Plot - Prey Vs Distance
plotAngleToPreyAndDistance <- function(datRenderHuntEvent,vDistToPrey_Fixed_FullRange,t)
{
  
  ##  Angle To Prey ##
  par(new = FALSE)
  par(mar=c(4,4,4,4))
  plot(t,vDistToPrey_Fixed_FullRange[1:NROW(t)]*DIM_MMPERPX,type='l',
       xlab="(msec)",
       ylab=NA,
       col="purple",
       main="Motion Relative Prey and Eye Angles",
       asp=1,
       lwd=2,ylim=c(0,5))
  axis(side = 2,col="purple",cex=1.2,lwd=2)
  Polarrfc <- colorRampPalette(rev(brewer.pal(8,'Dark2')));
  colR <- c(Polarrfc(NROW(tblPreyRecord) ) ,"#FF0000");
  n<-0
  ##Add Prey Angle On Separate Axis
  par(new=TRUE) ##Add To Path Length Plot But On Separate Axis So it Scales Nicely
  
  for (vAToPrey in lAngleToPrey)
  {
    l <- min(NROW(t),NROW(vAToPrey))
    n<-n+1; 
    plot((vAToPrey[1:l,1]-min(datRenderHuntEvent$frameN))/(Fs/1000),filtfilt(bf_eyes,vAToPrey[1:l,2]),type='l',axes=F,col=colR[n]
         ,xlab=NA,ylab=NA, ylim=c(-40,40))
  }
  axis(side = 4,col=colR[n])
  mtext(side = 4,cex=0.8, line = 2.2, expression('Angle To Prey'^degree), font=2 )
  mtext(side = 2,cex=0.8, line = 2.2, expression("Distance To Prey (mm)"), font=2 ) 
  
  legend("bottomleft",c(paste("Distance to Prey "),paste("Angle to Prey",names(vAToPrey)) ) , #,selectedPreyID
         col=c("purple",colR),cex=0.7,box.lwd =0,lty=1,lwd=2 )
  ###
} ## Plot Prey Angle And Distance 


## PLot The Relative Angle Of Fish Bearing to Prey Over Distance to Prey as a Polar Plot
##- Can Deal With Multiple Prey IDS, 
## The Is assumed to be at the centre of the polar plot 
## Colour Code According To Eye Vergence

## \Returns prey azimuth Vector
polarPlotAngleToPreyVsDistance <- function(datRenderHuntEvent,newPlot=TRUE)
{
  Range <- 80 ##300 Pixels Around the prey
  
  ### Plot Relative Angle To Each Prey ###
  vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)
    
  #display.brewer.all() to see avaulable options
  ##Choose Heat Map For white being Low (BG) Red High Vergence
  Polarrfc <-  colorRampPalette((brewer.pal(9,'YlOrRd' ))); ##Color Bling Friendly Pallet
  colR <- (c(Polarrfc(100 ))); ##Assume 80 Degrees Max EyeVergence
  #colR["alpha",] <- 110 ##Opacity
  
  relAngle <- list()
  
  #txtW <- strwidth(parse(text=paste("270", "^o ", sep=""))) ##Override as it fails When In Layout Mode
  txtW <- -0.1# strwidth(parse(text=paste("270", "^o ", sep="")))
  fgColor <- "white"
  if (newPlot)
  {
    plot(1,type='n',xlim=c(-(Range+4*txtW),(Range+4*txtW)) ,
         ylim=c(-(Range+4*txtW),(Range+4*txtW) ),
         main="Angle to Prey Vs Distance ",
         xlab=NA,ylab=NA)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(0,0,0.3,0.99))
    
    ## Make Range Circle Llines
    lines(c(0,0),c(0,Range*0.85) ,col=fgColor,lty=2,lwd=1) #V Line To 0
    txtW <- strwidth(parse(text=paste("270", "^o ", sep="")))/3
    text((Range+txtW/2)*cos(pi/180 * seq(0,-270,-90) + pi/2),
         (Range+txtW/2)*sin(pi/180 *seq(0,-270,-90) + pi/2),
         labels = parse(text=paste(seq(0,270,90), "^o ", sep="")) ,col=fgColor,cex=0.8,font=1.5)
    
    points(0,0,cex=0.8,col="blue")
    for (i in seq(0,Range,1/DIM_MMPERPX )  )
    {
      lines(i*cos(pi/180 * seq(0,360,1) ),i*sin(pi/180 * seq(0,360,1) ),col=fgColor)
      txtW <- strwidth(paste(as.character(i*DIM_MMPERPX),"",sep="") )/2
      txtH <- strheight(paste(as.character(i*DIM_MMPERPX),"",sep="") )/2
      ## Place the Distance Labels
      text(i*cos(pi/180 * -90 ),i*sin(pi/180 * -90 )-txtH,labels = paste(as.character(i*DIM_MMPERPX),"mm",sep="") ,
           col=fgColor,cex=0.7,font=1.8)
    }
    
    ##Plot Heat Map Legend
    x <- Range/2+(1:Range/2) ##Make Narrow 1/2 length bar
    points(x,rep(-70,NROW(x) ),pch=19,col=colR,cex=1.5)
    text(x[1],-78,labels=expression("0"^degree),col=fgColor,font=2.2)  ##0 V Angle
    txtW <- strwidth(parse(text=c(expression(),bquote( .(G_THRESHUNTVERGENCEANGLE)^degree))))/2 ##Text Width For Centre Aligment 
    segments(x[1]+ G_THRESHUNTVERGENCEANGLE/2-txtW,-72,x[1]+ G_THRESHUNTVERGENCEANGLE/2-txtW,-68) ##V Indicator Of Hunting Threshold
    text(x[1]+ G_THRESHUNTVERGENCEANGLE/2-txtW,-78,labels=c(expression(),bquote( .(G_THRESHUNTVERGENCEANGLE)^degree)),col=fgColor,font=2.0)  ##0 V Angle
    
    text(tail(x,1)-txtW,-78,labels=expression("80"^degree),col=fgColor,font=2.2)  ##0 V Angle
  }
  
  
  
  n <- 0
  for (f in vTrackedPreyIDs)
  {
    n<-n+1
    message(f)
    message(colR[n])
    if (is.na(f))
      next
    
    datRenderPrey <- datRenderHuntEvent[datRenderHuntEvent$PreyID == f,]
    ##Atan2 returns -180 to 180, so 1st add 180 to convert to 360, then sub the fishBody Angle, then Mod 360 to wrap in 360deg circle, then sub 180 to convert to -180 to 180 relative to fish heading angles
    #relAngle[[as.character(f)]] <- ( ((360+180/pi * atan2( datRenderHuntEvent$Prey_X-datRenderHuntEvent$posX,datRenderHuntEvent$posY - datRenderHuntEvent$Prey_Y)) - datRenderHuntEvent$BodyAngle) %% 360) -180
    
    EyeVergence <- datRenderPrey$LEyeAngle-datRenderPrey$REyeAngle
    EyeVergence[EyeVergence<1] <- 1 ##Fix Min Value To Belong to 20 Degrees V
    #points(relAngle[[as.character(f)]],datRenderPrey$frameN,type='b',cex=0.2,xlim=c(-180,180))
    
    ##Convert Frames To Seconds
    bearingRad = pi/180*(datRenderPrey$BodyAngle-90)##+90+180 - Body Heading
    posVX = datRenderPrey$posX -cos(bearingRad)*DIM_DISTTOMOUTH_PX
    posVY = datRenderPrey$posY+sin(bearingRad)*DIM_DISTTOMOUTH_PX
    ##For Rel Angle Use Bladder Centroid So As to minimize angle error
    relAngle[[as.character(f)]]  <- ( ( 180 +  180/pi * atan2(datRenderPrey$Prey_X -datRenderPrey$posX, datRenderPrey$posY - datRenderPrey$Prey_Y)) -datRenderPrey$BodyAngle    ) %% 360 - 180
  
    ##For Distance Use Estimated MouthPOint
    d <- sqrt(  (datRenderPrey$Prey_X -posVX )^2 + (datRenderPrey$Prey_Y - posVY)^2   ) 
    x <- (d)*cos(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    y <- (d)*sin(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    points(x,y,type='p',cex=0.2,xlim=c(-(Range),(Range) ) ,ylim=c(-(Range),(Range) ), main="",
           col=colR[EyeVergence]) ##Color Accourding To EyeVergence
    
  
    
  }
  
  return(relAngle)
} ## End of PlotAngleToPreyVsDistance 

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