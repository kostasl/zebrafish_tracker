

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
renderHuntEventPlayback <- function(datHuntEventMergedFrames,speed=1,saveToFolder=NA)
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
  EyeDist = 0.4/DIM_MMPERPX ##From Head Centre
  BodyArrowLength = 13
  LEyecolour = "#0C0CFF8A"
  REyecolour = "#FF00008A"

  #display.brewer.all() to see avaulable options
  Polarrfc <- colorRampPalette(rev(brewer.pal(8,'Dark2')));
  
  
  datHuntEventMergedFrames <- datHuntEventMergedFrames[datHuntEventMergedFrames$posX < frameWidth & datRenderHuntEvent$posY < frameHeight & !is.na(datHuntEventMergedFrames$frameN) ,]
 
  X11() ##Show On Screen
  startFrame <- min(datHuntEventMergedFrames$frameN)
  endFrame <- max(datHuntEventMergedFrames$frameN)

  for (i in seq(startFrame,endFrame,speed) )
  {
    
    tR = (startFrame: min( c(i,endFrame ) ) )
    ##Multiple Copies Of Fish Can Exist As its Joined the Food Records, when tracking more than one Food Item.
    ## Thus When Rendering the fish Choose one of the food items that appears in the current frame range
    
    
    datFishFrames <- datHuntEventMergedFrames[datHuntEventMergedFrames$frameN %in% tR,] ##in Range
    vTrackedPreyIDs <- unique(datFishFrames$PreyID)
    
    lastFrame <- i
    if (NROW(datFishFrames[datFishFrames$frameN == i,]) < 1)
      lastFrame <- max(datFishFrames$frameN)
     
    
    preyTargetID <- min(c(datFishFrames[datFishFrames$frameN == lastFrame,]$PreyID ) ) ##Choose A Prey ID found on the Last Frame The max Id F
    ##Now Isolate Fish Rec, Focus on Single Prey Item
    if (!is.na(preyTargetID))
      datFishFrames <- datFishFrames[datFishFrames$PreyID == preyTargetID ,]
    
    recLastFishFrame <- datFishFrames[datFishFrames$frameN == lastFrame,]
    
    
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
    text(X_FRAME[1] + 30,frameHeight-10,labels=paste("#",i,round(1000* (i-startFrame)/(Fs)),"msec" ) ,col="darkblue",cex=0.7)
    
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
        
        points(lastPreyFrame$Prey_X,frameWidth-lastPreyFrame$Prey_Y,col=colR[[nn]],pch=16,cex=lastPreyFrame$Prey_Radius/5)
        lines(rangePreyFrame$Prey_X,frameWidth-rangePreyFrame$Prey_Y,col="red")
        text(lastPreyFrame$Prey_X+5,frameWidth-lastPreyFrame$Prey_Y+10,labels=f,col="darkred",cex=0.5)
      }
    }
    
    dev.flush()
    if (!is.na(saveToFolder) )
    {
      dev.copy(png,filename=paste(saveToFolder,"/",sprintf("%05d", i) ,".png",sep=""), bg="white" );
      dev.off ();
    }
    
   } ##For Each Frame
  
  
}##RenderHunt Event


## PLot The Relative Angle Of Fish Bearing to Prey Over Time On  a Polar Plot - For Each Prey Of this Hunt Event
polarPlotAngleToPrey <- function(datRenderHuntEvent)
{
### Plot Relative Angle To Each Prey ###
vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)
Range <- ((max(datRenderHuntEvent$frameN) - min(datRenderHuntEvent$frameN) ) / G_APPROXFPS)+1
relAngle <- list()


txtW <- strwidth(parse(text=paste("270", "^o ", sep="")))
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
    
    d <- (datRenderPrey$frameN-min(datRenderHuntEvent$frameN)) / G_APPROXFPS
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
  
}

#

## PLot The Relative Angle Of Fish Bearing to Prey Over Distance to Prey as a Polar Plot
##- For Each Prey Of this Hunt Event
## The Is assumed to be at the centre of the polar plot 
polarPlotAngleToPreyVsDistance <- function(datRenderHuntEvent)
{
  ### Plot Relative Angle To Each Prey ###
  vTrackedPreyIDs <- unique(datRenderHuntEvent$PreyID)
  Range <- 80 ##300 Pixels Around the prey
  relAngle <- list()
  
  txtW <- strwidth(parse(text=paste("270", "^o ", sep="")))
  plot(1,type='n',xlim=c(-(Range+4*txtW),(Range+4*txtW)) ,ylim=c(-(Range+4*txtW),(Range+4*txtW) ),main="Angle to Prey Vs Distance ")
  
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
    
    d <- sqrt(  (datRenderPrey$Prey_X -datRenderPrey$posX )^2 + (datRenderPrey$Prey_Y -datRenderPrey$posY)^2   ) 
    x <- (d)*cos(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    y <- (d)*sin(2*pi-pi/180 * relAngle[[as.character(f)]] + pi/2)
    points(x,y,type='p',cex=0.2,xlim=c(-(Range),(Range) ) ,ylim=c(-(Range),(Range) ), main="",col=colR[n])
    
    points(0,0,cex=0.8,col="blue")
    for (i in seq(0,Range,1/DIM_MMPERPX )  )
    {
      lines(i*cos(pi/180 * seq(0,360,1) ),i*sin(pi/180 * seq(0,360,1) ),col="blue")
      txtW <- strwidth(paste(as.character(i*DIM_MMPERPX),"",sep="") )/2
      text(i*cos(pi/180 * 0 )+txtW,i*sin(pi/180 * 0 ),labels = paste(as.character(i*DIM_MMPERPX),"mm",sep="")   ,col="blue",cex=0.7)
    }
    
    lines(c(0,0),c(0,Range+Range/30) ,col="blue") 
    txtW <- strwidth(parse(text=paste("270", "^o ", sep="")))/2
    text((Range+txtW)*cos(pi/180 * seq(0,-270,-90) + pi/2)+Range/40,(Range+txtW)*sin(pi/180 *seq(0,-270,-90) + pi/2) ,labels = parse(text=paste(seq(0,270,90), "^o ", sep="")) ,col="blue",cex=0.8)
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