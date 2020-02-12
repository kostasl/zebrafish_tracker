## Organize Manuscript Figures ### 
#### Kostas Lagogiannis 2019 
## \brief Make a scipt clarifying the script files used to produce each figure Used in the MS 


##Fig 1  Epxperimental TimeLine: manually designed  Inkscape

### Fig 2A ####
## The kinematics was produced by selecting one of the figure produced
## from Hunt event analysis loop in : runHuntEpisodeAnalysis.r
## 

### Fig 2B ####
## 
#source("Stats/stat_HuntRateInPreyRange.R")
#source("Stats/stat_HuntDuration.R")

### Fig 3 ####
#source("Stats/stat_HuntEfficiency.r")

### Fig 4 ####
#source("DataLabelling/plotLabelledDataResults.R")
#source("Stats/stat_CaptureSpeedVsDistanceToPrey.R")

### Fig 5 ####
#source("Stats/stat_LinRegression_TurnVsBearing.R")

## Fig 6  clustering Speed, TurnRatio and Distance to Prey using 3D 2xGaussian mixture method####
#source("Stats/stat_ClusterCaptureSpeedVsUndershootAndDistance.r")

### Fig 7 Show Covariance using Gaussian 3D non clustering model aong wit ####
#source("Stats/stat_CaptureSpeedVsUndershootAndDistance.r")




library(tools)
library(RColorBrewer);
library("MASS");
library(extrafont) ##For F
library(mvtnorm)

library(ggplot2) ##install.packages("ggplot2")
library(ggExtra)##  install.packages("ggExtra") ##devtools::install_github("daattali/ggExtra").
library(cowplot)
library(ggpubr) ##install.packages("ggpubr")


source("config_lib.R")
#setEnvFileLocations("OFFICE") #OFFICE,#LAPTOP HOME

source("DataLabelling/labelHuntEvents_lib.r")

####################
#source("TrackerDataFilesImport.r")
### Hunting Episode Analysis ####


###Used for drawing contour in ggplot -
## Draw the model fit above the cluster points
getFastClusterGrid <- function(drawMCMC)
{
  ### Add the cluster contours ###
  xran <- seq(0,0.8,0.05) ##Distance Grid
  yran <- seq(0,70,1) ##Speed Grid
  
  ##Example Code for PLotting Inferred Multivariate Cluster
  nsteps <- NROW(drawMCMC$mID[,,1][1,])
  mat_cov_fast <- rowMeans(drawMCMC$cov[2,,,(nsteps-1000):nsteps,1],dim=2) ##Average over samples
  
  mat_mu <- rowMeans(drawMCMC$mu[,,(nsteps-1000):nsteps,1],dim=2)
  #valGrid <- matrix( expand.grid(distance=xran,speed=yran),nrow=NROW(xran),ncol=NROW(yran) )
  valGrid <- expand.grid(distance=xran,speed=yran)
  #matrix(valGrid$distance,ncol=NROW(xran))
  cluster_z <- mvtnorm::dmvnorm(valGrid,mean=mat_mu[2,],sigma=mat_cov_fast )
  cluster_fast <- cbind(valGrid,cluster_z, factor(rep(16,NROW(valGrid) ),levels=c(1,16),labels=c("slow","fast")  ) )
  
  names(cluster_fast) <- c("DistanceToPrey", "CaptureSpeed", "Density","Cluster")
  
  return(cluster_fast)
}

###Used for drawing contour in ggplot
getSlowClusterGrid <- function(drawMCMC)
{
  ### Add the cluster contours ###
  xran <- seq(0,0.8,0.05) ##Distance Grid
  yran <- seq(0,70,1) ##Speed Grid
  
  ##Example Code for PLotting Inferred Multivariate Cluster
  nsteps <- NROW(drawMCMC$mID[,,1][1,])
  mat_cov_slow <- rowMeans(drawMCMC$cov[1,,,(nsteps-1000):nsteps,1],dim=2) ##Average over samples
  
  mat_mu <- rowMeans(drawMCMC$mu[,,(nsteps-1000):nsteps,1],dim=2)
  #valGrid <- matrix( expand.grid(distance=xran,speed=yran),nrow=NROW(xran),ncol=NROW(yran) )
  valGrid <- expand.grid(distance=xran,speed=yran)
  #matrix(valGrid$distance,ncol=NROW(xran))
  cluster_z <- mvtnorm::dmvnorm(valGrid,mean=mat_mu[1,],sigma=mat_cov_slow )
  cluster_slow <- cbind(valGrid,cluster_z, factor(rep(1,NROW(valGrid) ),levels=c(1,16),labels=c("slow","fast")  ) )
  
  names(cluster_slow) <- c("DistanceToPrey", "CaptureSpeed", "Density","Cluster")
  
  return(cluster_slow)
}

## Used for PCA 
standardizeHuntData <- function(datCapStat)
{
  within( datCapStat,{
    ###Assume split in High Low Values is around mean
    Efficiency_norm <- (Efficiency-mean(Efficiency))/sd(Efficiency) 
    HuntPower_norm    <- (HuntPower-mean(HuntPower)) /sd(HuntPower)
    CaptureSpeed_norm <-(CaptureSpeed-mean(CaptureSpeed))/sd(CaptureSpeed)
    DistSpeed_norm <- (DistanceToPrey*CaptureSpeed -mean(DistanceToPrey*CaptureSpeed))/sd(DistanceToPrey*CaptureSpeed)
    DistanceToPrey_norm <- (DistanceToPrey-mean(DistanceToPrey))/sd(DistanceToPrey)
    Undershoot_norm    <- (Undershoot-1)/sd(Undershoot)
    ##Use Centre As The Mean Of The Most Efficient Hunters
    TimeToHitPrey_norm <-   (FramesToHitPrey/G_APPROXFPS - mean(datCapStat[datCapStat$Efficiency >0.5,]$FramesToHitPrey/G_APPROXFPS) ) /sd(FramesToHitPrey/G_APPROXFPS)
    DistUnder_norm  <- (DistanceToPrey_norm*Undershoot_norm - mean(DistanceToPrey_norm*Undershoot_norm)) /sd(DistanceToPrey_norm*Undershoot_norm)   
    
  })
}



datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated.rds",sep="")) ##Original :huntEpisodeAnalysis_FirstBoutData_Validated

#### Plot Raw Capture Data Indicating Low/High Speed Clustering for each
### Load Pre Calc RJAgs Model Results
##   stat_CaptSpeedVsDistance_RJags.RData ##stat_CaptSpeedCluster_RJags.RData
load(file =paste(strDataExportDir,"stat_CaptSpeedVsDistance_RJags.RData",sep=""))

#### LOAD Capture First-Last Bout hunting that include the cluster classification - (made in stat_CaptureSpeedVsDistanceToPrey)
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_NL_2Dclustered.rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_LL_2Dclustered.rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_DL_2Dclustered.rds",sep="")) 

datHuntLabelledEventsSB <- getLabelledHuntEventsSet()
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)

### PREAMP DONE####
##################












##############Clustered  Capture Speed Vs Turn Ratio #### 
#### GGPLOT VERSION ###

pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_clusterCaptureSpeedVsDistToPrey_NF.pdf",sep=""),width=7,height=7)
 #layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
 #par(mar = c(3.9,4.7,12,1))

  p_NF = ggplot( datCapture_NL, aes(DistanceToPrey, CaptureSpeed,color =Cluster,fill=Cluster))  +
    ggtitle(NULL) +
    theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),
          axis.text = element_text(family="Helvetica",face="bold", size=16),
          plot.margin = unit(c(1,1,1,1), "mm")) + fill_palette("jco") +
          theme( ##Add the Legend
            legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6)
          ) 
          
    
  
  p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_NL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
    scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") )
    
  contour_fast <- getFastClusterGrid(draw_NF) ## Draw the mvtnorm model fit contour
  contour_slow <- getSlowClusterGrid(draw_NF)
  
  
  
  p_NF = p_NF +
    geom_contour(contour_fast, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
    geom_contour(contour_slow, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
    scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
    scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 80)) 
    #theme_linedraw()

    ggMarginal(p_NF, x="DistanceToPrey",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
  
    
    ##Make Custom Marginal Plot
   #xplot <- ggdensity(datCapture_NL,"DistanceToPrey",  mapping=aes(x="DistanceToPrey",color=datCapture_NL$Cluster), fill = "Cluster") +
    #clean_theme()  + theme(plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
    #fill_palette("jco") ##scale_color_manual( values = c("#00AFBB", "#00AFBB" ) )+
   #yplot <- ggdensity(datCapture_NL, "CaptureSpeed",mapping=aes(x="CaptureSpeed",color=datCapture_NL$Cluster), fill = "Cluster")+
  #  rotate() + clean_theme() + theme(plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none")+
  #  fill_palette("jco") #
  
    ##Cordinates Run 0-1 From lower left 0,0 
   # ggdraw() +
  #  draw_plot(xplot, x = 0.055, y = 0.8, width = 0.74, height = 0.2) +
   # draw_plot(p_NF, x = 0, y = 0, width = 0.8, height = 0.8) +
    #draw_plot(yplot, x = 0.8 , y = 0.055, width = 0.2, height = 0.74)

dev.off()

pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_clusterCaptureSpeedVsDistToPrey_LF.pdf",sep=""),width=7,height=7)

  p_LF <- ggplot( datCapture_LL, aes(DistanceToPrey, CaptureSpeed,color =Cluster,fill=Cluster)) + ggtitle(NULL)  +
    theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),
          axis.text = element_text(family="Helvetica",face="bold", size=16),
          plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
    fill_palette("jco")
  
  p_LF <- p_LF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_LL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
    scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") )
  contour_fast <- getFastClusterGrid(draw_LF)
  contour_slow <- getSlowClusterGrid(draw_LF)
  p_LF = p_LF + geom_contour(contour_fast, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
                geom_contour(contour_slow, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
                scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
                scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 80))
  ggMarginal(p_LF ,
             x="DistanceToPrey",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 

dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_clusterCaptureSpeedVsDistToPrey_DF.pdf",sep=""),width=7,height=7)

  p_DF = ggplot( datCapture_DL, aes(DistanceToPrey, CaptureSpeed,color =Cluster,fill=Cluster)) + ggtitle(NULL)  +
    theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),
          axis.text = element_text(family="Helvetica",face="bold", size=16),
          plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
    fill_palette("jco")
  
  p_DF = p_DF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_DL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) + 
                scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) # scale_color_manual( values = c(colourHPoint[4],colourHPoint[1])  )
  contour_fast <- getFastClusterGrid(draw_DF)
  contour_slow <- getSlowClusterGrid(draw_DF)
  p_DF = p_DF + geom_contour(contour_fast, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
                geom_contour(contour_slow, mapping = aes(x = DistanceToPrey, y = CaptureSpeed, z = Density) ,linetype=2 ) +
                scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
                scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 80)) 
  ggMarginal(p_DF ,x="DistanceToPrey",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 

dev.off()


### Probability of Membership in High speed Cluster / 
pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_ProbOfFactCapture_ggplot.pdf",sep=""),width=7,height=7)
 
  dat2_NF <- rbind( data.frame(D=tail(draw_NF$pS[,,1],500),Group=rep("NF",500) ),
                    data.frame(D=tail(draw_LF$pS[,,1],500),Group=rep("LF",500) ),
                    data.frame(D=tail(draw_DF$pS[,,1],500),Group=rep("DF",500) ))
   par(mar=c(5,5,5,5))
    ggplot(dat2_NF, aes(x=D),group=Group ) +
    #geom_density( lwd=1.5,aes(linetype=label,colour=label) ) +  
    geom_line(stat="density",lwd=1.5,show.legend=T,aes(linetype=Group,colour=Group) ) +
    theme(legend.position = c(0.1, 0.8),legend.title=element_blank(),legend.key.width = unit(3, "line")  ) + ## No Legend
    scale_x_continuous(name= expression(paste("Probability of high speed capture  ["~p["s"]~"]" )), limits=c(0, 1),expand=c(0,0) ) +
    scale_y_continuous(name="Density function", limits=c(0, 15),expand=c(0,0)) 
    
    #  guides( size = guide_legend(order = 3) )
    
    #scale_color_manual(labels=dat2_NF$label,values=colourHLine) ##Change legend text
    #
  
  
  #plot_probM = plot_probM + geom_density(dat2_NF[dat2_NF$label=="NF",], mapping=aes(x=pS_LF,colour=colourHLine[2])) 
  #plot_probM + geom_density(dat_pS, mapping=aes(x=pS_DF,colour=colourHLine[3]))  +
   #             scale_color_manual(labels = c("T999", "T888","T88asd8"),values=colourHLine) ##Change legend text
  
#####
dev.off()
  
  
pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_clusterMembership.pdf",sep=""),width=7,height=7)
  par(mar = c(3.9,4.7,1,1))
  #### ## Probability Density of Strike capture ####
  plot(density(tail(draw_NF$pS[,,1],1000),pBw=0.05),col=colourLegL[1],xlim=c(0,1),ylim=c(0.4,10),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex )
  lines(density(tail(draw_LF$pS[,,1],1000)),col=colourLegL[2],lwd=3,lty=2)
  lines(density(tail(draw_DF$pS[,,1],1000)),col=colourLegL[3],lwd=3,lty=3)
  #lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)
  mtext(side = 1,cex=cex, line = lineXAxis, expression(paste(bold("Probability of high speed capture  ["~p["s"]~"]" ) ) ) ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, expression(bold("Density function" ) ) )
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NROW(datCapture_NL$DistanceToPrey))  ),
                    bquote(LF[""] ~ '#' ~ .(NROW(datCapture_LL$DistanceToPrey))  ),
                    bquote(DF[""] ~ '#' ~ .(NROW(datCapture_DL$DistanceToPrey))  )
                    #bquote(ALL ~ '#' ~ .(ldata_ALL$N)  ) 
         ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)
  
  
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig5_stat_meanDistanceOfFastCapture.pdf",sep=""),width=7,height=7)
  par(mar = c(3.9,4.7,1,1))
  #### ## Probability Density of Strike capture ####
  plot(density(tail(draw_NF$mu[2,1,,1],1000)),col=colourLegL[1],xlim=c(0,0.6),ylim=c(0.0,35),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex )
  lines(density(tail(draw_LF$mu[2,1,,1],1000)),col=colourLegL[2],lwd=3,lty=2)
  lines(density(tail(draw_DF$mu[2,1,,1],1000)),col=colourLegL[3],lwd=3,lty=3)
  #lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)
  mtext(side = 1,cex=cex, line = lineXAxis, expression(paste(bold("Estimated mean distance of high speed capture (mm)") )  ) ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, expression(bold("Density function")  ))
  
  # legend("topleft",
  #        legend=c(  expression (),
  #                   bquote(NF[""] ~ '#' ~ .(NROW(datCapture_NL$DistanceToPrey))  ),
  #                   bquote(LF[""] ~ '#' ~ .(NROW(datCapture_LL$DistanceToPrey))  ),
  #                   bquote(DF[""] ~ '#' ~ .(NROW(datCapture_DL$DistanceToPrey))  )
  #                   #bquote(ALL ~ '#' ~ .(ldata_ALL$N)  ) 
  #        ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
  #        col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)


dev.off()

##########UNDERSHOOT Vs Distance
#fig6.CaptureSpeed/fig6-stat_modelCaptureSpeedVsUndershootAndDistance_Valid.pdf

pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndDistance_NF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))

p_NF = ggplot( datCapture_NL, aes(Undershoot, DistanceToPrey ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_NL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_y_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
  scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 

ggMarginal(p_NF, x="Undershoot",y="DistanceToPrey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndDistance_LF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))

p_NF = ggplot( datCapture_LL, aes(Undershoot, DistanceToPrey ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_LL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_y_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
  scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 

ggMarginal(p_NF, x="Undershoot",y="DistanceToPrey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndDistance_DF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))

p_NF = ggplot( datCapture_DL, aes(Undershoot, DistanceToPrey ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_DL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_y_continuous(name="Distance to prey (mm)", limits=c(0, 0.8)) +
  scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 

ggMarginal(p_NF, x="Undershoot",y="DistanceToPrey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()

########UNdershoot - Speed ###


pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndSpeed_NF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))
  
  p_NF = ggplot( datCapture_NL, aes(Undershoot, CaptureSpeed ,color =Cluster,fill=Cluster)) +
    ggtitle(NULL) +
    theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
    fill_palette("jco") +
    theme( ##Add the Legend
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  
  
  p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_NL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
    scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
    scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 60)) +
    scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 
  
  p_NF = p_NF + geom_vline(xintercept = 1, linetype="dotted", color = "grey", size=1.0)
  
  ggMarginal(p_NF, x="Undershoot",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()


pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndSpeed_DF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))

p_DF = ggplot( datCapture_DL, aes(Undershoot, CaptureSpeed ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_DF = p_DF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_DL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 60)) +
  scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 

p_DF = p_DF + geom_vline(xintercept = 1, linetype="dotted", color = "grey", size=1.0)

ggMarginal(p_DF, x="Undershoot",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()



pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_UndershootAndSpeed_LF.pdf",sep=""),width=7,height=7)
#layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# ##Margin: (Bottom,Left,Top,Right )
#par(mar = c(3.9,4.7,12,1))

p_LF = ggplot( datCapture_LL, aes(Undershoot, CaptureSpeed ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_LF = p_LF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_LL$Cluster) ) +  xlim(0, 0.8) +  ylim(0, 80) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_y_continuous(name="Capture Speed (mm/sec)", limits=c(0, 60)) +
  scale_x_continuous(name="Turn ratio", limits=c(0, 2)) 
p_LF = p_LF + geom_vline(xintercept = 1, linetype="dotted", color = "grey", size=1.0)

ggMarginal(p_LF, x="Undershoot",y="CaptureSpeed", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()

plot(dens_dist_NF_all,xlim=c(0.0,0.5),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,5),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_dist_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_dist_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_dist_LF_all,xlim=c(0.0,0.5),col=colourLegL[2],lwd=4,lty=1)
lines(dens_dist_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_dist_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_dist_DF_all,xlim=c(0.0,0.5),col=colourLegL[3],lwd=4,lty=1)
lines(dens_dist_DF_fast,col=colourLegL[3],lwd=2,lty=2)
lines(dens_dist_DF_slow,col=colourLegE[3],lwd=2,lty=2)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(NROW(datCapture_NL))  ),
                  bquote(LF[""] ~ '#' ~ .(NROW(datCapture_LL))  ),
                  bquote(DF[""] ~ '#' ~ .(NROW(datCapture_DL))  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density ") )
mtext(side = 1,cex=cex, line = lineXAxis, expression(paste("Probability of high speed capture  ["~p["s"]~"]" ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)



###### Capture Speed  ###

lineAxis = 2.4
lineXAxis = 2.7
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.3,2,1))
##Make SPeed Density Of Each Cluster
dens_speed_NF_all <- density(datCapture_NL$CaptureSpeed)
dens_speed_NF_fast <- density(datCapture_NL$CaptureSpeed[lClustScore_NF$pchL == 16])
dens_speed_NF_slow <- density(datCapture_NL$CaptureSpeed[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_speed_LF_all <- density(datCapture_LL$CaptureSpeed)
dens_speed_LF_fast <- density(datCapture_LL$CaptureSpeed[lClustScore_LF$pchL == 16])
dens_speed_LF_slow <- density(datCapture_LL$CaptureSpeed[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_speed_DF_all <- density(datCapture_DL$CaptureSpeed)
dens_speed_DF_fast <- density(datCapture_DL$CaptureSpeed[lClustScore_DF$pchL == 16])
dens_speed_DF_slow <- density(datCapture_DL$CaptureSpeed[lClustScore_DF$pchL == 1])

## Plot Density Speed ##
plot(dens_speed_NF_all,xlim=c(0.0,60),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,0.1),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_speed_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_speed_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_speed_LF_all,xlim=c(0.0,60),col=colourLegL[2],lwd=4,lty=1,ylim=c(0,0.1))
lines(dens_speed_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_speed_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_speed_DF_all,xlim=c(0.0,60),col=colourLegL[3],lwd=4,lty=1,ylim=c(0,0.1))
lines(dens_speed_DF_fast,col=colourLegL[3],lwd=2,lty=2)
lines(dens_speed_DF_slow,col=colourLegE[3],lwd=2,lty=2)
####### END OF Speed ###

#'######### TURN RATIO ##########
##Make SPeed Density Of Each Cluster
dens_turn_NF_all <- density(datCapture_NL$Undershoot)
dens_turn_NF_fast <- density(datCapture_NL$Undershoot[lClustScore_NF$pchL == 16])
dens_turn_NF_slow <- density(datCapture_NL$Undershoot[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_turn_LF_all <- density(datCapture_LL$Undershoot)
dens_turn_LF_fast <- density(datCapture_LL$Undershoot[lClustScore_LF$pchL == 16])
dens_turn_LF_slow <- density(datCapture_LL$Undershoot[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
fracSlow_DF <- table(lClustScore_DF$pchL )[1]/NROW(lClustScore_DF$pchL)
dens_turn_DF_all <- density(datCapture_DL$Undershoot)
dens_turn_DF_fast <- density(datCapture_DL$Undershoot[lClustScore_DF$pchL == 16] )
dens_turn_DF_slow <- density(datCapture_DL$Undershoot[lClustScore_DF$pchL == 1])



## Plot TURN RATIO  ##
plot(dens_turn_NF_all,xlim=c(0.0,2),col=colourLegL[1],lwd=4,lty=1,ylim=c(0,3),
     main=NA,cex=cex,xlab=NA,ylab=NA)
lines(dens_turn_NF_fast,col=colourLegL[1],lwd=2,lty=2)
lines(dens_turn_NF_slow,col=colourLegE[1],lwd=2,lty=2)

plot(dens_turn_LF_all,xlim=c(0.0,2),col=colourLegL[2],lwd=4,lty=1,ylim=c(0,3))
lines(dens_turn_LF_fast,col=colourLegL[2],lwd=2,lty=2)
lines(dens_turn_LF_slow,col=colourLegE[2],lwd=2,lty=2)

plot(dens_turn_DF_all,xlim=c(0.0,2),col=colourLegL[3],lwd=4,lty=1,ylim=c(0,3))
lines(dens_turn_DF_fast$x,dens_turn_DF_fast$y*(1-fracSlow_DF),col=colourLegL[3],lwd=2,lty=2)
lines(dens_turn_DF_slow$x,dens_turn_DF_slow$y*fracSlow_DF,col=colourLegE[3],lwd=2,lty=2)
####### END OF Speed ###




pdf(file= paste(strPlotExportPath,"distal",strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.5,1,1))

plot(datCapture_NL$Undershoot, datCapture_NL$CaptureSpeed,col=colourLegL[1],pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datCapture_NL$CaptureSpeed ~ datCapture_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex)  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datCapture_LL$Undershoot, datCapture_LL$CaptureSpeed,col=colourLegL[2],pch=lClustScore_LF$pchL,
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datCapture_LL$CaptureSpeed ~ datCapture_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis-0.7, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datCapture_DL$Undershoot, datCapture_DL$CaptureSpeed,col=colourLegL[3],pch=lClustScore_DF$pchL,
     ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datCapture_DL$CaptureSpeed ~ datCapture_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
contour(densDL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


dev.off()



## EMPIRICAL - UNdeshoot vs Prey Distance 
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig7-UndershootDistanceCV_Distal_scatter.pdf",sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(4.5,4.3,0.5,1))

plot(datCapture_NL$Undershoot, datCapture_NL$DistanceToPrey,col=colourLegL[1],pch=lClustScore_NF$pchL,
     xlab=NA,ylab=NA,ylim=c(0,1.0),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datCapture_NL$DistanceToPrey ~ datCapture_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datCapture_LL$Undershoot, datCapture_LL$DistanceToPrey,col=colourLegL[2],pch=lClustScore_LF$pchL,
     ylim=c(0,1),xlim=c(0,2.0),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datCapture_LL$DistanceToPrey ~ datCapture_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
mtext(side = 2,cex=cex, line = 2.2, expression("Distance to prey  (mm) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datCapture_DL$Undershoot, datCapture_DL$DistanceToPrey,col=colourLegL[3],pch=lClustScore_DF$pchL,
     ylim=c(0,1.0),xlim=c(0,2),   xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datCapture_DL$DistanceToPrey ~ datCapture_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex) 

dev.off()


########## END oF CaptureSpeed vs Distance  ###




##Make SPeed Density Of Each Cluster
dens_dist_NF_all <- density(datCapture_NL$DistanceToPrey)
dens_dist_NF_fast <- density(datCapture_NL$DistanceToPrey[lClustScore_NF$pchL == 16])
dens_dist_NF_slow <- density(datCapture_NL$DistanceToPrey[lClustScore_NF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_LF_all <- density(datCapture_LL$DistanceToPrey)
dens_dist_LF_fast <- density(datCapture_LL$DistanceToPrey[lClustScore_LF$pchL == 16])
dens_dist_LF_slow <- density(datCapture_LL$DistanceToPrey[lClustScore_LF$pchL == 1])

##Make SPeed Density Of Each Cluster
dens_dist_DF_all <- density(datCapture_DL$DistanceToPrey)
dens_dist_DF_fast <- density(datCapture_DL$DistanceToPrey[lClustScore_DF$pchL == 16])
dens_dist_DF_slow <- density(datCapture_DL$DistanceToPrey[lClustScore_DF$pchL == 1])

outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 2.4
lineXAxis = 3.0
cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1


####################################################
##############    # GAPE-TIMING #    #############
##################################################
## These plots use the merged FirstBouts-CaptureData from : huntEpisodeAnalysis_FirstBoutData_wCapFrame_XX_clustered
##  That result after validation from code running validateCaptureStrikeData - which calculates time to hitPrey, speeds etc
## Once validated the data is re-merged with the clustering Bayssian results

pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_DistanceVsTimeToHitPrey_NF.pdf",sep=""),width=7,height=7)
p_NF = ggplot( datCapture_NL, aes( DistanceToPrey,FramesToHitPrey/G_APPROXFPS ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_NF = p_NF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_NL$Cluster) ) +  xlim(0, 0.6) +  ylim(0, 0.4) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.6)) +
  scale_y_continuous(name="Time to reach prey", limits=c(0, 0.4)) 

ggMarginal(p_NF, x="Distance to prey",y="Time to reach prey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()

## time to reach min dist to prey 
pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_DistanceVsTimeToHitPrey_LF.pdf",sep=""),width=7,height=7)
p_LF = ggplot( datCapture_LL, aes( DistanceToPrey,FramesToHitPrey/G_APPROXFPS ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_LF = p_LF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_LL$Cluster) ) +  xlim(0, 0.6) +  ylim(0, 0.4) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.6)) +
  scale_y_continuous(name="Time to reach prey", limits=c(0, 0.4)) 

ggMarginal(p_LF, x="Distance to prey",y="Time to reach prey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()


## time to reach min dist to prey 
pdf(file= paste(strPlotExportPath,"/stat/fig6_stat_DistanceVsTimeToHitPrey_DF.pdf",sep=""),width=7,height=7)
p_DF = ggplot( datCapture_DL, aes( DistanceToPrey,FramesToHitPrey/G_APPROXFPS ,color =Cluster,fill=Cluster)) +
  ggtitle(NULL) +
  theme(axis.title =  element_text(family="Helvetica",face="bold", size=16),plot.margin = unit(c(1,1,1,1), "mm"), legend.position = "none") +
  fill_palette("jco")

p_DF = p_DF + geom_point( size = 3, alpha = 0.6,aes(color =datCapture_DL$Cluster) ) +  xlim(0, 0.6) +  ylim(0, 0.4) +
  scale_color_manual( values = c("#00AFBB", "#E7B800", "#FC4E07") ) +
  scale_x_continuous(name="Distance to prey (mm)", limits=c(0, 0.6)) +
  scale_y_continuous(name="Time to reach prey", limits=c(0, 0.4)) 

ggMarginal(p_DF, x="Distance to prey",y="Time to reach prey", type = "density",groupColour = TRUE,groupFill=TRUE,show.legend=TRUE) 
dev.off()



#######################################
########## PCA  - FACTOR ANALYSIS ####
#######################################
# Check Correlation Of UNdershoot With Hunt POwer
##Take all expID from the successful hunt Events we have extracted hunt variables from 
vexpID <- list(LF = datTrackedEventsRegister[datCapture_LL$RegistarIdx,]$expID,
               NF=datTrackedEventsRegister[datCapture_NL$RegistarIdx,]$expID,
               DF=datTrackedEventsRegister[datCapture_DL$RegistarIdx,]$expID)

#datFishSuccessRate[datFishSuccessRate$expID %in% vexpID$LF, ]$HuntPower

## Add Exp ID Column - Signifying Which Larvae Executed the Capture Success Hunt- 
datCapture_LF_wExpID <- cbind(datCapture_LL,expID=vexpID$LF)
datCapture_NF_wExpID <- cbind(datCapture_NL,expID=vexpID$NF)
datCapture_DF_wExpID <- cbind(datCapture_DL,expID=vexpID$DF)


##Merge Hunt Power To Hunt-Capture Variables 
datMergedCapAndSuccess_LF <- merge(x=datCapture_LF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)
datMergedCapAndSuccess_NF <- merge(x=datCapture_NF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)
datMergedCapAndSuccess_DF <- merge(x=datCapture_DF_wExpID,y=datFishSuccessRate,by="expID",all.x=TRUE)

## Merge 
mergedCapDat <- rbind(datMergedCapAndSuccess_LF,datMergedCapAndSuccess_DF,datMergedCapAndSuccess_NF)

mergedCapDat$groupID <- as.factor(mergedCapDat$groupID)
groupLabels <- levels(mergedCapDat$groupID)
## Now Compile  Behaviour data Per Larvae
mergedCapDat$groupID <- as.numeric(mergedCapDat$groupID)
mergedCapDat_mod<-mergedCapDat ##Temp Copy
mergedCapDat_mod$expID <- as.numeric(as.character(mergedCapDat_mod$expID))
datHunterStat <- aggregate(mergedCapDat_mod,by=list(mergedCapDat_mod$expID),mean)

##Recover Group ID Factor
datHunterStat$groupIDF <- levels(datTrackedEventsRegister$groupID)[datHunterStat$groupID]

##Error Check Assert - Check IDs Have been Matched
  stopifnot(datHunterStat[datHunterStat$groupIDF == 'DL',]$expID %in% unique(datTrackedEventsRegister[datTrackedEventsRegister$groupID == 'DL',]$expID))
  stopifnot(datHunterStat[datHunterStat$groupIDF == 'LL',]$expID %in% unique(datTrackedEventsRegister[datTrackedEventsRegister$groupID == 'LL',]$expID))
###

datHunterStat <- standardizeHuntData(datHunterStat)
mergedCapDat <- standardizeHuntData(mergedCapDat)
# Show Stdandardized Efficiency Distribution 
#hist(datHunterStat$Efficiency_norm )


## Set Colours
require("graphics")
colClass <- c("#00AFBB", "#E7B800", "#FC4E07")
colEfficiency <- hcl.colors(12, alpha = 1, rev = FALSE) # heat.colors rainbow(12)
colFactrAxes <- hcl.colors(6,palette="RdYlBu")
colourGroup <- c(colourLegL[3],colourLegL[2],colourLegL[1])
### PCA ANalysis Of Variance - Finding the Factors That contribute to efficiency
## ##Make MAtrix
##Also CHeck OUt varimax and factanal

### PCA ANalysis Of Variance - Finding the Factors That contribute to efficiency
## ##Make MAtrix

######### Show PCA For Hunter '#####

datPCAHunter_norm <- data.frame( with(datHunterStat,{ #,'DL','NL' mergedCapDat$HuntPower < 5
  cbind(Efficiency=Efficiency_norm, #1
        #HuntPower, #2 ## Does not CoVary With Anyhting 
        #Group=groupID, #3
        DistanceToPrey=DistanceToPrey_norm, #4
        CaptureSpeed_norm, #5
        Undershoot_norm, #6
        DistSpeedProd=DistSpeed_norm, #7
        #DistUnderProd=DistUnder_norm, #8
        #SpeedUnderProd=SpeedUnder_norm, #9
        TimeToHitPrey=TimeToHitPrey_norm #10
        #Cluster=Cluster#11
  )                                   } )          )

  
  ###
  pca_Hunter_norm <- prcomp(datPCAHunter_norm,scale.=FALSE)
  summary(pca_Hunter_norm)
  pcAxis <- c(1,2,3)
  rawHd <- pca_Hunter_norm$x[,pcAxis]
  
  biplot(pca_Hunter_norm,choices=c(1,2))

  densNL <-  kde2d(rawHd[,1][datHunterStat$groupID == 3], rawHd[,2][datHunterStat$groupID == 3],n=80)
  densLL <-  kde2d(rawHd[,1][datHunterStat$groupID == 2], rawHd[,2][datHunterStat$groupID == 2],n=80)
  densDL <-  kde2d(rawHd[,1][datHunterStat$groupID == 1], rawHd[,2][datHunterStat$groupID == 1],n=80)
  
  
  scaleV <- 2
  scaleVE <- scaleV
  
  
  arrow_Efficiency <- c(scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
                        scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,
                        scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[3]]*pca_Hunter_norm$sdev[pcAxis[3]]^2)
  
  
  ##This plot is repeated in plotPCA_analysis 
pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntersBehaviourPC1_2_GroupColour_ALL.pdf",sep=""),width=7,height=7)
  ## bottom, left,top, right
  par(mar = c(5.9,4.3,2,1))
  
  plot(rawHd[,1], rawHd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       #col=colEfficiency[round(datHunterStat$Efficiency*10)], pch=pchLPCA[as.numeric(datHunterStat$groupID) ],
       col=colourGroup[datHunterStat$groupID ], pch=pchLPCA[as.numeric(datHunterStat$groupID)],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID)],
       #xlab="PC1",ylab="PC2",
       xlim=c(-4.2,4.2),ylim=c(-3.0,3.2),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  contour(densNL,add=TRUE,col=colourGroup[1],nlevels=4,lwd=2,lty= 1)
  contour(densLL,add=TRUE,col=colourGroup[2],nlevels=4,lwd=2,lty= 2)
  contour(densDL,add=TRUE,col=colourGroup[3],nlevels=4,lwd=2,lty= 3)
  
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[2,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[2,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.7*scaleV*pca_Hunter_norm$rotation[2,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       1.7*scaleV*pca_Hunter_norm$rotation[2,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[3,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[3,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],lwd=3,lty=1)
  text(1.2*scaleV*pca_Hunter_norm$rotation[3,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       0.1+0.8*scaleV*pca_Hunter_norm$rotation[3,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_Hunter_norm$rotation[4,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         -scaleV*pca_Hunter_norm$rotation[4,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col="black",lty=1,lwd=2)
  text(-1.9*scaleV*pca_Hunter_norm$rotation[4,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       -1.0*scaleV*pca_Hunter_norm$rotation[4,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col="black",labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[6,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[6,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=1,lwd=3)
  text(0.8*scaleV*pca_Hunter_norm$rotation[6,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       1.3*scaleV*pca_Hunter_norm$rotation[6,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection

  arrows(0,0,arrow_Efficiency[1],
         arrow_Efficiency[2],col="blue",lty=1,lwd=2)
  text(1.8*scaleV*pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       0.8*scaleV*pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,
       col="blue",labels="Efficiency")
  
  ###Heat Map Scale
    #posLeg <- c(3,-3) 
    #points(seq(posLeg[1],posLeg[1]+2,2/10),rep(posLeg[2],11),col=colEfficiency,pch=15,cex=3)
    #text(posLeg[1]-0.1,posLeg[2]+0.3,col="black",labels= prettyNum(min(datHunterStat$Efficiency),digits=1,format="f" ),cex=cex)
    #text(posLeg[1]+1,posLeg[2]+0.3,col="black",labels= prettyNum(max(datHunterStat$Efficiency)/2,digits=1,format="f" ),cex=cex)
    #text(posLeg[1]+2,posLeg[2]+0.3,col="black",labels= prettyNum(max(datHunterStat$Efficiency),digits=1,format="f" ),cex=cex)
    #max(mergedCapDat_filt$Efficiency)/2
  # 
  
  legend("topleft", legend=c(  expression (),
                               bquote(DF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 1, ]))  ),
                               bquote(LF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 2, ]))  ),
                               bquote(NF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 3, ]))  )
                               #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
        ),
         pch=pchLPCA,
         col=colourGroup)## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  ##legend("bottomright",legend=c("Slow","Fast"),fill=colClass, col=colClass,title="Cluster")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_Hunter_norm$sdev)
  pcEffVar <- ((pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]])^2 + (pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_Hunter_norm$rotation[1,][1:nComp]*pca_Hunter_norm$sdev[1:nComp])^2)

  #title(NA,sub=paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), 
  #                   " Coeff. variation:",prettyNum(sd(datHunterStat$Efficiency)/mean(datHunterStat$Efficiency) ,digits=2)) )
  message("Captured Variance ",prettyNum( 100*(pca_Hunter_norm$sdev[pcAxis[1]]^2 + pca_Hunter_norm$sdev[pcAxis[2]]^2) /sum( pca_Hunter_norm$sdev ^2),digits=3,format="f" ),"%" )
  message(paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(datHunterStat$Efficiency)/mean(datHunterStat$Efficiency) ,digits=2)))
  
dev.off()



### MAKE A 3D View ###
library(rgl)

open3d()##mergedCapDat$groupID
rgl::plot3d( x=rawHd[,1], z=rawHd[,2], y=rawHd[,3], col = colourGroup[datHunterStat$groupID ] , type = "s", radius = 0.3,
             xlab="PC1", zlab="PC2",ylab="PC3",
             xlim=c(-5.,5), ylim=c(-5,5), zlim=c(-8,8),
             box = FALSE ,aspect = TRUE
             #,expand = 1.5
)
##Add Efficiency Axis
arrow3d(c(0,0,0),2*arrow_Efficiency,n=5, type = "rotation", col = "blue",thickness=3,width=0.5)
###END PCA PLOT ##

################# Fig 7, 3D Model Plot Is in stat_3DLarvaGroupBehaviour ###

