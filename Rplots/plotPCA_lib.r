
## Used for PCA 
standardizeHuntData <- function(datCapStat)
{
  within( datCapStat,{
    ###Assume split in High Low Values is around mean
    Efficiency_norm <- (Efficiency-mean(Efficiency))/sd(Efficiency) 
    HuntPower_norm    <- (HuntPower-mean(HuntPower)) /sd(HuntPower)
    CaptureSpeed_norm <-(CaptureSpeed-mean(CaptureSpeed))/sd(CaptureSpeed)
    Undershoot_norm    <- (Undershoot-1)/sd(Undershoot)
    DistanceToPrey_norm <- (DistanceToPrey-mean(DistanceToPrey))/sd(DistanceToPrey)
    CaptureAttempts_norm <- (CaptureEvents-mean(CaptureEvents)) /sd(CaptureEvents)
    ##Use Centre As The Mean Of The Most Efficient Hunters
    TimeToHitPrey_norm <-   (FramesToHitPrey/G_APPROXFPS - mean(datCapStat[datCapStat$Efficiency >0.5,]$FramesToHitPrey/G_APPROXFPS) ) /sd(FramesToHitPrey/G_APPROXFPS)
    DistSpeedUnder_norm <- (DistanceToPrey*CaptureSpeed*Undershoot_norm -mean(DistanceToPrey*CaptureSpeed*Undershoot_norm))/sd(DistanceToPrey*CaptureSpeed*Undershoot_norm)
    DistUnder_norm  <- (DistanceToPrey_norm*Undershoot_norm - mean(DistanceToPrey_norm*Undershoot_norm)) /sd(DistanceToPrey_norm*Undershoot_norm)   
    DistSpeed_norm <- (DistanceToPrey*CaptureSpeed -mean(DistanceToPrey*CaptureSpeed))/sd(DistanceToPrey*CaptureSpeed)
    SpeedUnder_norm <- (CaptureSpeed*Undershoot_norm -mean(CaptureSpeed*Undershoot_norm))/sd(CaptureSpeed*Undershoot_norm)
    
  })
}


plotPCAPerHunter <- function(datHunterStat_norm,strfilename)
{
  
  ## Set Colours
  require("graphics")
  colClass <- c("#00AFBB", "#E7B800", "#FC4E07")
  
  colEfficiency <- hcl.colors(12, alpha = 1, rev = FALSE) # heat.colors rainbow(12)
  colFactrAxes <- hcl.colors(6,palette="RdYlBu")
  colAxis <- c("#00AFBB", "#E7B800", "#FC4E07",colFactrAxes[1],colFactrAxes[2],colFactrAxes[3])
  colourGroup <- c(colourLegL[1],colourLegL[2],colourLegL[3])
  pchLPCA <- c(16,17,15)
  
  datPCAHunter_norm <- data.frame( with(datHunterStat_norm,{ #,'DL','NL' mergedCapDat$HuntPower < 5
    cbind(Efficiency=Efficiency_norm, #1
          #HuntPower, #2 ## Does not CoVary With Anyhting 
          #Group=groupID, #3
          DistanceToPrey=DistanceToPrey_norm, #4
          CaptureSpeed_norm, #5
          Undershoot_norm, #6
          DistSpeedProd=DistSpeed_norm, #7
          #DistUnderProd=DistUnder_norm, #8
          #SpeedUnderProd=SpeedUnder_norm, #9
          TimeToHitPrey=TimeToHitPrey_norm, #10
          Attempts=CaptureAttempts_norm
          
          #Cluster=Cluster#11
    )                                   } )          )
  
  
  pca_Hunter_norm <- prcomp(datPCAHunter_norm,scale.=FALSE)
  summary(pca_Hunter_norm)
  pcAxis <- c(1,2,1)
  rawHd <- pca_Hunter_norm$x[,pcAxis]
  
  biplot(pca_Hunter_norm,choices=c(1,2))
  
  densNL <-  kde2d(rawHd[,1][datHunterStat_norm$groupID == 3], rawHd[,2][datHunterStat_norm$groupID == 3],n=80)
  densLL <-  kde2d(rawHd[,1][datHunterStat_norm$groupID == 2], rawHd[,2][datHunterStat_norm$groupID == 2],n=80)
  densDL <-  kde2d(rawHd[,1][datHunterStat_norm$groupID == 1], rawHd[,2][datHunterStat_norm$groupID == 1],n=80)
  
  
  pdf(file= paste(strPlotExportPath,strfilename,sep=""),width=7,height=7)
  
  xplotRange = xlim=c(-2,3)
  yplotRange = ylim=c(-2.0,3)
  
  ## bottom, left,top, right
  par(mar = c(4.3,4.3,2,1))
  
  plot(rawHd[,1], rawHd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       #col=colEfficiency[round(datHunterStat$Efficiency*10)], pch=pchLPCA[as.numeric(datHunterStat$groupID) ],
       col=colourGroup[datHunterStat_norm$groupID ], pch=pchLPCA[as.numeric(datHunterStat_norm$groupID)],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID)],
       #xlab="PC1",ylab="PC2",
       #xlim=c(-3.5,4.2),,
       xlim=xplotRange,ylim=yplotRange,
       xlab=NA,ylab=NA,
       asp=1,
       cex=cex/1.4,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  contour(densNL,add=TRUE,col=colourGroup[3],nlevels=4,lwd=2,lty= 2, xlim=xplotRange,ylim=yplotRange)
  contour(densLL,add=TRUE,col=colourGroup[2],nlevels=4,lwd=2,lty= 1, xlim=xplotRange,ylim=yplotRange)
  contour(densDL,add=TRUE,col=colourGroup[1],nlevels=4,lwd=2,lty= 3, xlim=xplotRange,ylim=yplotRange)
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[2,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[2,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[1],lwd=3)
  text(1.1*scaleV*pca_Hunter_norm$rotation[2,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       0.2*scaleV*pca_Hunter_norm$rotation[2,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[3,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[3,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[2],lwd=3,lty=1)
  text(1.2*scaleV*pca_Hunter_norm$rotation[3,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       0.1+0.8*scaleV*pca_Hunter_norm$rotation[3,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[2],labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_Hunter_norm$rotation[4,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         -scaleV*pca_Hunter_norm$rotation[4,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col="black",lty=1,lwd=2)
  text(-1.5*scaleV*pca_Hunter_norm$rotation[4,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       -1.0*scaleV*pca_Hunter_norm$rotation[4,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col="black",labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[6,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[6,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[3],lty=1,lwd=3)
  text(0.8*scaleV*pca_Hunter_norm$rotation[6,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       1.1*scaleV*pca_Hunter_norm$rotation[6,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[3],labels="t Prey")
  
  
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[4],lty=1,lwd=2)
  text(1.3*scaleV*pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       1.1*scaleV*pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,
       col=colAxis[4],labels="Efficiency")
  
  ##Attempts To Capture   Component Projection
  arrows(0,0,scaleV*pca_Hunter_norm$rotation[7,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
         scaleV*pca_Hunter_norm$rotation[7,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[5],lty=1,lwd=3)
  text(0.8*scaleV*pca_Hunter_norm$rotation[7,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
       1.1*scaleV*pca_Hunter_norm$rotation[7,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col=colAxis[5],labels="Attempts")
  
  
  
  ###Heat Map Scale
  #posLeg <- c(3,-3) 
  #points(seq(posLeg[1],posLeg[1]+2,2/10),rep(posLeg[2],11),col=colEfficiency,pch=15,cex=3)
  #text(posLeg[1]-0.1,posLeg[2]+0.3,col="black",labels= prettyNum(min(datHunterStat$Efficiency),digits=1,format="f" ),cex=cex)
  #text(posLeg[1]+1,posLeg[2]+0.3,col="black",labels= prettyNum(max(datHunterStat$Efficiency)/2,digits=1,format="f" ),cex=cex)
  #text(posLeg[1]+2,posLeg[2]+0.3,col="black",labels= prettyNum(max(datHunterStat$Efficiency),digits=1,format="f" ),cex=cex)
  #max(mergedCapDat_filt$Efficiency)/2
  # 
  
  legend("topleft", legend=c(  expression (),
                               bquote(NF[""] ~ '#' ~ .(NROW(datHunterStat_norm[datHunterStat_norm$groupID == 3, ]))  ),
                               bquote(LF[""] ~ '#' ~ .(NROW(datHunterStat_norm[datHunterStat_norm$groupID == 2, ]))  ),
                               bquote(DF[""] ~ '#' ~ .(NROW(datHunterStat_norm[datHunterStat_norm$groupID == 1, ]))  )
                               
                               #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
  ),
  pch=c(pchLPCA[1],pchLPCA[2],pchLPCA[3]),lty=c(3,1,2),
  col=c(colourGroup[1],colourGroup[2],colourGroup[3]) )## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  ##legend("bottomright",legend=c("Slow","Fast"),fill=colClass, col=colClass,title="Cluster")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_Hunter_norm$sdev)
  pcEffVar <- ((pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]])^2 + (pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_Hunter_norm$rotation[1,][1:nComp]*pca_Hunter_norm$sdev[1:nComp])^2)
  
  #title(NA,sub=paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), 
  #                   " Coeff. variation:",prettyNum(sd(datHunterStat$Efficiency)/mean(datHunterStat$Efficiency) ,digits=2)) )
  message("Captured Variance ",prettyNum( 100*(pca_Hunter_norm$sdev[pcAxis[1]]^2 + pca_Hunter_norm$sdev[pcAxis[2]]^2) /sum( pca_Hunter_norm$sdev ^2),digits=3,format="f" ),"%" )
  message(paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(datHunterStat_norm$Efficiency)/mean(datHunterStat_norm$Efficiency) ,digits=2)))
  
  dev.off()
  
  return (pca_Hunter_norm)
} ##End Of Plot PCA
