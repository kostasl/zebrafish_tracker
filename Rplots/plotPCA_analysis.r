### Oct 2019 : utility plot scripts produced during pCA analysis - Saved here for reference
#### PCA Analysis ##


## NOTE "Load required data as in main_GenerateMSFigures.r , Before Calling any plot section

#######################################
########## PCA  - FACTOR ANALYSIS ####
#######################################



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
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_NL_clustered.rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_LL_clustered.rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_DL_clustered.rds",sep="")) 

datHuntLabelledEventsSB <- getLabelledHuntEventsSet()
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)



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
colourGroup <- c(colourLegL[1],colourLegL[2],colourLegL[3])
pchLPCA <- c(16,17,15)



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



pca_Hunter_norm <- prcomp(datPCAHunter_norm,scale.=FALSE)
summary(pca_Hunter_norm)
pcAxis <- c(1,2,1)
rawHd <- pca_Hunter_norm$x[,pcAxis]

biplot(pca_Hunter_norm,choices=c(1,2))

densNL <-  kde2d(rawHd[,1][datHunterStat$groupID == 3], rawHd[,2][datHunterStat$groupID == 3],n=80)
densLL <-  kde2d(rawHd[,1][datHunterStat$groupID == 2], rawHd[,2][datHunterStat$groupID == 2],n=80)
densDL <-  kde2d(rawHd[,1][datHunterStat$groupID == 1], rawHd[,2][datHunterStat$groupID == 1],n=80)


###Change The Filter Here, Do PCA again and then Locate and plto group Specific
mergedCapDat_filt <- mergedCapDat #mergedCapDat[mergedCapDat$groupID == 'NL',]

datpolyFactor_norm <- data.frame( with(mergedCapDat_filt,{ #,'DL','NL' mergedCapDat$HuntPower < 5
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
        Cluster=Cluster#11
  )                                   } )          )


###
pca_norm <- prcomp(datpolyFactor_norm,scale.=FALSE)
summary(pca_norm)
pcAxis <- c(1,2,1)
rawd <- pca_norm$x[,pcAxis]

biplot(pca_norm,choices=c(1,2))


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
    
    contour(densNL,add=TRUE,col=colourGroup[3],nlevels=4,lwd=2,lty= 2)
    contour(densLL,add=TRUE,col=colourGroup[2],nlevels=4,lwd=2,lty= 1)
    contour(densDL,add=TRUE,col=colourGroup[1],nlevels=4,lwd=2,lty= 3)
    
    scaleV <- 2
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
    scaleVE <- scaleV
    arrows(0,0,scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[1]]*pca_Hunter_norm$sdev[pcAxis[1]]^2,
           scaleVE*pca_Hunter_norm$rotation[1,][pcAxis[2]]*pca_Hunter_norm$sdev[pcAxis[2]]^2,col="blue",lty=1,lwd=2)
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
                                 bquote(NF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 3, ]))  ),
                                 bquote(LF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 2, ]))  ),
                                 bquote(DF[""] ~ '#' ~ .(NROW(datHunterStat[datHunterStat$groupID == 1, ]))  )
                                 
                                 #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
    ),
    pch=c(pchLPCA[1],pchLPCA[2],pchLPCA[3]),lty=c(2,1,3),
    col=c(colourGroup[1],colourGroup[2],colourGroup[3]) )## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
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





###Change The Filter Here, Do PCA again and then Locate and plto group Specific
mergedCapDat_filt <- mergedCapDat #mergedCapDat[mergedCapDat$groupID == 'NL',]

datpolyFactor_norm <- data.frame( with(mergedCapDat_filt,{ #,'DL','NL' mergedCapDat$HuntPower < 5
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
        Cluster=Cluster#11
  )                                   } )          )


###
pca_norm <- prcomp(datpolyFactor_norm,scale.=FALSE)
summary(pca_norm)
pcAxis <- c(1,2,1)
rawd <- pca_norm$x[,pcAxis]

biplot(pca_norm,choices=c(1,2))


pchLPCA <- c(15,17,16)

pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC1_2_GroupColour_ALL.pdf",sep=""),width=7,height=7)
  ## bottom, left,top, right
  par(mar = c(4.2,4.3,1,1))
  
  plot(rawd[,1], rawd[,2],
       ######col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       #col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID) ],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       col=colourGroup[mergedCapDat_filt$groupID ], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.7*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.2*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],lwd=2,lty=3)
  text(0.8*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2,lwd=2)
  text(-1.5*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-1.0*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(0.8*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.2*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.0*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  # ##Heat Map Scale
  # points(seq(1,2,1/10),rep(-2,11),col=colEfficiency,pch=15,cex=3)
  # text(1,-1.8,col="black",labels= prettyNum(min(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  # text(1+0.5,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency)/2,digits=1,format="f" ),cex=cex)
  # text(2,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  # max(mergedCapDat_filt$Efficiency)/2
  # 
  
  legend("bottomleft",legend=c( paste0(groupLabels[3],' #',table(mergedCapDat_filt$groupID)[3]),
                                paste0(groupLabels[2],' #',table(mergedCapDat_filt$groupID)[2]),
                                paste0(groupLabels[1],' #',table(mergedCapDat_filt$groupID)[1]) 
  ),
  pch=c(pchLPCA[3],pchLPCA[2],pchLPCA[1]),
  col=c(colourGroup[3],colourGroup[2],colourGroup[1]) )## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  #legend("bottomright",legend=c("Slow","Fast"),fill=colClass, col=colClass,title="Cluster")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  
  title(NA,sub=paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(mergedCapDat_filt$Efficiency)/mean(mergedCapDat_filt$Efficiency) ,digits=2)) )
  message("Captured Variance ",prettyNum( (pca_norm$sdev[pcAxis[1]]^2 + pca_norm$sdev[pcAxis[2]]^2) /sum( pca_norm$sdev ^2),digits=3,format="f" ) )
  message(paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(mergedCapDat_filt$Efficiency)/mean(mergedCapDat_filt$Efficiency) ,digits=2)))

dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC1_2_EfficiencyColour_DF.pdf",sep=""),width=7,height=7)
  
  plot(rawd[,1], rawd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       #col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID) ],
       col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID)],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.8*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.5*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],lwd=2,lty=3)
  text(0.8*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2,lwd=2)
  text(-1.5*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-1.0*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(0.8*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.2*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.0+1.0*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  legend("bottomleft",legend=c("NF","LF","DF"),pch=c(pchLPCA[3],pchLPCA[2],pchLPCA[1]),
         col="black")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  legend("bottomright",legend=c("Slow","Fast"),fill=colClass, col=colClass,title="Cluster")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  title(NA,sub=paste("% of Efficiency Variance Explained : ",prettyNum( 100*pcEffVar/EffVar) ))
  

dev.off()








pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC1_2_EfficiencyColour_LF.pdf",sep=""),width=7,height=7)
  
  plot(rawd[,1], rawd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID) ],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchL[4+datpolyFactor_norm$Group],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+datpolyFactor_norm$Group],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.8*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.5*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lwd=2,lty=3)
  text(1.5*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  text(-0.6*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-1.4*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(-0.06*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.8*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.1*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  legend("bottomleft",legend=c("LF","NF","DF"),pch=c(pchLPCA[2],pchLPCA[3],pchLPCA[1]),
         col="black")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  ##Heat Map Scale
  points(seq(1,2,1/10),rep(-2,11),col=colEfficiency,pch=15,cex=3)
  text(1,-1.8,col="black",labels= prettyNum(min(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  text(1+0.5,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency)/2,digits=1,format="f" ),cex=cex)
  text(2,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  max(mergedCapDat_filt$Efficiency)/2
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  title(NA,sub=paste("% of Efficiency Variance Explained : ",prettyNum( 100*pcEffVar/EffVar) ))
  


dev.off()

###DF Specific Text Plot Text tuninng
pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC3_5_EfficiencyColour_DF.pdf",sep=""),width=7,height=7)
  
  plot(rawd[,1], rawd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID) ],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchL[4+datpolyFactor_norm$Group],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+datpolyFactor_norm$Group],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC3"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC5" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.8*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.5*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lwd=2,lty=3)
  text(1.5*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  text(-0.6*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.7*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(1.2*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.8*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.1*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  legend("bottomleft",legend=c("LF","NF","DF"),pch=c(pchLPCA[2],pchLPCA[3],pchLPCA[1]),
         col="black")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  ##Heat Map Scale
  points(seq(1,2,1/10),rep(-2,11),col=colEfficiency,pch=15,cex=3)
  text(1,-1.8,col="black",labels= prettyNum(min(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  text(1+0.5,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency)/2,digits=1,format="f" ),cex=cex)
  text(2,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  max(mergedCapDat_filt$Efficiency)/2
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  title(NA,sub=paste("% of Efficiency Variance Explained : ",prettyNum( 100*pcEffVar/EffVar) ))
  


dev.off()


pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC2_3_EfficiencyColour_NF.pdf",sep=""),width=7,height=7)
  
  plot(rawd[,1], rawd[,2],
       #col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID) ],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchL[4+datpolyFactor_norm$Group],
       #col=colourLegL[datpolyFactor_norm$Group], pch=pchL[4+datpolyFactor_norm$Group],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC2"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC3" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.8*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.5+1.5*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lwd=2,lty=3)
  text(1.5*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  text(-0.6*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-1.2*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(-0.2*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.8*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.1*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  legend("bottomleft",legend=c("LF","NF","DF"),pch=c(pchLPCA[2],pchLPCA[3],pchLPCA[1]),
         col="black")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  ##Heat Map Scale
  points(seq(1,2,1/10),rep(-2,11),col=colEfficiency,pch=15,cex=3)
  text(1,-1.8,col="black",labels= prettyNum(min(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  text(1+0.5,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency)/2,digits=1,format="f" ),cex=cex)
  text(2,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  max(mergedCapDat_filt$Efficiency)/2
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  title(NA,sub=paste("% of Efficiency Variance Explained : ",prettyNum( 100*pcEffVar/EffVar) ))
  


dev.off()

##%*%pca_norm$rotation[1,]
sum(pca_norm$rotation[1,]^2) #Unit vectors

### Calc the Projection of Efficiency on the variance captured by each  PC
s = 0
tIdx <- 3 ##Choose Which Variable to Calc The  COntribution to Variance from Each PC 
Ve <- vector()
for (i in 1:9)
{
  ##Total Efficiency Vairance 
  s = s + (pca_norm$sdev[i]^2)*sum( ( pca_norm$rotation[,i]*pca_norm$rotation[tIdx,])^2  ) ### * pca_norm$sdev[i
  ##Efficiency COntribution to Var Of Each PC 
  Ve[i] <- (pca_norm$sdev[i]^2)*sum( ( pca_norm$rotation[,i]*pca_norm$rotation[tIdx,])^2  )##sum( (pca_norm$rotation[,i]%*%pca_norm$rotation[1,])^2 * pca_norm$sdev[i])^2
}
##Now Choose a PC

##CHeck Contribution Of PC to Efficiency Variance
sum(Ve/s)
###PLot Relative COntrib To Variance 
plot((100*Ve/s) ) 



pdf(file= paste(strPlotExportPath,"/stat/stat_PCAHuntVariablesAndEfficiencyPC1_2_GroupColour_ALL.pdf",sep=""),width=7,height=7)
  ## bottom, left,top, right
  par(mar = c(4.2,4.3,1,1))
  
  plot(rawd[,1], rawd[,2],
       ######col=colClass[1+as.numeric(mergedCapDat$Undershoot > 1)], pch=pchL[4+datpolyFactor_norm$Group], 
       #col=colEfficiency[round(mergedCapDat_filt$Efficiency*10)], pch=pchL[4+as.numeric(mergedCapDat_filt$groupID) ],
       #col=colClass[as.numeric(mergedCapDat_filt$Cluster)], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       col=colourGroup[mergedCapDat_filt$groupID ], pch=pchLPCA[as.numeric(mergedCapDat_filt$groupID)],
       #xlab="PC1",ylab="PC2",
       xlim=c(-2.5,2.5),ylim=c(-2,2.5),
       xlab=NA,ylab=NA,
       cex=cex,cex.axis=cex ) #xlim=c(-4,4),ylim=c(-4,4)
  mtext(side = 1,cex=cex, line = lineXAxis,  "PC1"   ,cex.main=cex )
  mtext(side = 2,cex=cex, line = lineAxis, "PC2" ,cex.main=cex)
  
  
  scaleV <- 2
  ##Distance to Prey Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],lwd=3)
  text(0.7*scaleV*pca_norm$rotation[2,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.2*scaleV*pca_norm$rotation[2,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[1],labels="Distance")
  ##CaptureSpeed  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],lwd=2,lty=3)
  text(0.8*scaleV*pca_norm$rotation[3,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+0.8*scaleV*pca_norm$rotation[3,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[2],labels="Speed")
  
  ##Undershoot Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2)
  #  text(0.4*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.1*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,labels="Overshoot")
  arrows(0,0,-scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",lty=2,lwd=2)
  text(-1.5*scaleV*pca_norm$rotation[4,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,-1.0*scaleV*pca_norm$rotation[4,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="black",labels="Undershoot")
  
  ##TimeToHit Prey Prod Axis  Component Projection
  arrows(0,0,scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],lty=5,lwd=2)
  text(0.8*scaleV*pca_norm$rotation[6,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,1.2*scaleV*pca_norm$rotation[6,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col=colFactrAxes[6],labels="t Prey")
  
  ##DistXSpeed Prod Axis  Component Projection
  #arrows(0,0,scaleV*pca_norm$rotation[5,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleV*pca_norm$rotation[5,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="purple",lty=5)
  
  ##EFFICIENCY Prod Axis  Component Projection
  scaleVE <- scaleV
  arrows(0,0,scaleVE*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,scaleVE*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",lty=2,lwd=2)
  text(0.4*scaleV*pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]]^2,0.1+1.0*scaleV*pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]]^2,col="blue",labels="Efficiency")
  
  # ##Heat Map Scale
  # points(seq(1,2,1/10),rep(-2,11),col=colEfficiency,pch=15,cex=3)
  # text(1,-1.8,col="black",labels= prettyNum(min(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  # text(1+0.5,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency)/2,digits=1,format="f" ),cex=cex)
  # text(2,-1.8,col="black",labels= prettyNum(max(mergedCapDat_filt$Efficiency),digits=1,format="f" ),cex=cex)
  # max(mergedCapDat_filt$Efficiency)/2
  # 
  
  legend("bottomleft",legend=c( paste0(groupLabels[3],' #',table(mergedCapDat_filt$groupID)[3]),
                                paste0(groupLabels[2],' #',table(mergedCapDat_filt$groupID)[2]),
                                paste0(groupLabels[1],' #',table(mergedCapDat_filt$groupID)[1]) 
  ),
  pch=c(pchLPCA[3],pchLPCA[2],pchLPCA[1]),
  col=c(colourGroup[3],colourGroup[2],colourGroup[1]) )## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  #legend("bottomright",legend=c("Slow","Fast"),fill=colClass, col=colClass,title="Cluster")## c(colourLegL[2],colourLegL[3],colourLegL[1])) # c(colourH[3],colourH[2])
  
  #Percentage of Efficiency Variance Explained
  nComp <- length(pca_norm$sdev)
  pcEffVar <- ((pca_norm$rotation[1,][pcAxis[1]]*pca_norm$sdev[pcAxis[1]])^2 + (pca_norm$rotation[1,][pcAxis[2]]*pca_norm$sdev[pcAxis[2]])^2)
  EffVar <- sum((pca_norm$rotation[1,][1:nComp]*pca_norm$sdev[1:nComp])^2)
  
  title(NA,sub=paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(mergedCapDat_filt$Efficiency)/mean(mergedCapDat_filt$Efficiency) ,digits=2)) )
  message("Captured Variance ",prettyNum( (pca_norm$sdev[pcAxis[1]]^2 + pca_norm$sdev[pcAxis[2]]^2) /sum( pca_norm$sdev ^2),digits=3,format="f" ) )
  message(paste(" Efficiency variance captured: ",prettyNum( 100*pcEffVar/EffVar,digits=3), " Coeff. variation:",prettyNum(sd(mergedCapDat_filt$Efficiency)/mean(mergedCapDat_filt$Efficiency) ,digits=2)))

dev.off()






biplot(pca_norm,choices=c(1,2))


theta <- function (a,b){ return( (180/pi)   *acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )) }

pcAxis <- c(3,5)
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[1,pcAxis]) #Efficiency
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[2,pcAxis]) #Group
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[3,pcAxis]) #DistanceToPrey
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[4,pcAxis]) #CaptureSpeed_norm
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[5,pcAxis]) #Undershoot_norm
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[6,pcAxis]) #DistSpeedProd
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[7,pcAxis]) #DistUnderProd
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[8,pcAxis]) #SpeedUnderProd
theta(pca_norm$rotation[1,pcAxis], pca_norm$rotation[9,pcAxis]) #All
##PCA

theta(pca_norm$rotation[,1], pca_norm$rotation[,2])


library(rgl)

open3d()##mergedCapDat$groupID
rgl::plot3d( x=rawd[,1], z=rawd[,2], y=rawd[,3], col = colourLegL[datpolyFactor_norm$Group] , type = "s", radius = 0.5,
             xlab="PC1", zlab="PC2",ylab="PC3",
             xlim=c(-8.,8), ylim=c(-8,8), zlim=c(-8,8),
             box = FALSE ,aspect = TRUE
             #,expand = 1.5
)
###END PCA PLOT ##



Ei_LF_norm=eigen(cov(datpolyFactor_norm))
symnum(cov(datpolyFactor_norm))
Ei_LF_norm
## ##Make MAtrix
datpolyFactor <- with(mergedCapDat[mergedCapDat$groupID == 'NL',],{
  cbind(Efficiency, #1
        #HuntPower, # ## Does not CoVary With Anyhting 
        DistanceToPrey, #2
        CaptureSpeed, #3
        Undershoot, #4
        DistanceToPrey*CaptureSpeed, #5
        DistanceToPrey*Undershoot, #6
        CaptureSpeed*Undershoot, #7
        DistanceToPrey*CaptureSpeed*Undershoot #8
  )
  
})

Ei_NL<-eigen(cov(datpolyFactor))
Ei_NL

datpolyFactor <- with(mergedCapDat[mergedCapDat$groupID == 'LL',],{
  cbind(Efficiency, #1
        #HuntPower, # ## Does not CoVary With Anyhting 
        DistanceToPrey, #2
        CaptureSpeed, #3
        Undershoot, #4
        DistanceToPrey*CaptureSpeed, #5
        DistanceToPrey*Undershoot, #6
        CaptureSpeed*Undershoot, #7
        DistanceToPrey*CaptureSpeed*Undershoot #8
  )
})


##Ei_LF=eigen(cov(datpolyFactor))
pca_LL <- prcomp(datpolyFactor,scale.=TRUE)
summary(pca_LL)
biplot(pca_LL)

Ei_LF
plot(pca_LL)

## ##Make MAtrix
datpolyFactor <- with(mergedCapDat[mergedCapDat$groupID == 'LL',],{
  cbind(Efficiency, #1
        #HuntPower, # ## Does not CoVary With Anyhting 
        DistanceToPrey, #2
        CaptureSpeed, #3
        Undershoot, #4
        DistanceToPrey*CaptureSpeed, #5
        DistanceToPrey*Undershoot, #6
        CaptureSpeed*Undershoot, #7
        DistanceToPrey*CaptureSpeed*Undershoot #8
  )
  
})

Ei_DL <- eigen(cov(datpolyFactor))

### Print Eugen MAtrix
Ei_LL
Ei_NL
Ei_DL

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x=cov(datpolyFactor), col=col,symm = F)

#library(corrplot)

lm(Efficiency ~ (DistanceToPrey+CaptureSpeed+Undershoot)^3,data=mergedCapDat[mergedCapDat$groupID == 'NL',])




huntPowerColour <- rfc(8)
open3d()##mergedCapDat$groupID
rgl::plot3d( x=mergedCapDat$CaptureSpeed, z=mergedCapDat$DistanceToPrey, y=mergedCapDat$HuntPower, col = huntPowerColour[round(mergedCapDat$HuntPower)] , type = "s", radius = 1.3,
             xlab="Capture Speed (mm/sec)", zlab="Hunt Power",ylab="Distance to prey (mm)",
             xlim=c(0.,80), ylim=c(0,0.6), zlim=c(0,10),
             box = FALSE ,aspect = TRUE
             #,expand = 1.5
)

open3d()##mergedCapDat$groupID
rgl::plot3d( x=datMergedCapAndSuccess_LF$CaptureSpeed, z=datMergedCapAndSuccess_LF$DistanceToPrey, y=datMergedCapAndSuccess_LF$HuntPower, col = huntPowerColour[round(mergedCapDat$HuntPower)] , type = "s", radius = 1.3,
             xlab="Capture Speed (mm/sec)", ylab="Hunt Power",zlab="Distance to prey (mm)",
             xlim=c(0.,80), zlim=c(0,0.6), ylim=c(0,10),
             box = FALSE ,aspect = TRUE
             #,expand = 1.5
)


rgl::rgl.viewpoint(60,10)

