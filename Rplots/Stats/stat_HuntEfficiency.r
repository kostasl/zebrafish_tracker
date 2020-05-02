## stat Model Success Vs Failure - Plots Distribution For HuntRate and Success Rate
# Makes A model Assuming Poisson Rate lambda and Hunt Success Probability q
# Note: We can exclude the fish that produced no events when excluding eventID == 0- 
## UPDATE 13/11/19 : Converted from general hunt event rate estimation - to capture attempt estimation
## ** Efficiency plot reports attempts at target vs Success

source("TrackerDataFilesImport_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")
source("plotHuntStat_lib.r")

## The mixture Of Poisson Drawing from Gammaa, gives a negative binomial
## ID identifies the group ID of each larva
modelNBinom="model {
         
         for(i in 1:3) {
             r[i] ~ dgamma(1,1) ##
             q[i] ~ dunif(0.0,1)

             f[i] ~ dbeta(1,1) ##Prob of capture Fail
             p[i] ~ dbeta(1,1) ##Prob of capture success
             t[i] ~ dbeta(1,1) ##Prob Of Enganging With Prey Given HuntMode Is On
         }

         for(i in 1:NTOT){
             Events[i] ~  dnegbin(q[ID[i]],r[ID[i]] )
             TrackPrey[i] ~ dbinom(t[ID[i]],Events[i])
             Success[i] ~ dbinom(p[ID[i]],TrackPrey[i])
             Fail[i] ~ dbinom(f[ID[i]],TrackPrey[i])
         }
}"




## Model (Wrongly) Assumes Single Distribution Prior For Rates - Yet, fitting the group event rates shows
## that it is best to assume a mixture of hunt rates exist in each group, such that the group hunt rates appears neg Binom.
## wiki suggest gamma prior for lambda
##NOT USED###
# modelPoisson="model {
#          
#          for(i in 1:3) {
#              lambda[i] ~ dgamma(1,1) ##Suggested in wiki 
#              #lambda[i] ~ dexp(1) 
#              q[i] ~ dbeta(1,1) 
#              p[i] ~ dbeta(1,1) ##Prob of capture success
#              t[i] ~ dbeta(1,1) ##Prob Of Enganging With Prey Given HuntMode Is On
#          }
# 
#          for(i in 1:NTOT){
#              Events[i] ~ dpois(lambda[ID[i]])
#              TrackPrey[i] ~ dbinom(t[ID[i]],Events[i])
#              Success[i] ~ dbinom(q[ID[i]],TrackPrey[i])
#              Fail[i] ~ dbinom(p[ID[i]],TrackPrey[i])
#              
#          }
# }"


#strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
#strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged2"
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
#datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
##These Are Double/2nd Trials on LL, or Simply LL unpaired to any LE (Was checking Rates)
#vxCludeExpID <- c(4421,4611,4541,4351,4481,4501,4411)

##We Can Choose To Exclude The Fish That Produced No Hunting Events
#datHuntLabelledEventsSB <- datHuntLabelledEventsSB[  !(datHuntLabelledEventsSB$expID %in% vxCludeExpID) & 
#                                                     datHuntLabelledEventsSB$groupID %in% c("LL","NL","DL") ,]
##Load From Central Function
datHuntLabelledEventsSB <- getLabelledHuntEventsSet()

#saveRDS(datHuntLabelledEventsSB,file="pub/dat/datHuntLabelledEvents.rds")

datHuntLabelledEventsSB_LIVE <-  datHuntLabelledEventsSB[datHuntLabelledEventsSB$groupID %in% c("LL","DL","NL"),]
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB_LIVE)

tblResSB <- table(convertToScoreLabel(datHuntLabelledEventsSB_LIVE$huntScore),datHuntLabelledEventsSB_LIVE$groupID)

datFishSuccessRate$groupID <- factor(datFishSuccessRate$groupID)
strGroups <-levels(datFishSuccessRate$groupID)

NRecCount_DL <- table(datFishSuccessRate$groupID)["DL"]
NRecCount_NL <- table(datFishSuccessRate$groupID)["NL"]
NRecCount_LL <- table(datFishSuccessRate$groupID)["LL"]


datatest=list(Success=datFishSuccessRate$Success,
              Fail=datFishSuccessRate$Fails,
              TrackPrey=datFishSuccessRate$Fails+datFishSuccessRate$Success, ##Number of Prey Engangements
              Events=datFishSuccessRate$CaptureEvents, ######Includes No_Target - ie cases where stimuli Trigger Fish HuntMode But no Prey Tracking Seems to take place
              ID=as.numeric(datFishSuccessRate$groupID),
              NTOT=nrow(datFishSuccessRate));

varnames1=c("q","p","t","f","r")
burn_in=1000;
steps=20000;
thin=10;
chains=3

library(rjags)
strModelName = "model1.tmp"
fileConn=file(strModelName)
writeLines(modelNBinom,fileConn);
close(fileConn)

m=jags.model(file=strModelName,data=datatest,n.chains=chains);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames1)

## Save Samples 
save(draw,file =paste(strDataExportDir,"stat_huntefficiencyModel_RJags.RData",sep=""))
#######PLOT RESULTS
#X11()
colourH <- c("#0303E663","#03B30363","#E6030363")
colourD <- c("#0303E623","#03B30323","#E6030323")
colourL <- c("#0303E6AF","#03B303AF","#E60303AF")
#for(i in 1:3) hist(draw$q[i,,1],breaks=seq(0,0.5,0.01),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

#X11()
#for(i in 1:3) hist(draw$lambda[i,,1],breaks=seq(5,20,0.2),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

##Density in 2D of Success Vs Fail
#X11()


### Draw Distribution oF Hunt Rates - 
## for the exp draw (z= p/(1-p)) ## But it is the same for Rate Of Gamma Too / Or inverse for scale
plotsamples <- 1500
schain <-1:3
Range_ylim <- c(1,25)
cex <- 1.2

### It looks like Tha Gamma scale is the inverse:  gamma rate!
##13/11/19 Spotted difference in GammaRate q/(1-q) was inverse - now matches HuntRateInPreyRange model analysis
HEventHuntGammaRate_LL <-tail(draw$q[2,,schain],plotsamples)/(1-tail(draw$q[2,,schain],plotsamples));
HEventHuntGammaRate_DL <- tail(draw$q[1,,schain],plotsamples)/(1-tail(draw$q[1,,schain],plotsamples));      
HEventHuntGammaRate_NL <- tail(draw$q[3,,schain],plotsamples)/(1-tail(draw$q[3,,schain],plotsamples));      
HEventHuntGammaShape_LL <- tail(draw$r[2,,schain],plotsamples)
HEventHuntGammaShape_DL <- tail(draw$r[1,,schain],plotsamples)
HEventHuntGammaShape_NL <- tail(draw$r[3,,schain],plotsamples)

HEventSuccess_LL <-tail(draw$p[2,,schain],plotsamples)
HEventSuccess_DL <-tail(draw$p[1,,schain],plotsamples)
HEventSuccess_NL <-tail(draw$p[3,,schain],plotsamples)

###Mean Rates As Exp OF Gamma
MeanHuntRate_LL <- HEventHuntGammaShape_LL*1/HEventHuntGammaRate_LL
MeanHuntRate_DL <- HEventHuntGammaShape_DL*1/HEventHuntGammaRate_DL
MeanHuntRate_NL <- HEventHuntGammaShape_NL*1/HEventHuntGammaRate_NL


HConsumptionRate_NL <- HEventSuccess_NL*MeanHuntRate_NL
HConsumptionRate_LL <- HEventSuccess_LL*MeanHuntRate_LL
HConsumptionRate_DL <- HEventSuccess_DL*MeanHuntRate_DL

#'' Measure correlation between success probability and efficiency - 
chain <-1 ##use chain 1
HCovRateAndEfficiency_NL <- cor( tail(HEventSuccess_NL[,chain],plotsamples),tail(MeanHuntRate_NL[,chain],plotsamples) )
HCovRateAndEfficiency_LL <- cor( tail(HEventSuccess_LL[,chain],plotsamples),tail(MeanHuntRate_LL[,chain],plotsamples) )
HCovRateAndEfficiency_DL <- cor( tail(HEventSuccess_DL[,chain],plotsamples),tail(MeanHuntRate_DL[,chain],plotsamples) )


### Compare Rates of Capture Attempt (hunt rate~) difference 
message("Mean Capture Attempts LL:",prettyNum(mean(MeanHuntRate_LL),digits=3), " DL:",prettyNum(mean(MeanHuntRate_DL),digits=3), " NL:",prettyNum(mean(MeanHuntRate_NL),digits=3 ) ) 
message("Compare mean Hunt Rate LL-NL:", (mean(MeanHuntRate_LL)-mean(MeanHuntRate_NL) )  )
### Compare Hunt Efficiency difference 
HEventSuccess_LLVsNL <- HEventSuccess_LL - HEventSuccess_NL
HEventSuccess_LLVsDL <- HEventSuccess_LL - HEventSuccess_DL
HEventSuccess_NLVsDL <- HEventSuccess_NL - HEventSuccess_DL
message("Compare mean Hunt Efficiency LL-NL:", (mean(HEventSuccess_LLVsNL) )  )
message("Compare mean Hunt Efficiency LL-DL:", (mean(HEventSuccess_LLVsDL) )  )

P_LLgtNL <- length(HEventSuccess_LLVsNL[HEventSuccess_LLVsNL >  0])/length(HEventSuccess_LLVsNL)
P_LLgtDL <- length(HEventSuccess_LLVsDL[HEventSuccess_LLVsDL >  0])/length(HEventSuccess_LLVsDL)
P_NLgtDL <- length(HEventSuccess_NLVsDL[HEventSuccess_NLVsDL >  0])/length(HEventSuccess_NLVsDL)
message("Prob that LF group is More efficient than NF:",P_LLgtNL, " and DF:", P_LLgtDL)
message("Prob that NF is More Efficient than DF:",prettyNum(P_NLgtDL),digits=3 ) 
## Compare Capture Attempts
### Mean Rates As Exp OF Gamma
MeanHuntRate_LLvsNL <- MeanHuntRate_LL - MeanHuntRate_NL 
MeanHuntRate_LLvsDL <- MeanHuntRate_LL - MeanHuntRate_DL 
MeanHuntRate_NLvsDL <- MeanHuntRate_NL - MeanHuntRate_DL 
P_LLgtNL <- length(MeanHuntRate_LLvsNL[MeanHuntRate_LLvsNL >  0])/length(MeanHuntRate_LLvsNL)
P_LLgtDL <- length(MeanHuntRate_LLvsDL[MeanHuntRate_LLvsDL >  0])/length(MeanHuntRate_LLvsDL)
P_NLgtDL <- length(MeanHuntRate_NLvsDL[MeanHuntRate_NLvsDL >  0])/length(MeanHuntRate_NLvsDL)
message("Prob that LF group Attempts More than NF:",prettyNum(P_LLgtNL,3), " and DF:", prettyNum(P_LLgtDL,3), " NF>DF",prettyNum(P_NLgtDL,3) )

HConsumptionRate_NL <- HEventSuccess_NL*MeanHuntRate_NL
HConsumptionRate_LL <- HEventSuccess_LL*MeanHuntRate_LL
HConsumptionRate_DL <- HEventSuccess_DL*MeanHuntRate_DL



###MAIN OUTPUT PLOT ##

load(file =paste(strDataExportDir,"stat_huntefficiencyModel_RJags.RData",sep=""))

#strPlotName = paste(strPlotExportPath,"/stat/fig3-stat_HuntRateAndEfficiencyEstimationNegBin_Success.pdf",sep="")

strPlotName = paste(strPlotExportPath,"/stat/fig3-stat_HuntSuccessPieChart.pdf",sep="")
pdf(strPlotName,width=14,height=4.7,
    title="Hunting Success Labelled results ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
  #svg(filename=strPlotName,width=8,height=8)
  outer = FALSE
  line <- 2.6 ## SubFig Label Params
  lineGroupLabel <- line - 32 ##pie chart group label
  cex = 1.4
  adj  = 0.5
  padj <- -0
  las <- 1
  
  #layout(matrix(c(1,2,3,4,4,4,5,6,6), 3, 3, byrow = TRUE))
  layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(3.9,4.3,5.5,1))
  #colourL <-  c("#66C2A5","#B3B3B3") Fails Colourblind
  colPaired <- rev(brewer.pal(4,'Paired')) ## chosen for colorblindness
  colourL <- c(colPaired[1],colPaired[3]) 
  
  ## Get Number Of Larvae / 
  nlNL <- NROW(table(datHuntLabelledEventsSB[datHuntLabelledEventsSB$groupID == "NL",]$expID))
  nNL <- pieChartLabelledSuccessVsFails(tblResSB,"NL",colourL) #c(colourLegL[1],colourL[2]) 
  mtext(c(expression(),  bquote("NF ")),
        at="bottom",  outer=F,side=3,col="black",font=2,las=las,line=lineGroupLabel,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote(.(nNL) ~ "Hunt events ")),
        at="top",  outer=F,side=3,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote( .(nlNL) ~ "Larvae" )  ),
        at="top", outer=F,side=3,col="black",font=2,las=las,line=line-2,padj=padj,adj=adj,cex=cex)
  #mtext("A",at="topleft",outer=F,side=2,col="black",font=2,las=las,line=3,padj=-11,adj=0,cex=cex,cex.main=4)
  
  legend("bottomright",legend=c("Success","Fail"),
         fill=colourL, #c(colourLegL[1],colourL[2]),
         col = colourL,
         bg = "white",cex=cex+0.2,
         merge=FALSE,horiz=FALSE)
  
  
  ## Get Number Of Larvae
  nlLL <- NROW(table(datHuntLabelledEventsSB[datHuntLabelledEventsSB$groupID == "LL",]$expID))
  nLL <- pieChartLabelledSuccessVsFails(tblResSB,"LL",colourL)
  mtext(c(expression(),  bquote("LF ")),
        at="bottom",  outer=F,side=3,col="black",font=2,las=las,line=lineGroupLabel,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote(.(nLL) ~ "Hunt events ")),
        at="top",  outer=F,side=3,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote( .(nlLL) ~ "Larvae" )  ),
        at="top", outer=F,side=3,col="black",font=2,las=las,line=line-2,padj=padj,adj=adj,cex=cex)
  
  ## Get Number Of Larvae / 
  nlDL <- NROW(table(datHuntLabelledEventsSB[datHuntLabelledEventsSB$groupID == "DL",]$expID))
  ##Returns Number of Hunt Events
  nDL <- pieChartLabelledSuccessVsFails(tblResSB,"DL",colourL) #pieChartLabelledEvents(tblResSB,"DL")
  mtext(c(expression(),  bquote("DF ")),
        at="bottom",  outer=F,side=3,col="black",font=2,las=las,line=lineGroupLabel,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote(.(nDL) ~ "Hunt events ")),
        at="top",  outer=F,side=3,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex=cex)
  mtext(c(expression(),  bquote( .(nlDL) ~ "Larvae" )  ),
        at="top", outer=F,side=3,col="black",font=2,las=las,line=line-2,padj=padj,adj=adj,cex=cex)
  #text(x=1.4,y=-0.8,labels = "SB",cex=1.5)  

dev.off()


plotWidthIn <- 8
strPlotName = paste(strPlotExportPath,"/stat/fig3-stat_ModelCaptureRateAndEfficiency.pdf",sep="")
pdf(strPlotName,width=plotWidthIn,height=7,
    title="Hunting Success Baysian Estimation ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
  outer = FALSE
  line <- 2.6 ## SubFig Label Params
  lineGroupLabel <- line - 32 ##pie chart group label
  cex = 1.4
  adj  = 0.5
  padj <- -0
  las <- 1

  ##Margin: (Bottom,Left,Top,Right )
  #par(mar = c(5,6,2,3))
  par(mar = c(4.2,4.8,1.1,1))
  #layout(matrix(c(1,,2), 1, 2, byrow = TRUE))
  
####### Efficiency Inference Plot ## Taken From stat_SyccessVsFailModel.r ####
  nlevels <- 5
  zLL <- kde2d(c(HEventSuccess_LL[,schain]), c(MeanHuntRate_LL[,schain]),n=80)
  zNL <-  kde2d(c(HEventSuccess_NL[,schain]), c(MeanHuntRate_NL[,schain]),n=80)
  zDL <-  kde2d(c(HEventSuccess_DL[,schain]), c(MeanHuntRate_DL[,schain]),n=80)
  
  plot(HEventSuccess_DL, MeanHuntRate_DL,col=colourHPoint[3],ylim=Range_ylim,xlim=c(0,0.5),pch=pchL[3],
       main=NA, #"Bayesian Estimation for Hunt Rate and Efficiency",
       xlab=NA,#"Probability of Success q",
       ylab=NA,cex.main =cex,cex.axis=1.5 )#(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
  points(HEventSuccess_LL, MeanHuntRate_LL,col=colourHPoint[2],ylim=Range_ylim,xlim=c(0.1,0.5),pch=pchL[2])
  points(HEventSuccess_NL, MeanHuntRate_NL,col=colourHPoint[1],ylim=Range_ylim,xlim=c(0.1,0.5),pch=pchL[1])
  mtext(side = 1, cex=cex, line = line,expression(paste("Probability of success (",q,")" ) )  ) 
  mtext(side = 2, cex=cex, line = line, expression(paste("Estimated capture attempts / 10min" ) )  ) #(",lambda,")"
  #mtext("B",at="topleft",outer=F,side=2,col="black",font=2,las=las,line=4,padj=-11,adj=0,cex=cex,cex.main=4)
  
  contour(zDL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  contour(zLL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  
  contour(zNL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  #legend("topright", legend=paste(strGroups," n=",c(NRecCount_DL,NRecCount_LL,NRecCount_NL)),fill=colourL)
  legend("topright",cex=cex+0.2,
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NRecCount_NL)  ),
                    bquote(LF[""] ~ '#' ~ .(NRecCount_LL)  ),
                    bquote(DF[""] ~ '#' ~ .(NRecCount_DL)  )  ), #paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         pch=pchL, col=colourLegL)

dev.off()  

strPlotName = paste(strPlotExportPath,"/stat/fig3-stat_ModelConsumption.pdf",sep="")
pdf(strPlotName,width=plotWidthIn,height=7,
    title="Estimated Consumption  ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
  #par(mar = c(5,6,2,3))
  par(mar = c(4.2,4.7,1.1,1))
  #### Consumption
  plot(density(HConsumptionRate_NL),xlim=c(0,8),ylim=c(0,1.5),col=colourLegL[1],lwd=4,lty=1
       ,cex.main =cex,cex.axis=1.5, xlab=NA,ylab=NA,main=NA) #"Mean consumption per larva"
  lines(density(HConsumptionRate_LL),col=colourLegL[2],lwd=4,lty=2)
  lines(density(HConsumptionRate_DL),col=colourLegL[3],lwd=4,lty=3)
  legend("topright",cex=cex+0.2,
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NRecCount_NL)  ),
                    bquote(LF[""] ~ '#' ~ .(NRecCount_LL)  ),
                    bquote(DF[""] ~ '#' ~ .(NRecCount_DL)  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
         col=colourLegL,lty=c(1,2,3),lwd=4)
  mtext(side = 1, cex=cex, line = line, expression(paste("Estimated consumption (Prey/10min)  ") ))
  mtext(side = 2, cex=cex, line = line, expression("Density function") )
#  mtext("C",at="topleft",outer=F,side=2,col="black",font=2,las=las,line=4,padj=-11,adj=0,cex=cex,cex.main=4)
dev.off()


strPlotName = paste(strPlotExportPath,"/stat/fig3-stat_ecdf_HuntPower.pdf",sep="")
pdf(strPlotName,width=plotWidthIn,height=7,
    title="Hunting Success Baysian Estimation ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
  par(mar = c(4.2,4.7,1.1,1))
  #### Plot Hunt Power ####
  plotHuntPowerDataCDF(datHuntLabelledEventsSB)
  #mtext("D",at="topleft",outer=F,side=2,col="black",font=2,las=las,line=4,padj=-11,adj=0,cex=cex,cex.main=4)
  
dev.off()
#embed_fonts(strPlotName)


strPlotName = paste(strPlotExportPath,"/stat/fig3_stat_Efficiency_CDF.pdf",sep="")
pdf(strPlotName,width=plotWidthIn,height=7, #14*2/3
    title="Hunting Rate and Efficiency - Labelled results and Bayesian Inference ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
  
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.2,4.7,1.1,1))
  plotHuntEfficiencyDataCDF(datHuntLabelledEventsSB)

dev.off()
#### ##### # ## ## # # # # 





### Plot Covariance of Hunt Rate To Prob Of Success
fNL <- density((HEventSuccess_NL[,1]*MeanHuntRate_NL[,1]))
fLL <- density(HEventSuccess_LL[,1]*MeanHuntRate_LL[,1])
plot(density((HEventSuccess_LL[,1]*MeanHuntRate_LL[,1]) ),col=colourLegL[2],ylim=c(0,1))
lines(density((HEventSuccess_LL[,1]*MeanHuntRate_LL[,1])),col=colourLegL[2],lwd=2)
lines(density((HEventSuccess_NL[,1]*MeanHuntRate_NL[,1])),col=colourLegL[1])
lines(density((HEventSuccess_DL[,1]*MeanHuntRate_NL[,1])),col=colourLegL[3])

plot(density((HEventSuccess_NL[,1]*MeanHuntRate_LL[,1])),col=colourLegL[2])
lines(density((HEventSuccess_NL[,1]*MeanHuntRate_NL[,1])),col=colourLegL[1])
lines(density((HEventSuccess_NL[,1]*MeanHuntRate_DL[,1])),col=colourLegL[3])



strPlotName = paste(strPlotExportPath,"/stat/fig3S3_stat_HuntRateAndEfficiency_PDF.pdf",sep="")
pdf(strPlotName,width=7,height=7,
    title="Hunting Rate and Efficiency - Labelled results and Bayesian Inference ", ##on distribution of hunt rate parameter and probability of success, based on labelled data set
    onefile = TRUE,compress=FALSE) #col=(as.integer(filtereddatAllFrames$expID))
par(mar = c(4.2,4.7,1.1,1))

plotHuntEfficiencyDataPDF(datHuntLabelledEventsSB)
dev.off()
# 
# ######## ## # # # ## 
# ########################
# strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation_Fails.pdf",sep="")
# pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of Engaging with Prey and Failing, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
# plot(draw$p[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=Range_ylim,xlim=c(0,1),pch=19,
#      main="Bayesian Estimation for Hunt Rate and Efficiency (Fails)",
#      xlab="Probability of Failing p",
#      ylab=(expression(paste("Hunt Rate (",lambda,")" ) ) )  ) #paste("Hunt Rate", )
# points(draw$p[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=Range_ylim,xlim=c(0.1,1),pch=19)
# points(draw$p[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=Range_ylim,xlim=c(0.1,1),pch=19)
# legend("topright", legend=strGroups,fill=colourL)
# dev.off()
# 
# strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation_Track.pdf",sep="")
# pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of Engaging with Prey, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
# plot(draw$t[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=c(5,26),xlim=c(0,1),pch=19,
#      main="Bayesian Estimation for Hunt Rate and Tracking Efficiency / Engaging with prey ",
#      xlab="Probability of Tracking Prey t",
#      ylab=(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
# points(draw$t[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=c(5,26),xlim=c(0,1),pch=19)
# points(draw$t[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=c(5,26),xlim=c(0,1),pch=19)
# legend("topright", legend=strGroups,fill=colourL)
# dev.off()
# 
# 
# 
# ## Show  consumption estimates per fish , mean consumption ##
# strPlotName = paste(strPlotExportPath,"/stat/stat_HuntConsumptionEstimation.pdf",sep="")
# pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of Engaging with Prey, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
# 
# plot(density(HConsumptionRate_NL),xlim=c(0,8),ylim=c(0,1.5),col=colourLegL[1],lwd=3.5,lty=1,xlab=NA,ylab=NA,main=NA) #"Mean consumption per larva"
# lines(density(HConsumptionRate_LL),col=colourLegL[2],lwd=3.5,lty=2)
# lines(density(HConsumptionRate_DL),col=colourLegL[3],lwd=3.5,lty=3)
# legend("topright",
#        legend=c(  expression (),
#                   bquote(NF["e"] ~ '#' ~ .(NRecCount_NL)  ),
#                   bquote(LF["e"] ~ '#' ~ .(NRecCount_LL)  ),
#                   bquote(DF["e"] ~ '#' ~ .(NRecCount_DL)  )  ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
#        col=colourLegL,lty=c(1,2,3),lwd=3)
# mtext(side = 1,cex=1.3, line = 2.2, expression(paste("Estimated consumption (Prey/10min)  ") ))
# mtext(side = 2,cex=1.3, line = 2.2, expression("Density function") )
# 
# dev.off()
# #hist(draw$q[1,,1],breaks=seq(0,30,length=100),col=colourH[1])
# #hist(drawNL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[2])
# #hist(drawDL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[3])
# #legend(2,500,legend = c("LL","NL","DL"),fill=c(colourH[1],colourH[2],colourH[3]))
# 
# #hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
# #hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
# #hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

