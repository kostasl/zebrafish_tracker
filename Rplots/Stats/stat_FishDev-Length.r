## kostasl 11 Jun 2019
## stat Model for comparing larval development across the 3 rearing groups
## we assube 3 gaussians, one for each group, and  infer mean and std dev for each

source("config_lib.R")
source("TrackerDataFilesImport_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")


## The mixture Of Poisson Drawing from Gammaa, gives a negative binomial
modelNorm="model {
         
         for(i in 1:3) {
             mu[i]   ~ dunif(3,6 ) ## Mean length in mm
             prec[i] ~ dunif(0,50) ##dgamma(1,1) ##
             sigma[i] <- sqrt(1/prec[i]) 
        }

         for(i in 1:NTOT){
             LarvaLength[i] ~  dnorm(mu[groupID[i] ],prec[groupID[i]] )

         }
}"

##New Fish Size Record ###
#datFlatPxLength <- data.frame(
#  rbind(cbind(datFishLength$LF,which(strGroupID == "LL")),
#        cbind(datFishLength$NF,which(strGroupID == "NL")),
#        cbind(datFishLength$DF,which(strGroupID == "DL")))
#)
#datFlatPxLength <- cbind(datFlatPxLength,NA)
#names(datFlatPxLength) <- c("LengthPx","groupID","expID")
#datFlatPxLength <- datFlatPxLength[!is.na(datFlatPxLength$LengthPx), ]
#write.csv(datFlatPxLength,file= paste(strDataExportDir,"/FishLength_Updated.csv",sep=""))
#saveRDS(datFlatPxLength,file= paste(strDataExportDir,"/FishLength_Updated.rds",sep="") ) ## THis is the Processed Register File On 
#strProcDataFileName <- paste0(strDatDir, "/FishLength.RData") <- Original rec - no expid's
#load(file=strProcDataFileName) ##Save With Dataset Idx Identifier


################### LOAD CSV OF PX Larva lengths ################
## *NOTE: 11/6/19 I exported the old data in flat file with GroupID and expID,  **
## As the DF distrution looked wide, and rather bimodal, I measured more larvae to fill in the gaps. 
## The new data is in this csv, leaving the posibility to add more. 
## The scriptlet to run the labelling process on a set of expID is found in auxFunctions.r
datFlatPxLength <- read.csv(file= paste(strDataExportDir,"/FishLength_Updated2.csv",sep=""))
message(paste(" Loading Measured fish length in pixels data ... "))


## Extrach the groupID - GroupName Convention we have been using 
## Recalc First Bout Data based on Validated Info ###
datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
strGroupID <- levels(datTrackedEventsRegister$groupID)
# 
# strGroupName <- names(datFishLength)
# datFlatFrame <- data.frame(
#                 rbind(cbind(datFishLength$LF*DIM_MMPERPX,which(strGroupID == "LL")),
#                       cbind(datFishLength$NF*DIM_MMPERPX,which(strGroupID == "NL")),
#                       cbind(datFishLength$DF*DIM_MMPERPX,which(strGroupID == "DL")))
#                 )
# names(datFlatFrame) <- c("LengthPx","groupID")
# datFlatFrame <- datFlatFrame[!is.na(datFlatFrame$LengthPx), ]
# 


Jagsdata=list(groupID=datFlatPxLength$groupID,
              LarvaLength=datFlatPxLength$LengthPx*DIM_MMPERPX, 
              NTOT=nrow(datFlatPxLength));

varnames1=c("mu","sigma")
burn_in=1000;
steps=100000;
thin=10;
chains=3

library(rjags)
strModelName = "modelFishSize.tmp"
fileConn=file(strModelName)
writeLines(modelNorm,fileConn);
close(fileConn)

sizemodel=jags.model(file=strModelName,data=Jagsdata,n.chains=chains);
update(sizemodel,burn_in)
draw=jags.samples(sizemodel,steps,thin=thin,variable.names=varnames1)

##
dmodelSizeNF <- density(draw$mu[which(strGroupID == "NL"),,],bw=0.01)
dmodelSizeLF <- density(draw$mu[which(strGroupID == "LL"),,],bw=0.01)
dmodelSizeDF <- density(draw$mu[which(strGroupID == "DL"),,],bw=0.01)

dSizeNF <- density(datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "NL"),]$LengthPx*DIM_MMPERPX,bw=0.05)
dSizeLF <- density(datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "LL"),]$LengthPx*DIM_MMPERPX,bw=0.05)
dSizeDF <- density(datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "DL"),]$LengthPx *DIM_MMPERPX,bw=0.05)


#lines(dmodelSizeLF,col=colourHLine[2] , xlim=c(3,6),lty=2 ,lwd=2   )
#lines(dmodelSizeNF,col=colourHLine[1],lty=2,lwd=2)
#lines(dmodelSizeDF,col=colourHLine[3],lty=2,lwd=2)

xquant <- seq(0,6,diff(dSizeNF$x)[1])
XLIM <- c(3.5,5)
YLIM <- c(0,3)
pdistBW <- 2 ## mm/sec
strKern <- "gaussian"
ntail <- 100 #nrow(draw$mu[which(strGroupID == "NL"),,])*0.5
norm <- max(dSizeNF$y)

modelNorm <- max(dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),,],1),
                   sd=tail(draw$sigma[which(strGroupID == "NL"),,],1)))


## Normalized Model ?
sum(diff(dSizeNF$x)[1] * dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),ntail,],1),sd=tail(draw$sigma[which(strGroupID == "NL"),ntail,],1)))
##Normalized empirical density?
sum(diff(dSizeNF$x)[1]*dSizeNF$y)


## FIGURE CONTROL FOR LARVAL SIZE
#### PLOT Fitted Larva Length  ###
pdf(file= paste(strPlotExportPath,"/stat/fig2S1_stat_LarvalLengthsGaussian.pdf" ,sep=""),width = 14,height = 3.5)

layout(matrix(c(1,2,3,4),1,4, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.9,1,1))

line = 8.3 ## SubFig Label Params
lineAxis = 3.0
cex = 1.4
adj  = -3
padj <- -7.2
las <- 1
  
  plot(dSizeNF,col=colourHLine[1] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2 ,cex=cex,cex.axis=cex,cex.lab=cex)
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
         dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),ntail-i,],1),
                                  sd=tail(draw$sigma[which(strGroupID == "NL"),ntail-i,],1))
         ,type='l', col=colourHLine[1],lwd=1,lty=1)
  }
  lines(dSizeNF,col="black",lwd=4,lty=2)
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ) )
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  #mtext("K",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="NF",
         legend=c( paste0("",dSizeNF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[1]),lwd=c(3,1),lty=c(2,1),seg.len = 3,cex=cex ) 
  
  ## 2nd panel
  plot(dSizeLF,col=colourHLine[2] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2,cex=cex,cex.axis=cex,cex.lab=cex )
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
          dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "LL"),ntail-i,],1),
                sd=tail(draw$sigma[which(strGroupID == "LL"),ntail-i,],1)),
          type='l', col=colourHLine[2],lwd=1,lty=1)
  }
  lines(dSizeLF,col="black" , xlim=XLIM,ylim=YLIM ,lwd=4 ,type='l',xlab=NA,ylab=NA,lty=2 )
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  #mtext("L",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="LF",
         legend=c( paste0("",dSizeLF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[2]),lwd=c(3,1),lty=c(2,1),seg.len = 3,cex=cex) 
  
  
  plot(dSizeDF,col=colourHLine[3] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2,cex=cex,cex.axis=cex,cex.lab=cex )
  ##Show Model Fit
  for (i in 1:(ntail-1) )
  {
    lines(xquant,
          dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "DL"),ntail-i,],1),
                sd=tail(draw$sigma[which(strGroupID == "DL"),ntail-i,],1))
          ,type='l',col=colourHLine[3],lwd=1,lty=1)
  }
  lines(dSizeDF,col="black",xlim=XLIM,lwd=4,lty=2)
  
  mtext(side = 1,cex=cex, line = lineAxis, expression("Larval length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  #mtext("M",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj
  legend("topright",title="DF",
         legend=c( paste0("",dSizeDF$n, "# Data "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model  " )),
         col=c("black",colourLegL[3]),lwd=c(3,1),lty=c(2,1),seg.len = 3 ,cex=cex) 
  
  ##MPlot parameter Density
  plot(dmodelSizeNF,col=colourLegL[1] , xlim=XLIM,ylim=c(0,17),lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=1,cex=cex,cex.axis=cex,cex.lab=cex )
  lines(dmodelSizeLF,col=colourLegL[2] , xlim=XLIM ,lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=2 )
  lines(dmodelSizeDF,col=colourLegL[3] , xlim=XLIM ,lwd=4 ,type='l',main=NA,xlab=NA,ylab=NA,lty=3 )
  legend("topright", legend=c("NF","LF","DF"), 
         col=colourLegL,lty=c(1,2,3),lwd=3,seg.len=4,cex=cex)
  mtext(side = 1,cex=cex, line = lineAxis, expression("Estimated mean length (mm) " ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
  #mtext("N",at="topleft",outer=F,side=2,col="black",las=las,font=2,line=line,padj=padj,adj=adj,cex.main=cex,cex=cex) #las=las,line=line,padj=padj,adj=adj

dev.off()

#Statistical Tests
t.test(datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "NL"),]$LengthPx*DIM_MMPERPX,
       datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "DL"),]$LengthPx*DIM_MMPERPX)

t.test(datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "DL"),]$LengthPx*DIM_MMPERPX,
       datFlatPxLength[datFlatPxLength$groupID == which(strGroupID == "LL"),]$LengthPx*DIM_MMPERPX)
strGroupID
t.test(LengthPx*DIM_MMPERPX ~ groupID,data= datFlatPxLength[datFlatPxLength$groupID %in% c(1,3),])
