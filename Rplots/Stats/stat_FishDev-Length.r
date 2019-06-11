## kostasl 11 Jun 2019
## stat Model for comparing larval development across the 3 rearing groups
## we assube 3 gaussians, one for each group, and  infer mean and std dev for each

source("config_lib.R")
source("TrackerDataFilesImport_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")


## The mixture Of Poisson Drawing from Gammaa, gives a negative binomial
modelNorm="model {
         
         for(i in 1:3) {
             mu[i]   ~ dnorm(6,0.0001 )
             prec[i] ~ dunif(0,50) ##dgamma(1,1) ##
             sigma[i] <- sqrt(1/prec[i]) 
        }

         for(i in 1:NTOT){
             LarvaLength[i] ~  dnorm(mu[groupID[i] ],prec[groupID[i]] )

         }
}"


strProcDataFileName <- paste0(strDatDir, "/FishLength.RData")
message(paste(" Loading Measured fish length in pixels data ... "))
load(file=strProcDataFileName) ##Save With Dataset Idx Identifier

## Extrach the groupID - GroupName Convention we have been using 
## Recalc First Bout Data based on Validated Info ###
datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds",sep="") ) ## THis is the Processed Register File On 
strGroupID <- levels(datTrackedEventsRegister$groupID)

strGroupName <- names(datFishLength)
datFlatFrame <- data.frame(
                rbind(cbind(datFishLength$LF,which(strGroupID == "LL")),
                      cbind(datFishLength$NF,which(strGroupID == "NL")),
                      cbind(datFishLength$DF,which(strGroupID == "DL")))
                )

names(datFlatFrame) <- c("LengthPx","groupID")
datFlatFrame <- datFlatFrame[!is.na(datFlatFrame$LengthPx), ]
## Check OUt Length Hist
hist(datFlatFrame$LengthPx*DIM_MMPERPX ,xlim=c(2,6))

Jagsdata=list(groupID=datFlatFrame$groupID,
              LarvaLength=datFlatFrame$LengthPx, 
              NTOT=nrow(datFlatFrame));

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
dmodelSizeNF <- density(draw$mu[which(strGroupID == "NL"),,]*DIM_MMPERPX,bw=0.01)
dmodelSizeLF <- density(draw$mu[which(strGroupID == "LL"),,]*DIM_MMPERPX,bw=0.01)
dmodelSizeDF <- density(draw$mu[which(strGroupID == "DL"),,]*DIM_MMPERPX,bw=0.01)

dSizeNF <- density(datFishLength[!is.na(datFishLength$NF),]$NF*DIM_MMPERPX,bw=0.1)
dSizeLF <- density(datFishLength[!is.na(datFishLength$LF),]$LF*DIM_MMPERPX,bw=0.1)
dSizeDF <- density(datFishLength[!is.na(datFishLength$DF),]$DF*DIM_MMPERPX,bw=0.1)


#lines(dmodelSizeLF,col=colourHLine[2] , xlim=c(3,6),lty=2 ,lwd=2   )
#lines(dmodelSizeNF,col=colourHLine[1],lty=2,lwd=2)
#lines(dmodelSizeDF,col=colourHLine[3],lty=2,lwd=2)

xquant <- seq(0,120,1)
XLIM <- c(3,6)
YLIM <- c(0,1)
pdistBW <- 2 ## mm/sec
strKern <- "gaussian"
ntail <- nrow(draw$mu[which(strGroupID == "NL"),,])*0.2
norm <- max(dSizeNF$y)

modelNorm <- max(dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),,],1),
                   sd=tail(draw$sigma[which(strGroupID == "NL"),,],1)))

##sHOW dATA dENSITY
plot(dSizeLF$x,dSizeLF$y/norm,col=colourHLine[2] , xlim=XLIM,ylim=YLIM ,lwd=2 ,type='l',xlab=NA,ylab=NA )
lines(dSizeNF$x,dSizeNF$y/norm,col=colourHLine[1],lwd=2)
lines(dSizeDF$x,dSizeDF$y/norm,col=colourHLine[3],lwd=2)


#plot(xquant*DIM_MMPERPX,
#      dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),ntail,],1),
#            sd=tail(draw$sigma[which(strGroupID == "NL"),ntail,],1)),type='l',
#      col=colourH[1],lwd=1,xlim=XLIM,ylim=YLIM,xlab=NA,ylab=NA)

##Show Model Fit
for (i in 1:(ntail-1) )
{
  lines(xquant*DIM_MMPERPX,
       dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "NL"),ntail-i,],1),
                                sd=tail(draw$sigma[which(strGroupID == "NL"),ntail-i,],1))/modelNorm
       ,type='l', col=colourH[1],lwd=1)
}

##Show Model Fit
for (i in 1:(ntail-1) )
{
  lines(xquant*DIM_MMPERPX,
        dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "LL"),ntail-i,],1),
              sd=tail(draw$sigma[which(strGroupID == "LL"),ntail-i,],1))/modelNorm,
        type='l', col=colourH[2],lwd=1)
}


##Show Model Fit
for (i in 1:(ntail-1) )
{
  lines(xquant*DIM_MMPERPX,
        dnorm(xquant,mean=tail(draw$mu[which(strGroupID == "DL"),ntail-i,],1),
              sd=tail(draw$sigma[which(strGroupID == "DL"),ntail-i,],1))/modelNorm
        ,type='l',col=colourH[3],lwd=1)
}




#######PLOT RESULTS

