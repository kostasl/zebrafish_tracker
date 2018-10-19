##15-10-2018
### Model how Eye Angle reports distance from prey 
## Used to compare differences in distance estimation between rearing groups 


source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

#+ beta[3]*turn[i]

modelLin <- "model{

  # Likelihood
  for(i in 1:N){
    turn[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- beta[1] + beta[2]*bearing[i] 
  }

  # Prior for beta
  for(j in 1:2){
    beta[j] ~ dnorm(0,0.0001)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  sigma     <- 1/sqrt(inv.var)

}"

modelLin2 <- "model {
	for (i in 1:N){
		turn[i] ~ dnorm(turn.hat[i], tau)
		turn.hat[i] <- beta[1] + beta[2] * bearing[i]
	}
	beta[1] ~ dnorm(0, .0001)
	beta[2] ~ dnorm(0, .0001)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}"

##The Eye Angle Vs Distance Model
## Regression of an exponential Function for Eye Distance
modelExp  <- "model{
  phi_0 ~ dnorm(10,2) # Idle Eye Position
  phi_max ~ dnorm(15,5) # Max Eye Vergence Angle
  lambda ~ dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
  limDist <- max(distMax)
  u1 ~ dunif(0, limDist) ## End Hunt Distance - Close to prey
  u0 ~ dunif(u1, limDist) ##Start Hunt Distance -Far 
  
  # Likelihood
  for(i in 1:N){
    ##Make indicator if hunt event is within sampled Range 
    #if (u1 < distP[i]  & distP[i] < u0) {
    s[i] <- step( distP[i]-u1)*step(u0-distP[i]  ) 

    phi_hat[i] <- phi_0 + s[i] * phi_max* (1-exp(-lambda*(distMax[i] - distP[i] ) )) 
    phi[i] ~ dnorm(phi_hat[i],sigma[s[i]+1] ) ##choose sigma 

  }

  # Prior Sigma On Eye Angle when  In Or Out of hunt region 
  for(j in 1:2){
    #inv.var[j] ~ dgamma(0.01, 0.01)  ##Prior for inverse variance
    sigma[j] ~ dgamma(0.01, 0.01) ##Draw 
  }

}"


##The Eye Angle Vs Distance Model
## Regression of an exponential Function for EyeAngle Vs Distance Fitting Each Hunt Event
## Independently , and obtaining statistics over params of each fit
## H : Number of Hunt Events in Data
## N : vector of number of points in Hunt Event
modelExpInd  <- "model{
  ##Prior
  limDist <- max(distMax)

  # Priors 
  for(i in 1:max(hidx) ) {
    phi_0[i] ~ dnorm(10,2) # Idle Eye Position
    phi_max[i] ~ dnorm(15,5) # Max Eye Vergence Angle
    lambda[i] ~ dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
    u1[i] ~ dunif(0, limDist) ## End Hunt Distance - Close to prey
    u0[i] ~ dunif(u1[i], limDist) ##Start Hunt Distance -Far 

  # Sigma On Eye Angle when  In Or Out of hunt region 
   for(j in 1:2){
    sigma[i,j] ~ dgamma(0.01, 0.01) ##Draw 
    }
  }

  # Likelihood / Separate for each Hunt Event
  for(i in 1:N){
  ##Make indicator if hunt event is within sampled Range 
  #if (u1[hidx[i]] < distP[i]  & distP[i] < u0) 
  s[hidx[i],i] <- step( distP[i]-u1[hidx[i]])*step(u0[ hidx[i] ] - distP[i]  ) 
  
  phi_hat[ hidx[i],i] <- phi_0[hidx[i]] + s[hidx[i],i] * phi_max[hidx[i]]* (1-exp(-lambda[ hidx[i] ]*(distMax[i] - distP[i] ) )) 
  phi[i] ~ dnorm( phi_hat[ hidx[i],i], sigma[hidx[i],s[hidx[i],i]+1] ) ##choose sigma 
  }
}"


####Select Subset Of Data To Analyse

strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 

lEyeMotionDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData.rds",sep="") ) #Processed Registry on which we add )

datEyeVsPreyCombinedAll <-  data.frame( do.call(rbind,lEyeMotionDat ) )

strGroupID <- levels(datTrackedEventsRegister$groupID)


##Add The Empty Test Conditions
#strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
#datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
#datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
#datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]
lRegIdx <- list()
ldatsubSet <-list()
lRegIdx[["LL"]] <- unique(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL"),"RegistarIdx"])
## Get Event Counts Within Range ##
ldatREyePointsLL <- list()
ldatLEyePointsLL <- list()
ldatVEyePointsLL <- list()
lnDatLL          <- list()

##Do all this processing to add a sequence index To The hunt Event + make vergence angle INdex 
for (h in 1:NROW(lRegIdx[["LL"]]) )
{
  ldatsubSet[["LL"]] <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "LL") &
                                      datEyeVsPreyCombinedAll$RegistarIdx %in% lRegIdx[["LL"]][h],]
  ldatLEyePointsLL[[h]] <- cbind(ldatsubSet[["LL"]]$LEyeAngle,
                           as.numeric(ldatsubSet[["LL"]]$DistToPrey),
                           as.numeric(ldatsubSet[["LL"]]$DistToPreyInit ),
                           ldatsubSet[["LL"]]$RegistarIdx,
                           h)
  
  ldatREyePointsLL[[h]] <- cbind(ldatsubSet[["LL"]]$REyeAngle,
                           as.numeric(ldatsubSet[["LL"]]$DistToPrey),
                           as.numeric(ldatsubSet[["LL"]]$DistToPreyInit ),
                           ldatsubSet[["LL"]]$RegistarIdx,
                           h)
  
  ldatVEyePointsLL[[h]] <- cbind(vAngle=ldatsubSet[["LL"]]$LEyeAngle-ldatsubSet[["LL"]]$REyeAngle,
                           distToPrey=as.numeric(ldatsubSet[["LL"]]$DistToPrey),
                           initDistToPrey=as.numeric(ldatsubSet[["LL"]]$DistToPreyInit ),
                           RegistarIdx=ldatsubSet[["LL"]]$RegistarIdx,
                           seqIdx=h)
  
  ldatLEyePointsLL[[h]] <- ldatLEyePointsLL[[h]][!is.na(ldatLEyePointsLL[[h]][,2]),]
  ldatREyePointsLL[[h]] <- ldatREyePointsLL[[h]][!is.na(ldatREyePointsLL[[h]][,2]),]
  ldatVEyePointsLL[[h]] <- ldatVEyePointsLL[[h]][!is.na(ldatVEyePointsLL[[h]][,2]),]
  
  nDatLL[[h]] <- NROW(ldatVEyePointsLL[[h]]) ##Not Used Anymore
}

datVEyePointsLL <- data.frame( do.call(rbind,ldatVEyePointsLL ) ) 


lRegIdx[["NL"]] <- unique(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL"),"RegistarIdx"])
## Get Event Counts Within Range ##
ldatsubSet[["NL"]] <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "NL") &
                                                datEyeVsPreyCombinedAll$RegistarIdx %in% lRegIdx[["NL"]],]
vSeqIdx <- which(lRegIdx[["NL"]][lRegIdx[["NL"]] == ldatsubSet[["NL"]]$RegistarIdx  ]    ) 

datLEyePointsNL <- cbind(ldatsubSet[["NL"]]$LEyeAngle,
                         as.numeric(ldatsubSet[["NL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["NL"]]$DistToPreyInit ),
                          
                         )
datREyePointsNL <- cbind(ldatsubSet[["NL"]]$REyeAngle,
                         as.numeric(ldatsubSet[["NL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["NL"]]$DistToPreyInit ),
                         ldatsubSet[["NL"]]$RegistarIdx)
datVEyePointsNL <- cbind(ldatsubSet[["NL"]]$LEyeAngle-ldatsubSet[["NL"]]$REyeAngle,
                         as.numeric(ldatsubSet[["NL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["NL"]]$DistToPreyInit ),
                         ldatsubSet[["NL"]]$RegistarIdx)


datLEyePointsNL <- datLEyePointsNL[!is.na(datLEyePointsNL[,2]),]
datREyePointsNL <- datREyePointsNL[!is.na(datREyePointsNL[,2]),]
datVEyePointsNL <- datVEyePointsNL[!is.na(datVEyePointsNL[,2]),]



lRegIdx[["DL"]] <- unique(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL"),"RegistarIdx"])
## Get Event Counts Within Range ##
ldatsubSet[["DL"]] <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == "DL") &
                                                datEyeVsPreyCombinedAll$RegistarIdx %in% lRegIdx[["DL"]],]

datLEyePointsDL <- cbind(ldatsubSet[["DL"]]$LEyeAngle,
                         as.numeric(ldatsubSet[["DL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["DL"]]$DistToPreyInit ),
                         ldatsubSet[["DL"]]$RegistarIdx)
datREyePointsDL <- cbind(ldatsubSet[["DL"]]$REyeAngle,
                         as.numeric(ldatsubSet[["DL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["DL"]]$DistToPreyInit ),
                         ldatsubSet[["DL"]]$RegistarIdx)
datVEyePointsDL <- cbind(ldatsubSet[["DL"]]$LEyeAngle-ldatsubSet[["DL"]]$REyeAngle,
                         as.numeric(ldatsubSet[["DL"]]$DistToPrey),
                         as.numeric(ldatsubSet[["DL"]]$DistToPreyInit ),
                         ldatsubSet[["DL"]]$RegistarIdx)

datLEyePointsDL <- datLEyePointsDL[!is.na(datLEyePointsDL[,2]),]
datREyePointsDL <- datREyePointsDL[!is.na(datREyePointsDL[,2]),]
datVEyePointsDL <- datVEyePointsDL[!is.na(datVEyePointsDL[,2]),]


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
pchL <- c(16,2,4)
#
#Thse RC params Work Well to Smooth LF And NF
burn_in=100;
steps=4000;
thin=1;


##Larva Event Counts Slice

nDatNL <- NROW(datVEyePointsNL)
nDatDL <- NROW(datVEyePointsDL)

##Test limit data

#vsamplesLL <- sample (nDatLL,size=nDatLL)
#dataLL=list(phi=datVEyePointsLL[vsamplesLL,1],distP=datVEyePointsLL[vsamplesLL,2],N=NROW(vsamplesLL),distMax=datVEyePointsLL[vsamplesLL,3] );
#vsamplesNL <- sample (nDatNL,size=nDatNL)
#dataNL=list(phi=datVEyePointsNL[vsamplesNL,1],distP=datVEyePointsNL[vsamplesNL,2],N=NROW(vsamplesNL),distMax=datVEyePointsNL[vsamplesNL,3] );
#vsamplesDL <- sample (nDatDL,size=nDatDL)
#dataDL=list(phi=datVEyePointsDL[vsamplesDL,1],distP=datVEyePointsDL[vsamplesDL,2],N=NROW(vsamplesDL),distMax=datVEyePointsDL[vsamplesDL,3] );
## Subset Dat For Speed
datVEyePointsLL_Sub <- datVEyePointsLL[datVEyePointsLL$seqIdx %in% c(1,2),] 
dataLL=list(phi=datVEyePointsLL_Sub$vAngle,
            distP=datVEyePointsLL_Sub$distToPrey ,
            N=NROW(datVEyePointsLL_Sub),
            distMax=datVEyePointsLL_Sub$initDistToPrey,
            hidx=datVEyePointsLL_Sub$seqIdx );

varnames=c("u0","u1","phi_0","phi_max","lambda","sigma","s")



library(rjags)
fileConn=file("model.tmp")
#writeLines(modelGPV1,fileConn);
writeLines(modelExpInd,fileConn);
close(fileConn)

mLL=jags.model(file="model.tmp",data=dataLL);
#update(mLL,burn_in);update(mNL,burn_in);update(mDL,burn_in)
drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)


## compute 2D kernel density, see MASS book, pp. 130-131
nlevels <- 12
z <- kde2d(dataLL$distP, dataLL$phi, n=80)

## Plot the infered function
X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_LL_E.pdf",sep=""))
vX <- seq(0,5,by=0.01)
vY <- median(drawLL$phi_0 ) + median(drawLL$phi_max )*(1-exp(-  median(drawLL$lambda)*( mean(datLEyePointsLL[vsamplesLL,3]) - (vX) ) ) )
vY_u <- median(drawLL$phi_0 ) + median(drawLL$phi_max )*(1-exp(-quantile(drawLL$lambda[1,,1])[4]*( mean(datLEyePointsLL[vsamplesLL,3]) - (vX) ) ) )
vY_l <- median(drawLL$phi_0 ) + median(drawLL$phi_max )*(1-exp(- quantile(drawLL$lambda[1,,1])[2]*( mean(datLEyePointsLL[vsamplesLL,3]) - (vX) ) ) )
plot(dataLL$distP,dataLL$phi,pch=21,xlim=c(0,5),ylim=c(0,80),main="LL", bg=colourP[2],col=colourP[2],cex=0.5)
contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
lines( vX ,vY,xlim=c(0,5),ylim=c(0,55),type="l",col="red",lwd=3)
lines( vX ,vY_u,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
lines( vX ,vY_l,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
#dev.off()

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_LL_E.pdf",sep=""))
X11()
hist(drawLL$lambda[1,,1],main="LL")
#dev.off()
########################
## NL ###
mNL=jags.model(file="model.tmp",data=dataNL);
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)

## Plot the infered function NL

## compute 2D kernel density, see MASS book, pp. 130-131
nlevels <- 12
z <- kde2d(dataNL$distP, dataNL$phi, n=80)

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_NL_E.pdf",sep=""))
X11()
vX <- seq(0,5,by=0.01)
vY <- median(drawNL$phi_0 ) + median(drawNL$phi_max )*(1-exp(-  median(drawNL$lambda)*( mean(datLEyePointsNL[vsamplesNL,3]) - (vX) ) ) )
vY_u <- median(drawNL$phi_0 ) + median(drawNL$phi_max )*(1-exp(-quantile(drawNL$lambda[1,,1])[4]*( mean(datLEyePointsNL[vsamplesNL,3]) - (vX) ) ) )
vY_l <- median(drawNL$phi_0 ) + median(drawNL$phi_max )*(1-exp(- quantile(drawNL$lambda[1,,1])[2]*( mean(datLEyePointsNL[vsamplesNL,3]) - (vX) ) ) )
plot(dataNL$distP,dataNL$phi,pch=20,xlim=c(0,5),ylim=c(0,80),main="NL",col=colourP[3])
contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
lines( vX ,vY,xlim=c(0,5),ylim=c(0,55),type="l",col="red",lwd=3)
lines( vX ,vY_u,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
lines( vX ,vY_l,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_NL_E.pdf",sep=""))
hist(drawNL$lambda[1,,1],main="NL")
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawNL$u1[1,,1],breaks=50,xlim=c(0,7),col="red")
hist(drawNL$u0[1,,1],breaks=50,xlim=c(0,7),add=TRUE,col="red")
#dev.off()

############
### DL ###
mDL=jags.model(file="model.tmp",data=dataDL);
drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)


# Plot the infered function DL

## compute 2D kernel density, see MASS book, pp. 130-131
nlevels <- 12
z <- kde2d(dataDL$distP, dataDL$phi, n=80)

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_DL_E.pdf",sep=""))
X11()
vX <- seq(0,5,by=0.01)
vY <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(-  median(drawDL$lambda)*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
vY_u <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(-quantile(drawDL$lambda[1,,1])[4]*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
vY_l <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(- quantile(drawDL$lambda[1,,1])[2]*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
plot(dataDL$distP,dataDL$phi,pch=20,xlim=c(0,6),ylim=c(0,80),main="DL",col=colourP[1])
contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
lines( vX ,vY,xlim=c(0,5),ylim=c(0,55),type="l",col="red",lwd=3)
lines( vX ,vY_u,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
lines( vX ,vY_l,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
#dev.off()
X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_DL_E.pdf",sep=""))
hist(drawDL$lambda[1,,1],main="DL")
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_DL_E.pdf",sep=""))
hist(drawDL$u1[1,,1],breaks=50,xlim=c(0,7),col="red")
hist(drawDL$u0[1,,1],breaks=50,xlim=c(0,7),add=TRUE,col="red")
#dev.off()




X11()
hist(drawLL$sigma[2,,1],breaks=10000,xlim=c(0,2),col=colourH[1],
     xlab=paste(""),main=paste("During hunt Sigma  ") )

X11()
hist(drawLL$sigma[1,,1],breaks=100,col=colourH[1],
     xlab=paste(" "),main=paste("Outside hunt Sigma ") )


X11()
hist(drawLL$phi_max[1,,1])

X11()
hist(drawLL$phi_0[1,,1])

X11()
hist(drawLL$u0[1,,1],breaks=100,xlim=c(0,5))
hist(drawNL$u0[1,,1],breaks=100,xlim=c(0,5),add=TRUE)
hist(drawDL$u0[1,,1],breaks=100,xlim=c(0,5),add=TRUE)

X11()
hist(drawLL$u1[1,,1],breaks=50,xlim=c(0,5))
hist(drawNL$u1[1,,1],breaks=50,xlim=c(0,5),add=TRUE,col="red")
hist(drawDL$u1[1,,1],breaks=50,xlim=c(0,5),col="blue",add=TRUE)
#########


## Plot the infered function DL

pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_DL.pdf",sep=""))
X11()
vX <- seq(0,5,by=0.01)
vY <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(-  median(drawDL$lambda)*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
vY_u <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(-quantile(drawDL$lambda[1,,1])[4]*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
vY_l <- median(drawDL$phi_0 ) + median(drawDL$phi_max )*(1-exp(- quantile(drawDL$lambda[1,,1])[2]*( mean(datLEyePointsDL[vsamplesDL,3]) - (vX) ) ) )
plot(dataDL$distP,dataDL$phi,pch=20,xlim=c(0,6),ylim=c(0,55),main="DL")
lines( vX ,vY,xlim=c(0,5),ylim=c(0,55),type="l",col="red",lwd=3)
lines( vX ,vY_u,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
lines( vX ,vY_l,xlim=c(0,5),ylim=c(0,55),type="l",col="blue",lwd=2)
dev.off()
X11()
hist(drawDL$lambda[1,,1],main="DL")



ind = 100
##Save the Mean Slope and intercept
##quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2]
muLLa=mean(drawLL$beta[,(steps-ind):steps,1][1,]) 
muLLb=mean(drawLL$beta[,(steps-ind):steps,1][2,])
muNLa=mean(drawNL$beta[,(steps-ind):steps,1][1,])
muNLb=mean(drawNL$beta[,(steps-ind):steps,1][2,])
muDLa=mean(drawDL$beta[,(steps-ind):steps,1][1,])
muDLb=mean(drawDL$beta[,(steps-ind):steps,1][2,])
sig=mean(drawLL$sigma[,(steps-ind):steps,1])
###Plot Density of Slope
dLLb<-density(drawLL$beta[,(steps-ind):steps,1][2,])
dNLb<-density(drawNL$beta[,(steps-ind):steps,1][2,])
dDLb<-density(drawDL$beta[,(steps-ind):steps,1][2,])


pdf(file= paste(strPlotExportPath,"/stat/stat_densityolinregressionslope.pdf",sep=""))
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),main="Density Inference of Turn-To-Prey Slope ")
lines(dLLb,col=colourH[2],xlim=c(0.5,1.2),lwd=3,lty=2)
lines(dNLb,col=colourH[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

### PLot Scatter with regression lines with Conf intervals##
#X11()

pdf(file= paste(strPlotExportPath,"/stat/stat_TurnToPrey_LinearRegression.pdf",sep=""))
plot(lFirstBoutPoints[["DL"]][,1], lFirstBoutPoints[["DL"]][,2],
     main=paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[1] ,pch=pchL[1]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=1,lty=2);
#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(lFirstBoutPoints[["DL"]][,2] ~ lFirstBoutPoints[["DL"]][,1]),col=colourH[4],lwd=1.0) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(lFirstBoutPoints[["LL"]][,1], lFirstBoutPoints[["LL"]][,2],pch=pchL[2],col=colourP[2])
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(lFirstBoutPoints[["LL"]][,2] ~ lFirstBoutPoints[["LL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points(lFirstBoutPoints[["NL"]][,1], lFirstBoutPoints[["NL"]][,2],pch=pchL[3],col=colourP[3])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(lFirstBoutPoints[["NL"]][,2] ~ lFirstBoutPoints[["NL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muNLa,b=muNLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       , pch=pchL,col=colourL)

dev.off()





##Plot Densities Summary
sampLL <- coda.samples(mLL,                      variable.names=varnames,                      n.iter=20000, progress.bar="none")
sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL)
X11()
plot(sampNL)
X11()
plot(sampDL,main="DL")











## N : vector of number of points in Hunt Event
modelExpInd  <- "model{
##Prior

# Prior Sigma On Eye Angle when  In Or Out of hunt region 
for(i in 1:max(hidx)) {
for(j in 1:2){
#inv.var[j] ~ dgamma(0.01, 0.01)  ##Prior for inverse variance
sigma[i,j] ~ dgamma(0.01, 0.01) ##Draw 
}
}

# Likelihood / Separate for each Hunt Event
for(i in 1:N){
phi_0[hidx[i]] ~ dnorm(10,2) # Idle Eye Position
phi_max[hidx[i]] ~ dnorm(15,5) # Max Eye Vergence Angle
lambda[hidx[i]] ~ dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
limDist[hidx[i]] <- max(distMax)
u1[hidx[i]] ~ dunif(0, limDist[hidx[i]]) ## End Hunt Distance - Close to prey
u0[hidx[i]] ~ dunif(u1, limDist[hidx[i]]) ##Start Hunt Distance -Far 


##Make indicator if hunt event is within sampled Range 
#if (u1[hidx[i]] < distP[i]  & distP[i] < u0) 
s[hidx[i],i] <- step( distP[i]-u1[hidx[i]])*step(u0[ hidx[i] ]-distP[i]  ) 

phi_hat[hidx[i],i] <- phi_0[hidx[i]] + s[hidx[i],i] * phi_max[hidx[i]]* (1-exp(-lambda[ hidx[i] ]*(distMax[i] - distP[i] ) )) 
phi[hidx[i],i] ~ dnorm( phi_hat[hidx[i],i], sigma[s[hidx[i],i]+1] ) ##choose sigma 

}"

