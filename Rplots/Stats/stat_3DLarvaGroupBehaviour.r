### Kostas Lagogiannis 2019-06-24 

## 3D Gaussian Model for each group, to discover covariance structure in Undershoot to Distance/Speed 
## ******** No clustering Between slow and Fast Swims ****** 
## I made this to complement the Clustering Method, so as to characterize the overall covariance structure
## Update 17/10/19 : To Model Group Behaviour we need to model individual larvae behaviour and estimate differences between group behaviour
##                  ie do not use events to infer changes in group behaviour


library(rjags)
library(runjags)

source("config_lib.R")
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotCaptureSpeedFit <- function(datSpeed,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,70,1)
  XLIM <- c(0,60)
  YLIM <- c(0,0.15)
  pdistBW <- 2 ## mm/sec
  strKern <- "gaussian"
  #ntail <- NROW(drawMCMC$mu[1,2,,nchain])*0.10
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,2,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,2,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,2,ntail-i,nchain],1)),type='l',col=colourH[colourIdx],lty=2 )
  }
  ##Data
  lines(density(datSpeed$CaptureSpeed,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


## Plots the Data Density and the 2 Gaussians fititng high and low speed capture swims
plotUndeshootClusterFit <- function(datTurn,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(0,2,0.02)
  XLIM <- c(0,2)
  YLIM <- c(0,3)
  pdistBW <- 0.1 ## mm/sec
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  
  plot(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM,ylim=YLIM ,main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,1,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,1,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,1,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  
  lines(density(datTurn$Undershoot,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=XLIM )
  legend("topright",title="",
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 
  
}


plotDistanceClustFit <- function(datDist,drawMCMC,colourIdx,nchain = 1)
{
  xquant <- seq(-0.1,0.8,0.05)
  pdistBW <- DIM_MMPERPX ## Manuall annotation  error is at least 1 px error , so smoothing with this bw is relevant
  strKern <- "gaussian"
  ntail <- min(50,NROW(drawMCMC$mu[1,1,,1])*0.10)
  plot(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8),ylim=c(0,5) ,
       main=NA)
  for (i in 1:(ntail-1) )
  {
    lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[1,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[1,3,ntail-i,nchain],1)),type='l',col=colourHLine[colourIdx],lty=1 )
    #lines(xquant,dnorm(xquant,mean=tail(drawMCMC$mu[2,3,ntail-i,nchain],1),sd=tail(drawMCMC$sigma[2,3,ntail-i,nchain],1)),type='l',col=colourLegL[colourIdx],lty=2 )
  }
  lines(density(datDist$DistanceToPrey,bw=pdistBW,kernel=strKern),col="black",lwd=4,xlim=c(0,0.8) )
  legend("topright",title=NA,
         legend=c( paste("Data Density "), #(Bw:",prettyNum(digits=2, pdistBW ),")" ) ,
                   paste("Model low speed " ),
                   paste("Model high speed " )),
         col=c("black",colourR[4],colourLegL[colourIdx]),lwd=c(3,1,1),lty=c(1,1,2) ) 

}


initfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(#mID=c(rbinom(N,1,0.5)), 
#                                     sigma = matrix(c (  c(runif(1,min=0,max=0.1),runif(1,min=0,max=2)),
#s                                                         c(runif(1,min=0,max=0.1),runif(1,min=0,max=15))  ),nrow=2,byrow=T  ),
#                                     mu  = matrix(c (  c( rnorm(1,mean=1,sd=sqrt(1/10) ), rnorm(1,mean=8,sd=sqrt(1/2) ) ),
#                                                        c( rnorm(1,mean=1, sd=sqrt(1/10) ) , rnorm(1,mean=30, sd=sqrt(1/0.1) )    ) )
#                                                     ,nrow=2,byrow = T  ),
                                     ".RNG.name"="base::Super-Duper",
                                     ".RNG.seed"=round(runif(1,0,60000)) ),
                                     simplify=FALSE)
  return(initlist)
}



##  3D Gaussian Hierarchical  Model of Larvae Hunt Behaviour 
## Estimating Hunt Behaviour per Larvae before inferring mean group behaviour
strmodel3Variables_LarvaHuntBehaviour <- "
var x_rand[NLarv,3];

model {

##Draw capt speed from 2d gaussian
for (i in 1:N)
{
  ##Draw from gaussian model of Each Larva
  c[i,1:3] ~ dmnorm(mu[Lid[i],],prec[Lid[i], , ]) ## data in column 1 and 2
  
}



##Covariance matrix and its inverse -> the precision matrix
## for each Gaussian in the mixture - Single Gaussian  Here -
for  (l in 1:NLarv)
{
  prec[l,1:3,1:3] ~ dwish(R,3)
  cov[l,1:3,1:3]  <- inverse(prec[l,1:3,1:3])  
  
  ## Larva priors Are linked to the Group's Priors
  mu[l,1] ~ dnorm(muG[1,1], tG[1,1]) ##turn ratio
  mu[l,2] ~ dnorm(muG[1,2], tG[1,2]) ##cap speed
  mu[l,3] ~ dnorm(muG[1,3], tG[1,3]) ##Distance prey
  
  ## Synthesize data from the distribution for This Larva
  x_rand[l,] ~ dmnorm(mu[l,],prec[l,,])
  
}

### Make Group Priors 
for  (g in 1:1)
{
  muG[g,1] ~ dnorm(1, 0.01)T(0.0,2) ##turn ratio
  muG[g,2] ~ dnorm(25,0.01)T(0,) ##cap speed
  muG[g,3] ~ dnorm(0.1,0.1)T(0,) ##Distance prey
  
  tG[g,1] ~ dgamma(10,3)
  tG[g,2] ~ dgamma(4,25) ##sd < 20
  tG[g,3] ~ dgamma(10,3)
}

for(i in 1:3){
  for(j in 1:3){
    R[i,j] <- equals(i,j)*1e-4
  }

}

} "



strModelPDFFileName <- "/stat/UndershootAnalysis/fig7S1-stat_modelCaptureSpeedVsUndershootAndDistance_Valid.pdf"
strDataPDFFileName <- "/stat/UndershootAnalysis/fig7-UndershootCaptureSpeedCV_scatter_Valid.pdf"
strCaptSpeedDensityPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelCaptureSpeed_Valid.pdf"
strUndershootDensityPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelUndershoot_Valid.pdf"
strDistanceDensityPDFFileName <- "/stat/UndershootAnalysis/stat_modelDistance_Valid.pdf"
strModelCovarPDFFileName <- "/stat/UndershootAnalysis/fig7-stat_modelCaptureSpeedVsUndershootAndDistance_COVar.pdf"

datTrackedEventsRegister <- readRDS( paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_ToValidate.rds","",sep="") ) ## THis is the Processed Register File On 
#lMotionBoutDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_MotionBoutData_SetC.rds",sep="") ) #Processed Registry on which we add )
#lEyeMotionDat <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData_SetC",".rds",sep="")) #
lFirstBoutPoints <-readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_Validated.rds",sep="")) 

datTurnVsStrikeSpeed_NL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$NL[,"Turn"]/lFirstBoutPoints$NL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$NL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$NL[,"DistanceToPrey"],Validated= lFirstBoutPoints$NL[,"Validated"] )
datTurnVsStrikeSpeed_LL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$LL[,"Turn"]/lFirstBoutPoints$LL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$LL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$LL[,"DistanceToPrey"],Validated= lFirstBoutPoints$LL[,"Validated"] )
datTurnVsStrikeSpeed_DL <- data.frame( cbind(TurnRatio=lFirstBoutPoints$DL[,"Turn"]/lFirstBoutPoints$DL[,"OnSetAngleToPrey"],CaptureSpeed=lFirstBoutPoints$DL[,"CaptureSpeed"]),DistanceToPrey=lFirstBoutPoints$DL[,"DistanceToPrey"],Validated= lFirstBoutPoints$DL[,"Validated"] )

###Filter For Hunting Where Prey Is approached into strike distance, rather than Initial Prey Distance being within strike Distance
#datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$OnSetDistance > 0.6,]
#datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$OnSetDistance > 0.6,]
###Validated Only
replace(datTurnVsStrikeSpeed_NL$Validated, is.na(datTurnVsStrikeSpeed_NL$Validated), 0)
replace(datTurnVsStrikeSpeed_LL$Validated, is.na(datTurnVsStrikeSpeed_LL$Validated), 0)
replace(datTurnVsStrikeSpeed_DL$Validated, is.na(datTurnVsStrikeSpeed_DL$Validated), 0) 

###Validated Only
datTurnVsStrikeSpeed_NL <- datTurnVsStrikeSpeed_NL[datTurnVsStrikeSpeed_NL$Validated == 1, ]
datTurnVsStrikeSpeed_LL <- datTurnVsStrikeSpeed_LL[datTurnVsStrikeSpeed_LL$Validated == 1, ]
datTurnVsStrikeSpeed_DL <- datTurnVsStrikeSpeed_DL[datTurnVsStrikeSpeed_DL$Validated == 1, ]

datTurnVsStrikeSpeed_ALL <- rbind(datTurnVsStrikeSpeed_NL,datTurnVsStrikeSpeed_LL,datTurnVsStrikeSpeed_DL)


#### LOAD Capture First-Last Bout hunting that include the cluster classification - (made in stat_CaptureSpeedVsDistanceToPrey)
datCapture_NL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_NL_clustered.rds",sep="")) 
datCapture_LL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_LL_clustered.rds",sep="")) 
datCapture_DL <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData_wCapFrame_DL_clustered.rds",sep="")) 

#######################################
########## PCA  - FACTOR ANALYSIS ####
#######################################
# Check Correlation Of UNdershoot With Hunt POwer
##Take all expID from the successful hunt Events we have extracted hunt variables from 
vexpID <- list(LF = datTrackedEventsRegister[datCapture_LL$RegistarIdx,]$expID,
               NF=datTrackedEventsRegister[datCapture_NL$RegistarIdx,]$expID,
               DF=datTrackedEventsRegister[datCapture_DL$RegistarIdx,]$expID)

## Merge EXP ID
## Add Exp ID Column - Signifying Which Larvae Executed the Capture Success Hunt- 
datCapture_LF_wExpID <- cbind(datCapture_LL,expID=vexpID$LF,groupID=2)
datCapture_NF_wExpID <- cbind(datCapture_NL,expID=vexpID$NF,groupID=3)
datCapture_DF_wExpID <- cbind(datCapture_DL,expID=vexpID$DF,groupID=1)
datCapture_ALL_wExpID <- rbind(datCapture_LF_wExpID,datCapture_NF_wExpID,datCapture_DF_wExpID)
###Empirical Distribution
datHuntLarvaStat <- aggregate(datCapture_ALL_wExpID,by=list(datCapture_ALL_wExpID$expID),mean)

##
##
steps <-10000
nchains <- 3
nthin <- 5
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","cov","x_rand","muG","tG","NLarv") #Mixture Model
##Make Serial Larvae ID, that link each hunt event to an individual larva 
ldata_LF <- with(datCapture_LF_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=NROW(unique(expID)) ) }) ##Live fed
ldata_NF <- with(datCapture_NF_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=NROW(unique(expID))  ) }) ##Live fed
ldata_DF <- with(datCapture_DF_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Lid=as.numeric(as.factor(as.numeric(expID)) ) ,N=NROW(expID), NLarv=NROW(unique(expID))  ) }) ##Live fed
ldata_ALL <-with(datCapture_ALL_wExpID, {list(c=cbind(Undershoot,CaptureSpeed,DistanceToPrey),Lid=as.numeric(as.factor(as.numeric(expID)) ),Gid=groupID ,N=NROW(expID),NLarv=NROW(unique(expID))  ) }) ##Live fed list(c=datTurnVsStrikeSpeed_ALL,N=NROW(datTurnVsStrikeSpeed_ALL)) ##Dry fed


jags_model_LF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_LF, 
                         n.adapt = 100, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_LF$N))
update(jags_model_LF, 300)
draw_LF=jags.samples(jags_model_LF,steps,thin=nthin,variable.names=str_vars)

## Not Fed
jags_model_NF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_NF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_NF$N)) 
update(jags_model_NF,300)
draw_NF=jags.samples(jags_model_NF,steps,thin=nthin,variable.names=str_vars)

## Dry  Fed
jags_model_DF <- jags.model(textConnection(strmodel3Variables_LarvaHuntBehaviour), data = ldata_DF, 
                         n.adapt = 500, n.chains = nchains, quiet = F,inits=initfunct(nchains,ldata_DF$N))
update(jags_model_DF, 300)
draw_DF=jags.samples(jags_model_DF,steps,thin=nthin,variable.names=str_vars)


message("Mean LF Und:", prettyNum( mean(draw_LF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_LF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_LF$muG[,3,,1]), digits=3)
)


message("Mean NF Und:", prettyNum( mean(draw_NF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_NF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_NF$muG[,3,,1]), digits=3)
)

message("Mean DF Und:", prettyNum( mean(draw_DF$muG[,1,,]) , digits=3),
        " Speed : ",prettyNum( mean(draw_DF$muG[,2,,1]), digits=3),
        " Distance : ",prettyNum(mean(draw_DF$muG[,3,,1]), digits=3)
)

save(draw_NF,draw_LF,draw_DF,file = paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouModel_RJags.RData"))
## ALL  groups
#jags_model_ALL <- jags.model(textConnection(strmodel_capspeedVsUndershoot_Mixture), data = ldata_ALL, 
                            #n.adapt = 500, n.chains = 3, quiet = F)
#update(jags_model_ALL, 300)
#draw_ALL=jags.samples(jags_model_ALL,steps,thin=2,variable.names=str_vars)


#load(paste0(strDataExportDir,"stat_Larval3DGaussianBehaviouModel_RJags.RData"))

schain <- 5:10
stail <- 300
plot(density(tail(draw_LF$muG[,2,,1], 150) ),ylim=c(0,1),xlim=c(0,80))
lines(density(tail(draw_NF$muG[,2,,schain], 150)))
lines(density(tail(draw_DF$muG[,2,,schain], 150)))

plot(density(tail(draw_LF$muG[,1,,3], 150) ),ylim=c(0,1),xlim=c(0,2))
lines(density(tail(draw_LF$muG[,1,,1], 150) ),ylim=c(0,1),xlim=c(0,2))
lines(density(tail(draw_LF$muG[,1,,2], 150) ),ylim=c(0,1),xlim=c(0,2))

##Undershoot Posterior Group Vs INdividual Density
with(draw_LF,{
  plot(density(tail(muG[,1,,], 100) ),ylim=c(0,16),xlim=c(0,2),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,1,,],stail)),lty=2)
})

##Undershoot Posterior Group Vs INdividual Density
with(draw_DF,{
  plot(density(tail(muG[,1,,], 100) ),ylim=c(0,16),xlim=c(0,2),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,1,,],stail)),lty=2)
})




##Speed Posterior Group Vs INdividual Density
with(draw_LF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed LF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==2,]$CaptureSpeed,bw=2),col="blue",lwd=2)

##Speed Posterior Group Vs INdividual Density
with(draw_DF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed DF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
   lines(dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail)) )) ,col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==1,]$CaptureSpeed,bw=2),col="blue",lwd=2)

##Speed Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed NF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==3,]$CaptureSpeed,bw=2),col="blue",lwd=2)




##Distance Posterior Group Vs INdividual Density
with(draw_LF,{
  plot(density(tail(muG[,3,,1], 100) ),ylim=c(0,16),xlim=c(0,1),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,3,,1],stail)),lty=2)
})



##Covariance 
### Lid,Matrix i,Matrix j,Sample,chain
draw_LF$cov[,2,2,1,1]

##CONVERGENCE - CHECK
#Idx: Lid,Variable (1Under),Sample Row,Chain
plot(density(tail(draw_LF$muG[,2,,1],1000) ),type='l')
lines(density(tail(draw_LF$muG[,2,,2],1000) ),col="red")
lines(density(tail(draw_LF$muG[,2,,3],1000) ),col="blue")


plot(density(draw_LF$mu[9,1,,schain] ) )

### Estimated For Each Larva - Plot Group Population
##Plot Distance Density

plot(density(sapply(tail(draw_LF$mu[,3,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Distance to Prey",ylim=c(0,6)) ##Mean Group Undershoot From Mean Of Each Larva
ldist <- list()
for ( i in (1:28) )
  ldist[[i]]  <-sapply(tail(draw_LF$mu[i,3,,],stail),mean)
##Overlay The Density From The Estimated Mean Overshoot Of Each Larva
lines(density( unlist(lapply(ldist,mean) ) ) )

lines(density(sapply(tail(draw_NF$mu[,3,,],stail),mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
ldist <- list()
for ( i in (1:28) )
  ldist[[i]]  <-sapply(tail(draw_NF$mu[i,3,,],stail),mean)
##Overlay The Density From The Estimated Mean Overshoot Of Each Larva
lines(density( unlist(lapply(ldist,mean) ) ) )


lines(density(sapply(tail(draw_DF$mu[1,3,,],stail) ,mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
ldist <- list()
for ( i in (1:28) )
  ldist[[i]]  <-sapply(tail(draw_DF$mu[i,3,,],stail),mean)
##Overlay The Density From The Estimated Mean Overshoot Of Each Larva
lines(density( unlist(lapply(ldist,mean) ) ) )

##Plot Speed Density
plot(density(sapply(tail(draw_LF$mu[,2,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Capture Speed",ylim=c(0,0.1)) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,2,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,2,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

##Plot Undershoot Density / Mean Sample point Across larva 
plot(density(sapply(tail(draw_LF$mu[,1,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Turn ratio") ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,1,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,1,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

#Idx: Lid,Variable (2=SpeedDist Covar),Sample Row,Chain
## 
plot(sapply((draw_LF$rho[,2,,1]),mean),type='l')


### Estimate  densities  ###

load(paste0(strDataExportDir,"stat_CaptSpeedVsUndershootAndDistance_RJags.RData"))

nContours <- 6
ntail <- 1200 #NROW(draw_NF$mu[1,1,,1])*0.20



zLL <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,2,,],ntail)),n=180)
zNL <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,2,,],ntail)),n=180)
zDL <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,2,,],ntail)),n=180)
#zALL <- kde2d(c(tail(draw_ALL$mu[,1,,1],ntail)), c(tail(draw_ALL$mu[,2,,1],ntail)),n=80)

zLLD <- kde2d(c(tail(draw_LF$mu[,1,,],ntail)), c(tail(draw_LF$mu[,3,,],ntail)),n=180)
zNLD <- kde2d(c(tail(draw_NF$mu[,1,,],ntail)), c(tail(draw_NF$mu[,3,,],ntail)),n=180)
zDLD <- kde2d(c(tail(draw_DF$mu[,1,,],ntail)), c(tail(draw_DF$mu[,3,,],ntail)),n=180)


zLLS <- kde2d(c(tail(draw_LF$mu[,3,,],ntail)), c(tail(draw_LF$mu[,2,,],ntail)),n=180)
zNLS <- kde2d(c(tail(draw_NF$mu[,3,,],ntail)), c(tail(draw_NF$mu[,2,,],ntail)),n=180)
zDLS <- kde2d(c(tail(draw_DF$mu[,3,,],ntail)), c(tail(draw_DF$mu[,2,,],ntail)),n=180)

## Check out the covar coeffient , compare estimated densities
pBw   <- 0.02
## Strike Cluster Only (Fast speed) draw_NF$mu[,2,,]
## The Undershoot To Capt. Speed covar coefficient
dLLb_rhoUS<-density(tail(draw_LF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)  ## Undershoot-Speed Covar
dNLb_rhoUS<-density(tail(draw_NF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUS<-density(tail(draw_DF$rho[,1,,],ntail),kernel="gaussian",bw=pBw)
#The Speed - Distance covar coefficient
dLLb_rhoSD<-density(tail(draw_LF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoSD<-density(tail(draw_NF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoSD<-density(tail(draw_DF$rho[,2,,],ntail),kernel="gaussian",bw=pBw)
##The UNdershoot Distance covar 
dLLb_rhoUD<-density(tail(draw_LF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)
dNLb_rhoUD<-density(tail(draw_NF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)
dDLb_rhoUD<-density(tail(draw_DF$rho[,3,,],ntail),kernel="gaussian",bw=pBw)

#dALLb_rho<-density(tail(draw_ALL$rho[,,1],ntail),kernel="gaussian",bw=pBw)

save(dLLb_rhoSD,dNLb_rhoSD,dDLb_rhoSD,file = paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))
## ALL  


load(paste0(strDataExportDir,"stat_CaptSpeedVsDistance_Covariance_RJags.RData"))
###Check COnv
draw <- draw_NF
plot(draw$mu[1,1,,1],type='l',ylim=c(0,2),col=rfc(nchains)[1] )
lines(draw$mu[1,1,,2],type='l',ylim=c(0,2),col=rfc(nchains)[2] )
lines(draw$mu[1,1,,3],type='l',ylim=c(0,2),col=rfc(nchains)[3] )
lines(draw$mu[1,1,,4],type='l',ylim=c(0,2),col=rfc(nchains)[4] )
lines(draw$mu[1,1,,5],type='l',ylim=c(0,2),col=rfc(nchains)[5] )
##Get the synthesized data:

#plot(tail(draw_NF$x_rand[1,,1],ntail ),tail(draw_NF$x_rand[2,,1],ntail ),col=colourH[1])
#points(tail(draw_LF$x_rand[1,,1],ntail ),tail(draw_LF$x_rand[2,,1],ntail ),col=colourH[2])
#points(tail(draw_DF$x_rand[1,,1],ntail ),tail(draw_DF$x_rand[2,,1],ntail ),col=colourH[3])




### MAIN COVARIANCE PLOT  (Fast Cluster)##
###Show covariance In the High Speed Capture Cluster ##

pdf(file= paste0(strPlotExportPath,"/stat/fig7S1-stat_3Dmodel_TurnRatioVsSpeedAndDistance.pdf"),width=14,height=7,
    title=" 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
nContours <- 5
### Show Speed Fit ###
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 2.7
lineXAxis = 3.0
lineTitle = 0.5

cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1
par(mar = c(3.9,4.7,3.5,1))

#layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5),2,6, byrow = TRUE))
layout(matrix(c(1,1,2,2),1,4, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )



  ## Plot the mean of the 2D Models ##
  ##Collect Draws from all chains
  plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1], xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
  points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
  points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
  #points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])
  
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Capture speed (mm/sec)  " ))
  #mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
  
  contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
  contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
  contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
  
  legend("topright",
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                    bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                    bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                    #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
         ),
         pch=pchL, col=colourLegL,cex=cex)
  ###############
  
  
  ## Distance To Prey Vs Turn Ratio##
  plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0.5,1.5),ylim=c(0,0.5),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
  points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),col=colourHPoint[2],pch=pchL[2])
  points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),col=colourHPoint[3],pch=pchL[3])
  #points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])
  
  mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
  mtext(side = 2,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))
  
  contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
  contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
  contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
  contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)
  
  #mtext("B",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
dev.off()



      ### COVARIANCES 
      ## Plot the covariance Coefficients##
  
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_SpeedVsTurn_Covar.pdf"),width=7,height=7,
      title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
     ##Speed TO Distance Covariance Coeff      
  
      plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,5),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)
      #lines(dALLb_rhoUS,col=colourLegL[4],lwd=3,lty=4)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Capture speed and turn ratio "  ))
      #mtext("C",at="topleft",outer=outer,side=2,col="black",font=2  ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
      
      ###############
dev.off()      
      
      
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_SpeedVsDistance_Covar.pdf"),width=7,height=7,
      title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
      ##Speed TO Distance Covariance Coeff    
      plot(dNLb_rhoSD,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,5),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoSD,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoSD,col=colourLegL[3],lwd=3,lty=3)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Capture speed and distance to prey"  ))
      
      
      
      
      legend("topleft",
             legend=c(  expression (),
                        bquote(NF[""] ),
                        bquote(LF[""] ),
                        bquote(DF[""]   )
                        
             ),
             lty=c(1,2,3),lwd=3,fill=colourLegL,cex=cex)
      #mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
dev.off()     
      
      
pdf(file= paste0(strPlotExportPath,"/stat/stat_3dmodel_TurnVsDistance_Covar.pdf"),width=7,height=7,
          title="Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
      ##Speed TO Distance Covariance Coeff
      plot(dNLb_rhoUD,col=colourLegL[1],xlim=c(-0.5,0.5),lwd=3,lty=1,ylim=c(0,5),
           main=NA, #"Density Inference of Turn-To-Prey Slope ",
           xlab=NA,ylab=NA,cex=cex,cex.axis=cex) #expression(paste("slope ",gamma) ) )
      lines(dLLb_rhoUD,col=colourLegL[2],lwd=3,lty=2)
      lines(dDLb_rhoUD,col=colourLegL[3],lwd=3,lty=3)
      mtext(side = 1,cex=cex, line = lineXAxis, expression("Covariance coefficient"  ))
      mtext(side = 2,cex=cex, line = lineAxis, expression("Density function " ))
      mtext(side = 3,cex=cex, line = lineTitle, expression("Turn ratio and distance to prey"  ))
      
#      mtext("E",at="topleft",outer=F,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)
dev.off()


### 3D density figure of Means ##
#library(MASS)
#library(plotly)
#den3d <- kde2d(x, y)
#persp(zLLD, box=FALSE)
## the new part:
#plot_ly(x=zLLD$x, y=zLLD$y, z=zLLD$z) %>% add_surface()

library( rgl )
#library(plot3D)
# Static chart
ntail <- 500
# datMu3D <-  data.frame( cbind.data.frame( TurnR=as.numeric(tail(draw_NF$mu[,1,,1],ntail)),CSpeed=tail(draw_NF$mu[,2,,1],ntail),Dist=tail(draw_NF$mu[,3,,1],ntail),col=colourHL[1])  )
# datMu3D <- rbind(datMu3D,
#                  data.frame( cbind.data.frame( TurnR=tail(draw_LF$mu[,1,,1],ntail),CSpeed=tail(draw_LF$mu[,2,,1],ntail),Dist=tail(draw_LF$mu[,3,,1],ntail),col=colourHL[2])  ))
# datMu3D <- rbind(datMu3D,
#                  data.frame( cbind.data.frame( TurnR=tail(draw_DF$mu[,1,,1],ntail),CSpeed=tail(draw_DF$mu[,2,,1],ntail),Dist=tail(draw_DF$mu[,3,,1],ntail),col=colourHL[3])  ))


datMu3D <-  data.frame( cbind.data.frame(
                        TurnR=tail(draw_LF$muG[,1,,1], ntail) ,
                        CSpeed=tail(draw_LF$muG[,2,,1], ntail),
                        Dist=tail(draw_LF$muG[,3,,1], ntail), col=colourHL[2])  )

datMu3D <- rbind(datMu3D,
                  data.frame( cbind.data.frame(
                   TurnR=tail(draw_NF$muG[,1,,1], ntail) ,
                   CSpeed=tail(draw_NF$muG[,2,,1], ntail),
                   Dist=tail(draw_NF$muG[,3,,1], ntail), col=colourHL[3])  )
              )

datMu3D <- rbind(datMu3D,
                 data.frame( cbind.data.frame(
                   TurnR=tail(draw_DF$muG[,1,,1], ntail) ,
                   CSpeed=tail(draw_DF$muG[,2,,1], ntail),
                   Dist=tail(draw_DF$muG[,3,,1], ntail), col=colourHL[1])  )
                 
                 )


open3d()
rgl::plot3d( x=datMu3D$TurnR, y=datMu3D$CSpeed, z=datMu3D$Dist, col = datMu3D$col, type = "s", radius = 1.3,
             #xlab="Turn Ratio", ylab="Capture Speed (mm/sec)",zlab="Distance to prey (mm)",
             xlab="Turn Ratio", ylab="Speed",zlab="Distance",
             #xlim=c(0.2,1.8), ylim=c(5,60), zlim=c(0,0.8),
             box = FALSE ,aspect = TRUE
             #,expand = 1.5
             )

rgl::rgl.viewpoint(60,10)
rgl::rgl.viewpoint(45,45)

#decorate3d(
#           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
#           top = TRUE, aspect = FALSE, expand = 1.03,cex=cex)

###Warning External Editing of PDF may fail. a Ghostscript conversion can fix this: 
#Use : gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/screen -dNOPAUSE -dBATCH  -dQUIET -sOutputFile=output.pdf fig7_Modelballs3D_SpeedVsTurn_view.pdf
rgl::rgl.postscript( paste0(strPlotExportPath,"/fig7_Modelballs3D_Perspective_TurnVsSpeed_view.pdf"),fmt="pdf",drawText = FALSE )
rgl::rgl.postscript( paste0(strPlotExportPath,"/fig7_Modelballs3D_Perspective_TurnVsSpeed_view.tex"),fmt="tex",drawText = TRUE )
rgl::rgl.snapshot( paste0(strPlotExportPath,"/fig7_Modelballs3D_DistVsTurn_view"),fmt="png" )

## I use the tex Axis doc to combine 3d Fig with axis text. I then compile a pdf, which I import into inkScape (using Cairo to rasterize image), so I can adjust size, add axis labels etc etc.
## END OF 3D plot Messing with exporting
## 




#####
####################################  Summary FIGURE Of Each Pair + A Covariance Of Distance To Turn Ratio  #####
## PLot Model / Means and covariance ##
## Open Output PDF 
pdf(file= paste(strPlotExportPath,strModelPDFFileName,sep=""),width=14,height=7,
    title="A 3D statistical model for Capture Strike speed / Turn Ratio / Distance to Prey")

### Show Speed Fit ###
outer = FALSE
line = 1 ## SubFig Label Params
lineAxis = 2.7
lineXAxis = 3.0
cex = 1.4
adj  = 3.5
padj <- -8.0
las <- 1


layout(matrix(c(1,2,3,4),2,2, byrow = TRUE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.7,2,1))

## Plot the mean of the 2D Models ##
##Collect Draws from all chains
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1], xlim=c(0.5,1.5),ylim=c(10,50),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Capture Speed (mm/sec)  " ))
mtext("A",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)

contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNL, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)

legend("topright",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
                  ),
       pch=pchL, col=colourLegL,cex=cex)
###############

## Distance To Prey Vs Speed ##
plot(tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0,0.5),ylim=c(10,50),ylab=NA,xlab=NA ,cex=cex,cex.axis=cex )
points(tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))
mtext(side = 2,cex=cex, line = lineAxis, expression(" Capture Speed (mm/sec)" ))

contour(zDLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNLS, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


#legend("topleft",
#       legend=c(  expression (),
#                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
#                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
#                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
#                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
#      ),
#       pch=pchL, col=colourLegL)
mtext("B",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)


## Distance To Prey Vs Turn Ratio##
plot(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),col=colourHPoint[1],pch=pchL[1],  xlim=c(0.5,1.5),ylim=c(0,0.5),ylab=NA,xlab=NA,cex=cex,cex.axis=cex  )
points(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),col=colourHPoint[2],pch=pchL[2])
points(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),col=colourHPoint[3],pch=pchL[3])
#points(tail(draw_ALL$mu[2,1,,1],ntail),tail(draw_DF$mu[2,2,,1],ntail),col=colourH[4],pch=pchL[4])

mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
mtext(side = 2,cex=cex, line = lineAxis, expression("Distance to prey (mm)  " ))

contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col="black",lwd=1)
contour(zDLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL [3],lty=2)
contour(zLLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[2],lty=2)
contour(zNLD, drawlabels=FALSE, nlevels=nContours,add=TRUE,col=colourLegL[1],lty=2)#contour(zALL, drawlabels=FALSE, nlevels=nContours,add=TRUE)


#legend("topright",
#       legend=c(  expression (),
#                  bquote(NF["e"] ~ '#' ~ .(ldata_NF$N)  ),
#                  bquote(LF["e"] ~ '#' ~ .(ldata_LF$N)  ),
#                  bquote(DF["e"] ~ '#' ~ .(ldata_DF$N)  )
#                  #, bquote(All ~ '#' ~ .(ldata_ALL$N)  )
#       ),
#       pch=pchL, col=colourLegL)
mtext("C",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)



### Show Covar Of Undershoot to Distance  Membership
plot(dNLb_rhoUD,col=colourLegL[1],xlim=c(-1,1),ylim=c(0.4,10),lwd=4,lty=1,main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dLLb_rhoUD,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUD,col=colourLegL[3],lwd=3,lty=3)
#lines(density(draw_ALL$pS),col=colourLegL[4],lwd=3,lty=4)

legend("topleft",
       legend=c(  expression (),
                  bquote(NF[""] ~ '#' ~ .(ldata_NF$N)  ),
                  bquote(LF[""] ~ '#' ~ .(ldata_LF$N)  ),
                  bquote(DF[""] ~ '#' ~ .(ldata_DF$N)  )
                  #,bquote(ALL ~ '#' ~ .(ldata_ALL$N)  )
       ), ##paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)

mtext(side = 2,cex=cex, line = lineAxis, expression("Density ") )
mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Turn ratio to distance covariance " ) )  )
#mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=cex)
mtext("D",at="topleft",outer=outer,side=2,col="black",font=2      ,las=1,line=line,padj=padj,adj=3,cex.main=cex,cex=cex)


dev.off()

##################################################
####################################################

### Show Covar Of Undershoot to Speed  Membership
plot(dNLb_rhoUS,col=colourLegL[1],xlim=c(-1,1),ylim=c(0.4,10),lwd=4,lty=1,main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex )
lines(dLLb_rhoUS,col=colourLegL[2],lwd=3,lty=2)
lines(dDLb_rhoUS,col=colourLegL[3],lwd=3,lty=3)


##### Individual Rand Vars Fit ###
## Capture Speeds ##
pdf(file= paste(strPlotExportPath,strCaptSpeedDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
npchain<-3
plotCaptureSpeedFit(datTurnVsStrikeSpeed_NL,draw_NF,1,npchain)
title(main="Model capture Speed")
plotCaptureSpeedFit(datTurnVsStrikeSpeed_LL,draw_LF,2,npchain)
plotCaptureSpeedFit(datTurnVsStrikeSpeed_DL,draw_DF,3,npchain)

dev.off()

## Undershoot  ##
pdf(file= paste(strPlotExportPath,strUndershootDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotUndeshootClusterFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model undershoot on 1st turn to prey")
plotUndeshootClusterFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotUndeshootClusterFit(datTurnVsStrikeSpeed_DL,draw_DF,3)
dev.off()

## Distance ##
pdf(file= paste(strPlotExportPath,strDistanceDensityPDFFileName ,sep=""))
par(mar = c(3.9,4.3,1,1))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
plotDistanceClustFit(datTurnVsStrikeSpeed_NL,draw_NF,1)
title(main="Model distance from prey prior to capture")
plotDistanceClustFit(datTurnVsStrikeSpeed_LL,draw_LF,2)
plotDistanceClustFit(datTurnVsStrikeSpeed_DL,draw_DF,3)

dev.off()
## plot 
##plot(xquant,dnorm(xquant,mean=tail(draw_NF$mu[2,2,,1],1),sd=tail(draw_NF$sigma[2,2,,1],1)),type='l',col=colourH[1],lty=1 )



## SLow Clust
clust <- 2
###Undershoot-Speed Covar
plot(density(draw_LF$sigma[clust,1,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,1,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,1,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,1,,1]),
     col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,1,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,1,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))

###Speed Distance
plot(density(draw_LF$sigma[clust,3,,1]*draw_LF$sigma[clust,2,,1]*draw_LF$rho[clust,2,,1]),
     col=colourLegL[2],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_NF$sigma[clust,3,,1]*draw_NF$sigma[clust,2,,1]*draw_NF$rho[clust,2,,1]),
      col=colourLegL[1],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))
lines(density(draw_DF$sigma[clust,3,,1]*draw_DF$sigma[clust,2,,1]*draw_DF$rho[clust,2,,1]),
      col=colourLegL[3],xlim=c(-1.0,1),lwd=3,lty=1,ylim=c(0,4))


#mcmc_samples <- coda.samples(jags_model, c("mu", "rho", "sigma", "x_rand"),                             n.iter = 5000)

  ######################################################################
 ###         PLOT EMPIRICAL                              ##############
#####################################################################
###        UNdershoot Vs Capture speed               ###
densNL <-  kde2d(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,n=80)
densLL <-  kde2d(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,n=80)
densDL <-  kde2d(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,n=80)

covLL <- cov( 1/datTurnVsStrikeSpeed_LL$Undershoot,datTurnVsStrikeSpeed_LL$CaptureSpeed)
covDL <- cov( 1/datTurnVsStrikeSpeed_DL$Undershoot,datTurnVsStrikeSpeed_DL$CaptureSpeed)
covNL  <- cov( 1/datTurnVsStrikeSpeed_NL$Undershoot,datTurnVsStrikeSpeed_NL$CaptureSpeed)


pdf(file= paste(strPlotExportPath,"distal",strDataPDFFileName,sep=""))
layout(matrix(c(1,2,3),3,1, byrow = FALSE))
##Margin: (Bottom,Left,Top,Right )
par(mar = c(3.9,4.5,1,1))

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,60),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$CaptureSpeed ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
contour(densNL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex)  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$CaptureSpeed,col=colourLegL[2],
     ylim=c(0,60),xlim=c(0,2),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$CaptureSpeed ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
contour(densLL, drawlabels=FALSE, nlevels=7,add=TRUE,col=colourL[4],lty=2,lwd=1)
mtext(side = 2,cex=cex, line = lineAxis-0.7, expression("Capture Speed (mm/sec) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$CaptureSpeed,col=colourLegL[3],ylim=c(0,60),xlim=c(0,2),
     xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$CaptureSpeed ~ datTurnVsStrikeSpeed_DL$Undershoot)
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

plot(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$DistanceToPrey,col=colourLegL[1],
     xlab=NA,ylab=NA,ylim=c(0,1.0),xlim=c(0,2),main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_NL$DistanceToPrey ~ datTurnVsStrikeSpeed_NL$Undershoot)
abline(lFit,col=colourLegL[1],lwd=3.0) ##Fit Line / Regression
legend("topright",
       legend=paste("NF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex )  #prettyNum(digits=3, cov(datTurnVsStrikeSpeed_NL$Undershoot, datTurnVsStrikeSpeed_NL$CaptureSpeed)

plot(datTurnVsStrikeSpeed_LL$Undershoot, datTurnVsStrikeSpeed_LL$DistanceToPrey,col=colourLegL[2],
     ylim=c(0,1),xlim=c(0,2.0),xlab=NA,ylab=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_LL$DistanceToPrey ~ datTurnVsStrikeSpeed_LL$Undershoot)
abline(lFit,col=colourLegL[2],lwd=3.0) ##Fit Line / Regression
mtext(side = 2,cex=cex, line = 2.2, expression("Distance to prey  (mm) " ))
legend("topright",
       legend=paste("LF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ),cex=cex ) 


plot(datTurnVsStrikeSpeed_DL$Undershoot, datTurnVsStrikeSpeed_DL$DistanceToPrey,col=colourLegL[3],
     ylim=c(0,1.0),xlim=c(0,2),   xlab=NA,ylab=NA,main=NA,cex=cex)
lFit <- lm(datTurnVsStrikeSpeed_DL$DistanceToPrey ~ datTurnVsStrikeSpeed_DL$Undershoot)
abline(lFit,col=colourLegL[3],lwd=3.0) ##Fit Line / Regression
mtext(side = 1,cex=cex, line = lineXAxis, expression("Turn ratio ["~gamma~"]" ))
legend("topright",
       legend=paste("DF int.:",prettyNum(digits=3,lFit$coefficients[1])," slope: ",prettyNum(digits=3,lFit$coefficients[2])  ) ,cex=cex) 

dev.off()


### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig5S1-DetectionAngleDensity.pdf",sep=""))
  ### Show Speed Fit ###
  outer = FALSE
  line = 1 ## SubFig Label Params
  lineAxis = 2.7
  lineXAxis = 3.0
  cex = 1.4
  adj  = 3.5
  padj <- -8.0
  las <- 1
  
  ##Margin: (Bottom,Left,Top,Right )
  par(mar = c(4.5,4.3,0.5,1))
  
  
  plot(density(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[1],xlim=c(-120.0,120),lwd=3,lty=1,main=NA,xlab=NA,ylab=NA)
  lines(density(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[2],xlim=c(-120.0,120),lwd=3,lty=2)
  lines(density(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],bw=10),col=colourLegL[3],xlim=c(-120.0,120),lwd=3,lty=3,main=NA)
  
  legend("topleft",
         legend=c(  expression (),
                    bquote(NF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$NL))  ),
                    bquote(LF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$LL))  ),
                    bquote(DF[""] ~ '#' ~ .(NROW(lFirstBoutPoints$DL))  )
         ), 
         col=colourLegL,lty=c(1,2,3,4),lwd=3,cex=cex)
  
  mtext(side = 2,cex=cex, line = lineAxis, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection  " ) )  )
  #mtext("B",at="topleft",outer=outer,side=2,col="black",font=2,las=las,line=line,padj=padj,adj=adj,cex.main=c
dev.off()


### Plot Initial Distance to Prey Vs Final Distance Scatter
plot(lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],lFirstBoutPoints$LL[,"DistanceToPrey"],col=colourLegL[2],pch=pchL[2])
points(lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],lFirstBoutPoints$NL[,"DistanceToPrey"],col=colourLegL[1],pch=pchL[1])
points(lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],lFirstBoutPoints$DL[,"DistanceToPrey"],col=colourLegL[3],pch=pchL[3])
############# 


### Onset/ DETECTION Angle Density supplementary angle figure
pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/fig5S2-DetectionAngleVsDistance_scatter.pdf",sep=""))
  plot(lFirstBoutPoints$LL[,"OnSetAngleToPrey"],lFirstBoutPoints$LL[,"OnSetDistanceToPrey"],col=colourLegL[2],pch=pchL[2],xlim=c(-120.0,120),lwd=2,lty=1,main=NA,xlab=NA,ylab=NA)
  points(lFirstBoutPoints$NL[,"OnSetAngleToPrey"],lFirstBoutPoints$NL[,"OnSetDistanceToPrey"],col=colourLegL[1],pch=pchL[1],lwd=2)
  points(lFirstBoutPoints$DL[,"OnSetAngleToPrey"],lFirstBoutPoints$DL[,"OnSetDistanceToPrey"],col=colourLegL[3],pch=pchL[3],lwd=2)
  mtext(side = 2,cex=cex, line = lineAxis, expression("Prey distance upon detection (mm)") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey azimuth upon detection (deg)  " ) )  )
dev.off()
##

# 
# library(plot3D)
# dev.off()
# layout(matrix(c(1,2,3),1,3, byrow = FALSE))
# 
# 
# pdf(file= paste(strPlotExportPath,"/stat/UndershootAnalysis/stat_UndershootSpeedDistance_3Dmodelplot_AB.pdf",sep=""))
# 
# scatter3D(tail(draw_NF$mu[,1,,],ntail),tail(draw_NF$mu[,3,,],ntail),tail(draw_NF$mu[,2,,],ntail),grid=10,col = colourH[1],
#           ticktype = "detailed",theta=0,phi=0,box=,type= c("shade", "wire", "dots"),zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),xlab="Undershoot", ylab="Distance",zlab="Speed")
# 
# scatter3D(tail(draw_LF$mu[,1,,],ntail),tail(draw_LF$mu[,3,,],ntail),tail(draw_LF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[2],
#           zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)
# 
# scatter3D(tail(draw_DF$mu[,1,,],ntail),tail(draw_DF$mu[,3,,],ntail),tail(draw_DF$mu[,2,,],ntail),ticktype = "detailed",col = colourH[3],
#           zlim =range(draw_LF$mu[,,,]),xlim =c(0.5,1.5),ylim=c(0,0.5),add=T)
# 
# dev.off()
###
#install.packages("plotly")

############# Plot Position Of Prey Prior Capture Bout 

pdf(file= paste(strPlotExportPath,"/PreyPositionPriorCapture_Validated.pdf",sep=""))

plotCaptureBoutPreyPositions()
dev.off()