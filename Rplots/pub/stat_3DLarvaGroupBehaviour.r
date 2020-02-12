### Kostas Lagogiannis 2019-06-24 
## Learning steers the ontogeny of an efficient hunting sequence in zebrafish larvae.

## 3D Multivariate Gaussian Model capturing group behaviour
## ******** No clustering Between slow and Fast Swims ****** 
## I made this to complement the Clustering Method, so as to characterize the overall covariance structure


library(rjags)
library(here) ###Used ti Set correct paths

#library(runjags)
setwd(here())
#set_here(path="./pub") 
source("common_lib.R")
#source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
#source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
#source("HuntingEventAnalysis_lib.r")


initfunct <- function(nchains,N)
{
  initlist <- replicate(nchains,list(
                                     ".RNG.name"="base::Super-Duper",
                                     ".RNG.seed"=round(runif(1,0,60000)) ),
                                     simplify=FALSE)
  return(initlist)
}



## Plot Function for Covariance Density of Mean Group Behaviour, given by averaging the covariance of individual larvae
plotModelCovCoeff <- function(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail)
{
  
  ##Covariance 
  nsam <- NROW(draw_LF$muG[,3,,1])
  ylimR <- c(0,4)
  cov_LF <- (draw_LF$cov[,Ci,Cj,(nsam-ntail):nsam,]/sqrt(draw_LF$cov[,Ci,Ci,(nsam-ntail):nsam,]*draw_LF$cov[,Cj,Cj,(nsam-ntail):nsam,]) )
  cov_NF <- (draw_NF$cov[,Ci,Cj,(nsam-ntail):nsam,]/sqrt(draw_NF$cov[,Ci,Ci,(nsam-ntail):nsam,]*draw_NF$cov[,Cj,Cj,(nsam-ntail):nsam,]) )
  cov_DF <- (draw_DF$cov[,Ci,Cj,(nsam-ntail):nsam,]/sqrt(draw_DF$cov[,Ci,Ci,(nsam-ntail):nsam,]*draw_DF$cov[,Cj,Cj,(nsam-ntail):nsam,]) )
  
  ##Average Over Columns/Samples - Produce Distribution of Mean Cov Of Group -  Across Samples (Vector With n mean points)
  ##What is the Members Cov On Average?Distibution of E[rho_g] = Sum(rho_i,NLarvae)/NLarvae
  plot( density( apply(cov_LF,2,"mean"),
                 from=-1,to=1,n=200,bw=0.1),xlim=c(-0.5,0.5) ,col=colourLegL[2],lwd=3,main="",xlab="",ylab="",ylim=ylimR,lty=1)
  lines(density( apply(cov_NF,2,"mean"),
                 from=-1,to=1,n=200,bw=0.1),xlim=c(-0.5,0.5),col=colourLegL[1],lwd=3,lty=2 )
  lines(density( apply(cov_DF,2,"mean"),
                 from=-1,to=1,n=200,bw=0.1),xlim=c(-0.5,0.5),col=colourLegL[3],lwd=3,lty=3 )
}


## Returns estaimes of each larvae behaviour According the Model
getEstimatesPerLarva <- function(drawG,stail)
{
  ldist <- list()
  lSpeed <- list()
  lTurnRatio <- list()
  ##Iterate Through each Larva and get a mean estimate of behaviour according to model
  for ( i in (1:head(as.numeric(drawG$NLarv),1)) )
  {
    lTurnRatio[[i]]  <- sapply(tail(drawG$mu[i,1,,],stail),mean)
    ldist[[i]]       <- sapply(tail(drawG$mu[i,3,,],stail),mean)
    lSpeed[[i]]      <- sapply(tail(drawG$mu[i,2,,],stail),mean)
  }
  ##Overlay The Density From The Estimated Mean Overshoot Of Each Larva
  mEstDistToPrey <- unlist(lapply(ldist,mean) )
  mEstSpeed      <- unlist(lapply(lSpeed,mean) ) 
  mEstTurnRatio  <- unlist(lapply(lTurnRatio,mean) ) 
  return(cbind(Undershoot=mEstTurnRatio, DistanceToPrey=mEstDistToPrey,CaptureSpeed=mEstSpeed))  
  
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
  muG[g,1] ~ dnorm(1, 8)T(0.0,2) ##turn ratio
  muG[g,2] ~ dnorm(25,0.001)T(0,) ##cap speed
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

##
##
steps <-15000
nchains <- 7
nthin <- 5
#str_vars <- c("mu","rho","sigma","x_rand") #Basic model 
str_vars <- c("mu","cov","x_rand","muG","tG","NLarv","Lid") #Mixture Model

getwd()
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
##Contains Serial Larvae ID, that links each hunt event to an individual larva 
ldata_LF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_LF.rds") )
ldata_NF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_NF.rds") )
ldata_DF  <-  readRDS(file=paste0("dat/huntEpisodeDataMergedWithLarvalSuccess_DF.rds") )
datHuntLarvaStat <- readRDS(file=paste0("dat/LarvaEmpiricalMeanHuntBehaviour.rds"))

### RUN JAGS MODEL ###
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
####### END OF RUN MODELS ##

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



schain <- 5:10 ## Chains used for data visualiazation
stail <- 300 ## Number Of  Chain Samples to Use for Plots - from the end of the chain


##CONVERGENCE - CHECK
##Each of draw_XX arrays is structured  as
# Idx: Lid,Variable (1 turnratio,2-speed,3-distance),Sample Row,Chain
plot(density(tail(draw_LF$muG[,3,,1],1000) ),type='l',main="Distance")
for (c in schain)
  lines(density(tail(draw_LF$muG[,3,,c],1000) ),col="red")

plot(density(tail(draw_LF$muG[,2,,1],1000) ),type='l',main="Speed")
for (c in schain)
  lines(density(tail(draw_LF$muG[,2,,c],1000) ),col="red")


################################## ###################
### 3D Display of Group Behaviour MODEL            ##
###################################################
  library( rgl )
  ntail <- 30
  
  ##Prepare Data
  datMu3D <-  data.frame( cbind.data.frame(
    TurnR=tail(draw_LF$muG[,1,,1], ntail) ,
    CSpeed=tail(draw_LF$muG[,2,,1], ntail),
    Dist=tail(draw_LF$muG[,3,,1], ntail), col=colourLegL[2],pch=pchL[4],group="LF" )  )
  
  datMu3D <- rbind(datMu3D,
                   data.frame( cbind.data.frame(
                     TurnR=tail(draw_NF$muG[,1,,1], ntail) ,
                     CSpeed=tail(draw_NF$muG[,2,,1], ntail),
                     Dist=tail(draw_NF$muG[,3,,1], ntail), col=colourLegL[1],pch=pchL[6]),group="NF"  )
  )
  
  datMu3D <- rbind(datMu3D,
                   data.frame( cbind.data.frame(
                     TurnR=tail(draw_DF$muG[,1,,1], ntail) ,
                     CSpeed=tail(draw_DF$muG[,2,,1], ntail),
                     Dist=tail(draw_DF$muG[,3,,1], ntail), col=colourLegL[3],pch=pchL[5],group="DF")  )
  )
  datMu3D$col  <- as.character( datMu3D$col)
  
  ##Open Window And Plot
  open3d()
  bbox <- par3d('bbox') 
  rgl::plot3d( x=datMu3D$TurnR, y=datMu3D$CSpeed, z=datMu3D$Dist, col = datMu3D$col, type = "s", radius = 1.3,
               #xlab="Turn Ratio", ylab="Capture Speed (mm/sec)",zlab="Distance to prey (mm)",
               xlab="Turn Ratio", ylab="Speed",zlab="Distance",
               xlim=c(0.5,1.5), ylim=c(10,50), zlim=c(0,0.8),
               box = TRUE ,aspect = TRUE,axes=FALSE
               #,expand = 1.5
  )
  box3d()
  title3d(main=NULL)
  rgl::axis3d('x+-',at=c(0.5,0.8,1,1.2,1.5))
  rgl::axis3d('z-+',at=seq(0.0,0.7,len=8))
  rgl::axis3d('y+-',at=seq(50,10,len=5),labels=rev(seq(10,50,len=5))) 
  
  
  rgl::rgl.viewpoint(0,-60,fov=35,type = c("userviewpoint") )
  rgl::rgl.viewpoint(0,0)
  


#################################################
## Speed Posterior Group Vs Individual Density
##############################################
## Live Fed
with(draw_LF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Compare Data to Model - Speed LF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==2,]$CaptureSpeed,bw=2),col="blue",lwd=3,lty=2)

##Dry Fed DF
##Speed Posterior Group Vs INdividual Density
with(draw_DF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Compare Data to Model  DF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
   lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail)) )) ,col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==1,]$CaptureSpeed,bw=2),col="blue",lwd=2,lty=2)

##Not Fed
##Speed Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,2,,], 150) ),ylim=c(0,1),xlim=c(0,60),lwd=2,col="red",main="Speed NF")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,2,,],stail)),lty=2)
  ###Show Inferred Distribution
  lines(1:100,dnorm(1:100,mean=mean( tail(muG[,2,,3],stail)),sd=sqrt(mean( tail( 1/tG[,2,,1],stail))) ),col="purple",lwd=4)
})
##Compare To Empirical - Change group DF,LF,NF-- V Good Match!
lines(density(datHuntLarvaStat[datHuntLarvaStat$groupID==3,]$CaptureSpeed,bw=2),col="blue",lwd=2,lty=2)


##########
##Distance Posterior Group Vs INdividual Density
with(draw_NF,{
  plot(density(tail(muG[,3,,1], 100) ),ylim=c(0,16),xlim=c(0,1),lwd=2,col="red")
  for ( i in (1:NLarv[1] ) )
    lines( density( tail( mu[i,3,,1],stail)),lty=2)
})

### Compare Group Model To Density Obtain through Mean Estimated Behaviour For Each Larva
### Obtain Estimated Mean Values For Each Larva & Plot Group Population
## Plot Distance Density
plot(density(sapply(tail(draw_LF$mu[,3,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Distance to Prey",ylim=c(0,6)) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_LF <- getEstimatesPerLarva(draw_LF,stail)
lines(density( unlist(lapply(lModelEst_LF[,"DistanceToPrey"],mean) ) ) )

lines(density(sapply(tail(draw_NF$mu[,3,,],stail),mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_NF <- getEstimatesPerLarva(draw_NF,stail)
lines(density( unlist(lapply(lModelEst_NF[,"DistanceToPrey"],mean) ) ) )

lines(density(sapply(tail(draw_DF$mu[1,3,,],stail) ,mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lModelEst_DF <- getEstimatesPerLarva(draw_DF,stail)
lines(density( unlist(lapply(lModelEst_DF[,"DistanceToPrey"],mean) ) ) )


## Plot Speed Density
plot(density(sapply(tail(draw_LF$mu[,2,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Capture Speed",ylim=c(0,0.1)) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,2,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,2,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

## Plot Undershoot Density / Mean Sample point Across larva 
plot(density(sapply(tail(draw_LF$mu[,1,,],stail),mean)),col=colourLegL[2] ,lwd=2,main="Turn ratio") ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_NF$mu[,1,,],stail),mean)),col=colourLegL[3] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva
lines(density(sapply(tail(draw_DF$mu[,1,,],stail) ,mean)),col=colourLegL[1] ,lwd=2) ##Mean Group Undershoot From Mean Of Each Larva

##"Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey"
  Ci <- 2
  Cj <- 3
  plotModelCovCoeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail)
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Prey distance to capture speed covariance coeff." ) )  )


#"Covariance in 3D statistical model for Capture Strike speed / Undershoot Ratio / Distance to Prey")
  Ci <- 1
  Cj <- 3
  plotModelCovCoeff(Ci,Cj,draw_LF,draw_NF,draw_DF,ntail)
  mtext(side = 2,cex=cex, line = lineAxis,padj=padj, expression("Density function") )
  mtext(side = 1,cex=cex, line = lineAxis, expression(paste("Turn-ratio to prey distance covariance coeff." ) )  )

