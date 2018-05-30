### GP Process Estimation Of Hunt Rate Vs Prey Density Using Bayesian Inference Model
myplot_res<- function(ind,qq=0.05){
  
  
  plot(foodlevelsLL,countsLL,col=colourH[1],)
  points(foodlevelsNL,countsNL,col=colourH[2])
  points(foodlevelsDL,countsDL,col=colourH[3])
  
  muLL=apply(drawLL$lambda[,(steps-ind):steps,1],1,mean)
  muNL=apply(drawNL$lambda[,(steps-ind):steps,1],1,mean)
  muDL=apply(drawDL$lambda[,(steps-ind):steps,1],1,mean)
  
  lines(foodlevelsLL,muLL,col=colourH[1],lwd=4)
  lines(foodlevelsNL,muNL,col=colourH[2],lwd=4)
  lines(foodlevelsDL,muDL,col=colourH[3],lwd=4)
  
  band=apply(drawLL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsLL,rev(foodlevelsLL)),c(band[1,],rev(band[2,])),col=colourH[1])
  
  band=apply(drawNL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsNL,rev(foodlevelsNL)),c(band[1,],rev(band[2,])),col=colourH[2])
  
  band=apply(drawDL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
  polygon(c(foodlevelsDL,rev(foodlevelsDL)),c(band[1,],rev(band[2,])),col=colourH[3])
  
}


modelGPV1="model {
  # Likelihood

for(i in 1:N){
  n[i] ~ dnorm(lambda[i],eps)
  #n[i] ~ dpois(lambda[i])
}

eps~dexp(10)
lambda ~ dmnorm(Mu, Sigma.inv)
#n ~ dmnorm(Mu, Sigma.inv)
Sigma.inv <- inverse(Sigma)

# Set up mean and covariance matrix
for(i in 1:N) {
  Mu[i] <- alpha
  Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)

  for(j in (i+1):N) {
    #Sigma[i,j] <- pow(tau,2) * exp( - rho * pow(food[i] - food[j], 2) )
     Sigma[i,j] <- pow(tau,2) * exp( - 0.5*rho^2 * pow(food[i] - food[j], 2) )
    Sigma[j,i] <- Sigma[i,j]
  }
}

alpha ~ dnorm(0,1e-4)T(0,) 
#tau ~ dnorm(tauRange,1e-1)T(0,)
rho = rhoMax

tau  ~ dgamma(tauRange,0.5) 
#rho ~ dunif(0,rhoMax)

}"

##Fixed Rho Model - From Gio 19-05-18 - We fixed rho since estimated both tau and rho made it very sensitive and noisy/non-smooth#
##He uses this to estimate Params for GP, and then constructs a classic GP with the estimated parameters
modelGPFixedRho="model {
  # Likelihood
  
  n ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:N) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
  
    for(j in (i+1):N) {
      Sigma[i,j] <- pow(tau,2) * exp( - 0.5*rho^2 * pow(food[i] - food[j], 2) )
      Sigma[j,i] <- Sigma[i,j]
    }
  }
 
  alpha=0 
  tau0 ~ dgamma(2,0.2) 
  tau  ~ dgamma(2,.2) 
  rho = 0.1
  
}"



#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

load("out/setn-12-D-5-16-datHuntStat.RData")

## Get Event Counts Within Range ##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE) )
datHuntVsPreyL <- rbind(datHuntVsPreyLL,datHuntVsPreyLE)

datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )
datHuntVsPreyNE <- cbind(datHuntStat[,"vHInitialPreyCount"]$NE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NE) )
datHuntVsPreyN <- rbind(datHuntVsPreyNL,datHuntVsPreyNE)

datHuntVsPreyN <- datHuntVsPreyN[!is.na(datHuntVsPreyN[,1]),]


datHuntVsPreyDL <- cbind(datHuntStat[,"vHInitialPreyCount"]$DL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DL) )
datHuntVsPreyDE <- cbind(datHuntStat[,"vHInitialPreyCount"]$DE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$DE) )
datHuntVsPreyD <- rbind(datHuntVsPreyDL,datHuntVsPreyDE)
##Remove NA 
datHuntVsPreyD <- datHuntVsPreyD[!is.na(datHuntVsPreyD[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,100)
colourH <- c(rgb(0.01,0.7,0.01,0.5),rgb(0.9,0.01,0.01,0.5),rgb(0.01,0.01,0.9,0.5),rgb(0.00,0.00,0.0,1.0))

##Thse RC params Work Well to Smooth LF And NF
tauRangeA =40
rhoMaxA = 0.01
Noise = 18 ##The Gaussian Noise Term

burn_in=100;
steps=1000;
thin=1;


##Larva Event Counts Slice
datSliceLL <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2], ]
datSliceNL <- datHuntVsPreyN[datHuntVsPreyN[,1] >= preyCntRange[1] & datHuntVsPreyN[,1] <= preyCntRange[2], ]
datSliceDL <- datHuntVsPreyD[datHuntVsPreyD[,1] >= preyCntRange[1] & datHuntVsPreyD[,1] <= preyCntRange[2], ]
nDatLL <- length(datSliceLL[,1])
nDatNL <- length(datSliceNL[,1])
nDatDL <- length(datSliceDL[,1])


nEventsLL2=as.numeric(datSliceLL[seq(1,nDatLL),2])
nEventsNL2=as.numeric(datSliceNL[seq(1,nDatNL),2])
nEventsDL2=as.numeric(datSliceDL[seq(1,nDatDL),2])

nFoodLL2=as.numeric(datSliceLL[seq(1,nDatLL),1])
nFoodNL2=as.numeric(datSliceNL[seq(1,nDatNL),1])
nFoodDL2=as.numeric(datSliceDL[seq(1,nDatDL),1])


dataLL2=list(n=nDatLL,N=length(nDatLL),food=as.integer(nFoodLL2) );
dataNL2=list(n=nDatNL,N=length(nDatNL),food=as.integer(nFoodNL2) ) ;
dataDL2=list(n=nDatDL,N=length(nDatDL),food=as.integer(nFoodDL2) ) ;


foodlevelsLL=dataLL2[dataLL2$food >= preyCntRange[1] & dataLL2$food <= preyCntRange[2]]$food
ordLL=order(foodlevelsLL)
foodlevelsLL=foodlevelsLL[ordLL]
countsLL=nEventsLL2[ordLL]

foodlevelsNL=dataNL2[dataNL2$food >= preyCntRange[1] & dataNL2$food <= preyCntRange[2]]$food
ordNL=order(foodlevelsNL)
foodlevelsNL=foodlevelsNL[ordNL]
countsNL=nEventsNL2[ordNL]

foodlevelsDL=dataDL2[dataDL2$food >= preyCntRange[1] & dataDL2$food <= preyCntRange[2]]$food
ordDL=order(foodlevelsDL)
foodlevelsDL=foodlevelsDL[ordDL]
countsDL=nEventsDL2[ordDL]


dataLL=list(n=countsLL,food=foodlevelsLL,N=length(countsLL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataNL=list(n=countsNL,food=foodlevelsNL,N=length(countsNL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);
dataDL=list(n=countsDL,food=foodlevelsDL,N=length(countsDL),tauRange=tauRangeA,rhoMax=rhoMaxA,tau0=Noise);

varnames=c("tau","rho","alpha","lambda")


library(rjags)
fileConn=file("model.tmp")
writeLines(modelGPV1,fileConn);
close(fileConn)

mLL=jags.model(file="model.tmp",data=dataLL);
mNL=jags.model(file="model.tmp",data=dataNL);
mDL=jags.model(file="model.tmp",data=dataDL);
update(mLL,burn_in);update(mNL,burn_in);update(mDL,burn_in)


drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)

strPlotName <- paste("plots/stat_HuntEventRateVsPrey_GPEstimate-tauMax",tauRangeA,".pdf",sep="-")
pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 

X11()
myplot_res(100)
  
dev.off()

