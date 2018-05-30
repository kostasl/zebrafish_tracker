model="model {
  # Likelihood
  
  for(i in 1:N){
	 n[i] ~ dnorm(lambda[i],eps)
  }

  eps~dexp(10)
  lambda ~ dmnorm(Mu, Sigma.inv)
  #n ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:N) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(tau, 2)+0.01
  
    for(j in (i+1):N) {
      Sigma[i,j] <- pow(tau,2) * exp( - rho * pow(food[i] - food[j], 2) )
      Sigma[j,i] <- Sigma[i,j]
    }
  }
 
  alpha ~ dnorm(0,1e-4)T(0,) 
  tau ~ dnorm(tauRange,1e-1)T(0,)
  rho ~ dunif(0,rhoMax)
  
}"

##Thse RC params Work Well to Smooth LF And NF
tauRangeA =10000
rhoMaxA = 0.01

burn_in=100;
steps=1000;
thin=1;

#nLL1=read.table("Stats/mcmc/testLL.dat")$Freq
#nNL1=read.table("Stats/mcmc/testNL.dat")$Freq
#nDL1=read.table("Stats/mcmc/testDL.dat")$Freq

load("out/setn-12-D-5-16-datHuntStat.RData")

## Get Event Counts Within Range ##
datHuntVsPreyLL <- cbind(datHuntStat[,"vHInitialPreyCount"]$LL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL) )
datHuntVsPreyLE <- cbind(datHuntStat[,"vHInitialPreyCount"]$LE , as.numeric(datHuntStat[,"vHLarvaEventCount"]$LE) )
datHuntVsPreyL <- rbind(datHuntVsPreyLL,datHuntVsPreyLE)

datHuntVsPreyL <- datHuntVsPreyL[!is.na(datHuntVsPreyL[,1]),]


### Cut And Examine The data Where There Are Between L and M rotifers Initially
preyCntRange <- c(0,120)
colourH <- c(rgb(0.01,0.7,0.01,0.5),rgb(0.9,0.01,0.01,0.5),rgb(0.01,0.01,0.9,0.5),rgb(0.00,0.00,0.0,1.0))

##Larva Event Counts Slice
datSliceLL <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2],2 ]
nDatLL <- length(datSliceLL)

nEventsLL2=as.numeric(datSliceLL[seq(1,nDatLL,1)])

dataLL2=list(n=nDatLL,NTOT=length(nDatLL),food=as.integer(datHuntVsPreyL[,1]));


foodlevelsLL=dataLL2[dataLL2$food >= preyCntRange[1] & dataLL2$food <= preyCntRange[2]]$food
ordLL=order(foodlevelsLL)
foodlevelsLL=foodlevelsLL[ordLL]
countsLL=nEventsLL2[ordLL]


dataLL=list(n=countsLL,food=foodlevelsLL,N=length(countsLL),tauRange=tauRangeA,rhoMax=rhoMaxA);

varnames=c("tau","eps","rho","alpha","lambda","tauRange","rhoMax")


library(rjags)
fileConn=file("model.tmp")
writeLines(model,fileConn);
close(fileConn)

mLL=jags.model(file="model.tmp",data=dataLL);
update(mLL,burn_in);

drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
ind = 100
plot_res<- function(ind,qq=0.05){

  X11()
	plot(foodlevelsLL,countsLL,col=colourH[1])

	muLL=apply(drawLL$lambda[,(steps-ind):steps,1],1,mean)

	lines(foodlevelsLL,muLL,col=colourH[1],lwd=4)

	band=apply(drawLL$lambda[,(steps-ind):steps,1],1,quantile,probs=c(qq,1-qq))
	polygon(c(foodlevelsLL,rev(foodlevelsLL)),c(band[1,],rev(band[2,])),col=colourH[1])
}

