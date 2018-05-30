### Do Regrasssion of HuntRates Vs Initial Prey using Gaussian Process
## The parameters of GP Covariance matrix are estimated using Bayesian inference of the same GP and then these are used
## to Calculate the GP (Inverse the Matrix) and Plot

SE <- function(Xi,Xj, rho,tau) tau^2*exp(-0.5*(Xi - Xj) ^ 2 * rho^2)
covC <- function(X, Y, rho,tau) outer(X, Y, SE, rho,tau)

plot_res<- function(ind,drawY,Xn,Yn,colour='red ',qq=0.05){
  
  ord=order(Xn)
  Xn=Xn[ord] ##Place Points In order so we can draw the Polygon Bands
  Yn=Yn[ord]
  
  points(Xn,Yn,col=colour)
  x_predict=seq(1,80,0.1)
  Ef=matrix(NA,ncol=length(x_predict),nrow=ind)
  for(j in 1:ind){
    i=steps-j+1
    cov_xx_inv=solve(covC(Xn,Xn,drawY$rho[1,i,1],drawY$tau[1,i,1])+diag(drawY$tau0[1,i,1],length(Xn)))
    Ef[j,] <- covC(x_predict, Xn,drawY$rho[1,i,1],drawY$tau[1,i,1]) %*% cov_xx_inv %*% Yn 
  }
  mu=apply(Ef,2,mean)
  band=apply(Ef,2,quantile,probs=c(qq,1-qq))
  lines(x_predict,mu,lwd=2,col=colour,xlim=c(0,max(x_predict) ) )
  polygon(c(x_predict,rev(x_predict)),c(band[1,],rev(band[2,])),col=colour)
}




model="model {
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
  tau0 ~ dgamma(40,0.5) 
  tau  ~ dgamma(40,0.5) 
  rho = 0.01
  
}"

library(rjags)
fileConn=file("model.tmp")
writeLines(model,fileConn);
close(fileConn)

burn_in=100;
steps=1000;
thin=1;
varnames=c("tau","rho","alpha","tau0")

#load("data/setn-12-D-5-16-datHuntStat.RData")

load("./Stats/data/setn-12-D-5-16-datHuntStat.RData")

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

##Larva Event Counts Slice
datSliceLL <- datHuntVsPreyL[datHuntVsPreyL[,1] >= preyCntRange[1] & datHuntVsPreyL[,1] <= preyCntRange[2], ]
datSliceNL <- datHuntVsPreyN[datHuntVsPreyN[,1] >= preyCntRange[1] & datHuntVsPreyN[,1] <= preyCntRange[2], ]
datSliceDL <- datHuntVsPreyD[datHuntVsPreyD[,1] >= preyCntRange[1] & datHuntVsPreyD[,1] <= preyCntRange[2], ]
nDatLL <- length(datSliceLL[,1])
nDatNL <- length(datSliceNL[,1])
nDatDL <- length(datSliceDL[,1])

nFoodLL2=as.numeric(datSliceLL[,1]); nEventsLL2=as.numeric(datSliceLL[,2]); 
nFoodNL2=as.numeric(datSliceNL[,1]); nEventsNL2=as.numeric(datSliceNL[,2]);
nFoodDL2=as.numeric(datSliceDL[,1]); nEventsDL2=as.numeric(datSliceDL[,2]);


m<-list() ##model List
draw<-list()

##LL

#ord=order(nFoodLL2)
#nFoodLL2=nFoodLL2[ord]
#nEventsLL2=nEventsLL2[ord]
data=list(n=nEventsLL2,food=nFoodLL2,N=length(nEventsLL2));


m[["LL"]]=jags.model(file="model.tmp",data=data);
update(m[["LL"]],burn_in)
draw[["LL"]]=jags.samples(m[["LL"]],steps,thin=thin,variable.names=varnames)



##NL
#ord=order(nFoodNL2)
#nFoodNL2=nFoodNL2[ord]
#nEventsNL2=nEventsNL2[ord]
data=list(n=nEventsNL2,food=nFoodNL2,N=length(nEventsNL2));


m[["NL"]]=jags.model(file="model.tmp",data=data);
update(m[["NL"]],burn_in)
draw[["NL"]]=jags.samples(m[["NL"]],steps,thin=thin,variable.names=varnames)

##DL
#ord=order(nFoodDL2)
#nFoodDL2=nFoodDL2[ord]
#nEventsDL2=nEventsDL2[ord]
data=list(n=nEventsDL2,food=nFoodDL2,N=length(nEventsDL2));

m[["DL"]]=jags.model(file="model.tmp",data=data);
update(m[["DL"]],burn_in)
draw[["DL"]]=jags.samples(m[["DL"]],steps,thin=thin,variable.names=varnames)

ind = 10

colourH <- c(rgb(0.01,0.7,0.01,0.5),rgb(0.9,0.01,0.01,0.5),rgb(0.01,0.01,0.9,0.5),rgb(0.00,0.00,0.0,1.0))


strPlotName <- paste("plots/stat_HuntEventRateVsPrey_GPEstimate-tauLL",round(mean(draw[["LL"]]$tau)),".pdf",sep="-")
pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 

plot(nFoodLL2,nEventsLL2,col=colourH[1],main = "GP Regression Of HuntRate Vs Initial Prey Count ",
     ylab="Number of Hunt Events",xlab="Initial Tracker-Estimated Prey Count",
     sub=paste("GP tau:",format(mean(draw[["LL"]]$tau),digits=4 ),"tau0:",format(mean(draw[["LL"]]$tau0),digits=4 ) ,"rho:",format(mean(draw[["LL"]]$rho),digits=4 ) )  
               )

legend(35,75,legend = c(paste("LL #",nDatLL),paste("NL #",nDatNL),paste("DL #",nDatDL)),fill=colourH)


plot_res(ind,draw[["LL"]],nFoodLL2,nEventsLL2,colourH[1],0.05)

plot_res(ind,draw[["NL"]],nFoodNL2,nEventsNL2,colourH[2],0.05)

#plot(nFoodDL2,nEventsDL2,col="blue")

plot_res(ind,draw[["DL"]],nFoodDL2,nEventsDL2,colourH[3],0.05)

dev.off()
#points(datSliceLL[,1],datSliceLL[,2],col="black")


