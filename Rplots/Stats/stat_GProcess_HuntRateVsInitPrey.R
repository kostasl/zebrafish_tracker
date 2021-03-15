### Do Regresssion of HuntRates Vs Initial Prey using Gaussian Process
## The parameters of GP Covariance matrix are estimated using Bayesian inference of the same GP and then these are used
## to Calculate the GP (Inverse the Matrix) and Plot
## Note this File Is Included by R NoteBook ForagingState Analysis_huntEvents.Rmd

source("HuntingEventAnalysis_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")


model="model {
  # Likelihood
  
  n ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:N) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
  
    for(j in (i+1):N) {
      #Sigma[i,j] <- pow(tau,2) * exp( - 0.5* pow((food[i] - food[j])*rho, 2) )
      ##exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
      Sigma[i,j] <-   pow(tau,2)*exp( - pow(rho*(food[i] - food[j]), 2) )
        
      Sigma[j,i] <- Sigma[i,j]
    }
  }
 
  alpha=0 
  tau0 ~ dgamma(250,1) 
  tau  ~ dgamma(250,1) 
  rho ~  dgamma(1,1/3) 
}"

library(rjags)
fileConn=file("model.tmp")
writeLines(model,fileConn);
close(fileConn)


## Prepare Data - 
preyCntRange <- c(0,80) ## Prey Density Range to Include in Model
#### LOAD And prepare MODEL DATA #### 
datHuntLabelledEventsKL <- getLabelledHuntEventsSet() # readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
##Clear Warningss : assign("last.warning", NULL, envir = baseenv()
## Link Evoked and Spontaneous Trajectories
datHuntStat <- makeHuntStat(datHuntLabelledEventsKL)

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


#### Run Model  #### 
# Saves Sampls to RData file and returns Samples and Data - so they can be plotted
inferGPModel_HuntRateVsPreyDensity <- function (burn_in=140,steps=10000,thin=2)
{
  varnames=c("tau","rho","alpha","tau0")
  m<-list() ##model List
  draw<-list()
  data <- list()
  ##LL
  #ord=order(nFoodLL2)
  #nFoodLL2=nFoodLL2[ord]
  #nEventsLL2=nEventsLL2[ord]
  data[["LF"]]=list(n=nEventsLL2,food=nFoodLL2,N=length(nEventsLL2));
  
  m[["LF"]]=jags.model(file="model.tmp",data=data[["LF"]]);
  update(m[["LF"]],burn_in)
  draw[["LF"]]=jags.samples(m[["LF"]],steps,thin=thin,variable.names=varnames)
  ##NL
  #ord=order(nFoodNL2)
  #nFoodNL2=nFoodNL2[ord]
  #nEventsNL2=nEventsNL2[ord]
  data[["NF"]]=list(n=nEventsNL2,food=nFoodNL2,N=length(nEventsNL2));
  
  m[["NF"]]=jags.model(file="model.tmp",data=data[["NF"]]);
  update(m[["NF"]],burn_in)
  draw[["NF"]]=jags.samples(m[["NF"]],steps,thin=thin,variable.names=varnames)
  
  ##DL
  #ord=order(nFoodDL2)
  #nFoodDL2=nFoodDL2[ord]
  #nEventsDL2=nEventsDL2[ord]
  data[["DF"]]=list(n=nEventsDL2,food=nFoodDL2,N=length(nEventsDL2));
  
  m[["DF"]]=jags.model(file="model.tmp",data=data[["DF"]]);
  update(m[["DF"]],burn_in)
  draw[["DF"]]=jags.samples(m[["DF"]],steps,thin=thin,variable.names=varnames)
  
  save(draw,data,m,steps,thin,
       file=paste0(strDataExportDir,"/jags_FoodDensityVsHuntRate_GP2.RData"))
  
  return (list(draw,data))
}


#### Plot Model ####
SE <- function(Xi,Xj, rho,tau) tau^2*exp(-(Xi - Xj) ^ 2 * rho^2)
covC <- function(X, Y, rho,tau) outer(X, Y, SE, rho,tau)

plot_res<- function(ind,drawY,Xn,Yn,colour='red ',qq=0.05,pPch=16){
  
  ord=order(Xn)
  Xn=Xn[ord] ##Place Points In order so we can draw the Polygon Bands
  Yn=Yn[ord]
  
  points(Xn,Yn,col=colour, pch=pPch)
  x_predict=seq(preyCntRange[1],preyCntRange[2],1)
  Ef=matrix(NA,ncol=length(x_predict),nrow=ind)
  for(j in 1:ind){
    #i=steps/thin-j+1
    i=length(drawY$tau[]) - j+1
    ##print(i)
    cov_xx_inv=solve(covC(Xn,Xn,drawY$rho[i],drawY$tau[i])+diag(drawY$tau0[i],length(Xn)))
    Ef[j,] <- covC(x_predict, Xn,drawY$rho[i],drawY$tau[i]) %*% cov_xx_inv %*% Yn 
  }
  mu=apply(Ef,2,mean)
  sd=apply(Ef,2,sd)
  
  #band=apply(Ef,2,quantile,probs=c(qq,1-qq))
  band1= mu + 2*sd
  band2= mu - 2*sd
  lines(x_predict,mu,lwd=3,col=colour,xlim=c(0,max(x_predict) ) )
  #polygon(c(x_predict,rev(x_predict)),c(band[1,],rev(band[2,])),col=colour)
  polygon(c(x_predict,rev(x_predict)),c(band1,rev(band2)),col=colour)
}



colourH <- c(rgb(0.01,0.7,0.01,0.5),rgb(0.9,0.01,0.01,0.5),rgb(0.01,0.01,0.9,0.5),rgb(0.00,0.00,0.0,1.0))


ind = 10

#
load(file=paste0(strDataExportDir,"/jags_FoodDensityVsHuntRate_GP2.RData"))
Rho <-format( mean(draw$LF$rho),digits=2)
tau <- format( mean(draw$LF$tau),digits=2)
#strPlotName <- paste("plots/stat_HuntEventRateVsPrey_GPEstimate-tauLL",round(mean(draw[["LL"]]$tau)),".pdf",sep="-")
strPlotName <-  paste(strPlotExportPath,"/stat_HuntEventRateVsPrey_GioGPEstimate-tauMax",tau,"-Rho",Rho,".pdf",sep="")
pdf(strPlotName,width=8,height=8,title="GP Function of Hunt Rate Vs Prey") 
par(mar = c(4.1,4.8,3,1))
  
  plot(data$LF$food,data$LF$n,col=colourH[1],
       main = NA,
       ylab="Number of Hunt Events in 10 min",
       xlab="Initial Prey Density (Rotifers/10ml)",
       cex=1.4,
       cex.axis = 1.7,
       cex.lab = 1.7,
       xlim = c(1,60),##preyCntRange,
       #log="x",
       pch=pointTypeScheme$LL,
       #sub=paste("GP tau:",format(mean(draw[["LF"]]$tau),digits=4 ),
      #           "tau0:",format(mean(draw[["LF"]]$tau0),digits=4 ) ,
      #           "rho:",format(mean(draw[["LF"]]$rho),digits=4 ) )  
  )
  
  legend("topright",legend = c(paste("LF #",data$LF$N),paste("NF #",data$NF$N ),paste("DF #",data$DF$N)),
         col=c(colourDataScheme[["LF"]]$Evoked,colourDataScheme[["NF"]]$Evoked,colourDataScheme[["DF"]]$Evoked),
         pch=c(pointTypeScheme$LL,pointTypeScheme$NL,pointTypeScheme$DL ),cex=1.5 )
  
  
  plot_res(ind,draw[["LF"]],data$LF$food,data$LF$n, colourH[1],0.05,pointTypeScheme$LL)
  
  plot_res(ind,draw[["NF"]],data$NF$food,data$NF$n,colourH[2],0.05,pointTypeScheme$NL)
  
  #plot(nFoodDL2,nEventsDL2,col="blue")
  
  plot_res(ind,draw[["DF"]],data$DF$food,data$DF$n,colourH[3],0.05,pointTypeScheme$DL)
  
  #legend(5,700,legend = c(paste("LL #",nDatLL),paste("NL #",nDatNL),paste("DL #",nDatDL)),fill=colourH)

dev.off()
#points(datSliceLL[,1],datSliceLL[,2],col="black")




## Plot - Compare initial Prey Density Between Rearing Groups experiments ###
strCumPlotName <-  paste(strPlotExportPath,"/fig2S2-InitPreyCount_CDF.pdf",sep="")
pdf(strCumPlotName,width=8,height=8,title="Compare prey density testing conditions between groups") 
  
  par(mar = c(3.9,4.7,2,1))
  plot(ecdf(foodlevelsNL),xlim=c(0,60),lwd=4,lty=1,col=colourLegL[1],main=NA,xlab=NA,ylab=NA,cex=cex,cex.axis=cex,pch=pchL[1])
  lines(ecdf(foodlevelsLL),xlim=c(0,60),lwd=4,lty=2,pch=pchL[2],col=colourLegL[2],cex=cex)
  lines(ecdf(foodlevelsDL),xlim=c(0,60),lwd=4,lty=3,pch=pchL[3],col=colourLegL[3],cex=cex)
  mtext(side = 1,cex=cex, line = 2.7, expression("Initial prey count in ROI (Tracker estimate)" ))
  mtext(side = 2,cex=cex, line = 2.2, expression(" Cumulative distribution " ))
  
  legend("bottomright",pch=pchL,cex=cex,
         legend = c(paste("NF #",nDatNL),paste("LF #",nDatLL),paste("DF #",nDatDL)),col=colourLegL)

dev.off()


## Do Significance tests - Pairwise - 
preyLevelsPerGroup <-( rbind(cbind(as.integer(datHuntStat[,"vHInitialPreyCount"]$LL),"LF"),
                             cbind(as.integer(datHuntStat[,"vHInitialPreyCount"]$DL),"DF"),
                             cbind(as.integer(datHuntStat[,"vHInitialPreyCount"]$NL),"NF")
                             , cbind(as.integer(rbinom(60, 20, .5)),"TF") ## Just Ranodm Binomial
))

preyLevelsPerGroup <- data.frame(preyCount=as.integer(preyLevelsPerGroup[,1]),group=as.factor(preyLevelsPerGroup[,2]) )
pairwise.t.test (preyLevelsPerGroup$preyCount ,preyLevelsPerGroup$group, pool.sd = TRUE,paired=FALSE)
### Added NS diffe in prey ouunt to Fig Supp
