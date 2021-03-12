## \author kostas lagogiannsi 2021
##  Estimates the hidden function Mean squeader Displacement vs prey Density for each group - Uses Non-parametric Gaussian Process with Bayesian Inference
## Runs Parallel Models Exploring Prior Parameter Space
library(rjags)

library(foreach)
library(doParallel)

source("HuntingEventAnalysis_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")

# Make Model with Parametrizable with Fixed Rho
# \return FileName Where New Model Was Saved
model <- function(tauShape,tauRate,const_rho) 
{
  ## Show A little sample of what the correlation function looks like
  plot(  ( rgamma(1,shape=tauShape,rate=tauRate)^2) * exp( - (const_rho*(seq(-50,50,1)) )^2 )  , 
         log="y",ylim=c(0.01,100000) ) 
  
  strMdl <- paste("model {
    # Likelihood
    
    MSD ~ dmnorm(Mu, Sigma.inv)
    Sigma.inv <- inverse(Sigma)
    
    # Set up mean and covariance matrix
    for(i in 1:N) {
      Mu[i] <- alpha
      Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
    
      for(j in (i+1):N) {
        #Sigma[i,j] <- pow(tau,2) * exp( - 0.5* pow((prey[i] - prey[j])*rho, 2) )
        ##exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
        
        Sigma[i,j] <-   pow(tau,2) * exp( - pow(rho*(prey[i] - prey[j]), 2) )
        Sigma[j,i] <- Sigma[i,j]
      }
    }
   
    alpha=0 
    rho =  ",const_rho," #dgamma(1,1)  ## Covariance Width / Links neighbouring Densities Between 
    tau0 ~ dgamma(",tauShape,",",tauRate,") ##Covariance Peak Strength
    tau  ~ dgamma(",tauShape,",",tauRate,")  
    
  }")
  
  ## Save to File And Return File Name
  modelFileName <-paste0("model-tauS",tauShape,"R",tauRate,"-rho",const_rho,".tmp")
  
  fileConn=file(modelFileName)
  writeLines(strMdl,fileConn);
  close(fileConn)
  
  return(modelFileName)
}



#### Run Model  #### 
# Saves Sampls to RData file and returns Samples and Data - so they can be plotted
inferGPModel_MSDVsPreyDensity <- function (burn_in=140,steps=10000,dataSamples=120,thin=2,modelFileName)
{
  library(rjags)
  
  modelData <- list()
  m<-list() ##model List
  draw<-list()
  
  strRG <- list(LF=c("LE","LL"),NF=c("NE","NL"),DF=c("DE","DL")) #LF=c("LE","LL")
  for (strG in names(strRG) )
  {
    message(strG) 
    ## INcrease Sampling of the more rare High Density  
    datDispersionG_Low <- datDispersion[!is.na(datDispersion$PreyCount) &
                                          !is.na(datDispersion$MSD) &
                                          datDispersion$groupID %in% strRG[[strG]] &
                                          datDispersion$PreyCount <= preyCntRange[2]/3, ]
    idx_L = sample(1:NROW(datDispersionG_Low),
                   size=min(dataSamples/2,NROW(datDispersionG_Low)), replace=FALSE )
    
    datDispersionG_High <- datDispersion[!is.na(datDispersion$PreyCount) &
                                           !is.na(datDispersion$MSD) &
                                           datDispersion$groupID %in% strRG[[strG]] &
                                           datDispersion$PreyCount > preyCntRange[2]/3 & 
                                           datDispersion$PreyCount <= preyCntRange[2], ]
    idx_H = sample(1:NROW(datDispersionG_High),
                   size=min(dataSamples/2,NROW(datDispersionG_High)), replace=FALSE )
    
    ##Combine Low/High Range Samples
    ## Sample MSD and PreyDensity From Each Group equally at both ranges (Reduce error in Estimate High Range)
    datDispersion_Subset <- rbind(datDispersionG_High[idx_H,],datDispersionG_Low[idx_L,])
    ## Correct Prey Count Noise in Empty Condition
    datDispersion_Subset[datDispersion_Subset$groupID == strRG[[strG]][1],"PreyCount" ] = 0
    
    # Show Data Points
    plot(datDispersion_Subset$PreyCount,datDispersion_Subset$MSD)
    
    #### LOAD And prepare MODEL DATA #### 
    assign("last.warning", NULL, envir = baseenv())
    
    modelData[[strG]] = list(prey=datDispersion_Subset$PreyCount, 
                             MSD=as.numeric(datDispersion_Subset$MSD),
                             N = NROW(datDispersion_Subset) )
    
    varnames = c("tau","rho","alpha","tau0")
    ## MODEL INIT 
    m[[strG]] = jags.model(file=modelFileName,data=modelData[[strG]], n.chains = 3,
                           inits = inits_func);
    
    update(m[[strG]],burn_in)
    
    draw[[strG]] = jags.samples(m[[strG]],steps,thin=thin,variable.names=varnames)
    
  }  
  
  message("Save JAGS results to file:",paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD_s",steps,"-",modelFileName,".RData"))
  save(draw,modelData,m,steps,thin,
       file=paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD",modelFileName,".RData"))
  
  return (list(draw=draw,data=modelData))
}


## Model Init Randomizer
inits_func <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      tau0 = rgamma(1, 1, rate=1/20),
      tau = rgamma(1, 1, rate=1/20),
      #rho = rgamma(1, 1, rate=1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}



#plot(dgamma(1:80,shape=1,rate=1),main="rho")
#plot(dgamma(1:80,shape=100,rate=1),main="tau")
##Sample the Covariance Kernel



modelFileName <- vector()
modelFileName[1] <-model(10,1,0.025)
modelFileName[2] <-model(50,1,0.025)
modelFileName[3] <-model(150,1,0.025)
modelFileName[4] <-model(250,1,0.025)
modelFileName[5] <-model(10,1,0.015)
modelFileName[6] <-model(50,1,0.015)
modelFileName[7] <-model(150,1,0.015)
modelFileName[8] <-model(250,1,0.015)

t = 2
## Prepare Data - 
preyCntRange <- c(0,80) ## Prey Density Range to Include in Model
message(paste(" Loading Dispersion Dat List to Analyse... "))
datDispersion <- loadDispersionData(FALSE,t)  
datHuntLabelledEventsSBMerged <- getLabelledHuntEventsSet()

#inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=100,thin=2,modelFileName[2])
#inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=100,thin=2,modelFileName[2])

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

## Run Across Conditions
vretM = foreach (t_model= modelFileName,.combine=c) %dopar%
{
    
    inferGPModel_MSDVsPreyDensity(burn_in=150,steps=2000,dataSamples=350,thin=2,t_model)
    
}

## Select One to Plot from 
retM <- vretM[1]
modelData <- retM$modelData
draw <- retM$draw



#### Plot Model ####
SE <- function(Xi,Xj, rho,tau) tau^2*exp(-(Xi - Xj) ^ 2 * rho^2)
covC <- function(X, Y, rho,tau) outer(X, Y, SE, rho,tau)

plot_res<- function(ind,drawY,Xn,Yn,colour='red ',qq=0.05,pPch=16,chain=1){
  
  ord=order(Xn)
  Xn=Xn[ord] ##Place Points In order so we can draw the Polygon Bands
  Yn=Yn[ord]
  
  points(Xn,Yn,col=colour, pch=pPch)
  x_predict=seq(preyCntRange[1],preyCntRange[2],1)
  Ef=matrix(NA,ncol=length(x_predict),nrow=ind)
  for(j in 1:ind){
    #i=steps/thin-j+1
    i=NROW(drawY$tau[,,chain]) - j+1
    ##print(i)
    cov_xx_inv=solve(covC(Xn,Xn,drawY$rho[,i,chain],drawY$tau[,i,chain])+diag(drawY$tau0[,i,chain],length(Xn)))
    Ef[j,] <- covC(x_predict, Xn,drawY$rho[,i,chain],drawY$tau[,i,chain]) %*% cov_xx_inv %*% Yn 
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
tauRangeA <- 253
Rho <-1
ind = 10

# 
strSuffix <- modelFileName[2]
#load(file=paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD",strSuffix,".RData"))

retM <- inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=300,thin=2,modelFileName[7])
draw <- retM[[1]]
modelData <- retM[[2]]
plotPDFOutput(modelData,draw,modelFileName[7])

#inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=100,thin=2,modelFileName[2])
plotPDFOutput <- function(modelData,draw,modelFileName)
{
  plot_Chain= 3
  #strPlotName <- paste("plots/stat_HuntEventRateVsPrey_GPEstimate-tauLL",round(mean(draw[["LL"]]$tau)),".pdf",sep="-")
  strPlotName <-  paste(strPlotExportPath,"/stat_MSDVsPreyN",modelData$LF$N,modelFileName,".pdf",sep="")
  pdf(strPlotName,width=8,height=8,title="GP Function of MSD Vs Prey Density") 
  par(mar = c(4.1,4.8,3,1))
  
  plot(modelData$LF$prey,modelData$LF$MSD,col=colourH[1],
       main = NA,
       ylab="Mean squared displacement (mm/ 2 sec)",
       xlab="Initial Prey Density (Rotifers/10ml)",
       cex=1.4,
       cex.axis = 1.7,
       cex.lab = 1.5,
       ylim = c(0,20),##preyCntRange,
       xlim = c(1,65),##preyCntRange,
       #log="x",
       pch=pointTypeScheme$LL,
       #sub=paste("GP tau:",format(mean(draw[["LF"]]$tau),digits=4 ),
       #           "tau0:",format(mean(draw[["LF"]]$tau0),digits=4 ) ,
       #           "rho:",format(mean(draw[["LF"]]$rho),digits=4 ) )  
  )
  
  legend("topright",legend = c(paste("LF #",modelData$LF$N),paste("NF #",modelData$NF$N ),paste("DF #",modelData$DF$N)),
         col=c(colourDataScheme[["LF"]]$Evoked,colourDataScheme[["NF"]]$Evoked,colourDataScheme[["DF"]]$Evoked),
         pch=c(pointTypeScheme$LL,pointTypeScheme$NL,pointTypeScheme$DL ),cex=1.5 )
  
  
  plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=plot_Chain)
  #plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=2)
  #plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=3)
  
  plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL)
  #plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL,chain=2)
  #plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL,chain=3)
  
  plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL)
  #plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL,chain=2)
  #plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL,chain=3)
  dev.off()
}



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

