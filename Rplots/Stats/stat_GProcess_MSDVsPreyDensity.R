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
modelFixedRho <- function(tauShape,tauRate,const_rho) 
{
  ## Show A little sample of what the correlation function looks like
  par(mfrow=c(2,1))
  plot(  ( rgamma(1,shape=tauShape,rate=tauRate)^2) * exp( - (const_rho*seq(-50,50,1) )^2 )  , 
         log="y",ylim=c(0.01,100000),type="l" )
  plot(dgamma(1:80,shape=tauShape,rate=tauRate),main="tau",type='l')
  
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

# Make Model with Parametrizable with Fixed Rho
# \return FileName Where New Model Was Saved
modelVarRho <- function(tauShape,tauRate,rhoShape,rhoRate) 
{
  ## Show A little sample of what the correlation function looks like
  par(mfrow=c(3,1))
  muRho <- mean(dgamma(1:1000,shape=rhoShape,rate=rhoRate))
  plot(  ( rgamma(1,shape=tauShape,rate=tauRate)^2) * exp( - (muRho*seq(-50,50,1))^2) , 
         log="y", ylim=c(0.01,1000000), type="l" ) 
  plot(dgamma(1:100,shape=tauShape,rate=tauRate),main="tau",type='l')
  plot(dgamma(1:100,shape=rhoShape,rate=rhoRate),main=paste("Rho mu:",muRho),type='l')
  
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
    rho ~  dgamma(",rhoShape,",",rhoRate,") # Covariance Width / Links neighbouring Densities Between 
    tau0 ~ dgamma(",tauShape,",",tauRate,") # Covariance Peak Strength
    tau  ~ dgamma(",tauShape,",",tauRate,")  
    
  }")
  
  ## Save to File And Return File Name
  modelFileName <-paste0("model-tauS",tauShape,"R",tauRate,"-rhoS",rhoShape,"R",rhoRate,".tmp")
  
  fileConn=file(modelFileName)
  writeLines(strMdl,fileConn);
  close(fileConn)
  
  return(modelFileName)
}


getUnifPreyCountSample <- function(datDispersionG,dataSamples,dispBreaks)
{
  library(dplyr)
  #Set Factor Indicating Bin Of Prey Density
  datDispersionG$dispRange <- cut(datDispersionG[ ,"PreyCount" ],breaks=dispBreaks,include.lowest = FALSE,right=FALSE)
  
  ##Sample Equally From Each bin so we obtain more uniform MSD data samples across densities
  iSampleLoan <- 0 ##If samples Missing from A bin, then increase # samples for next bin
  datDispersion_Subset <- data.frame()
  ## Sample MSD and PreyDensity From Each Group equally at both ranges (Reduce error in Estimate High Range)
  for (pdbin in levels(datDispersionG$dispRange) )
  {
    #print(pdbin)
    #idxs <- rownames()
    datdatDispersionGBin <- datDispersionG[datDispersionG$dispRange %in% pdbin,]
    ## Take Random Subset + add any samples missed from the previous bin
    datDispersionBin_Subset <- sample_n(datdatDispersionGBin,size=iSampleLoan+round(min(dataSamples,NROW(datdatDispersionGBin))/length(dispBreaks)) ) 
    
    if (NROW(datDispersion_Subset) >0 ){
      datDispersion_Subset <- rbind(datDispersion_Subset,datDispersionBin_Subset)
    }
    else
    {
      datDispersion_Subset <-datDispersionBin_Subset
    }
    
    ##If Samples Where not enough within Bin
    iSampleLoan <- iSampleLoan +  round(dataSamples/length(dispBreaks)) -NROW(datDispersionBin_Subset)
  }
  
  return(datDispersion_Subset)
}

#### Run Model  #### 
# Saves Sampls to RData file and returns Samples and Data - so they can be plotted
# A total of n=dataSamples   Are taken Across preydensities equally - by bining  in widths of 5 
inferGPModel_MSDVsPreyDensity <- function (burn_in=140,steps=10000,dataSamples=120,thin=2,modelFileName,inits_func = inits_func_fixRho)
{
  library(rjags)
  
  modelData <- list()
  m<-list() ##model List
  draw<-list()
  ## Prey Count Bins - Max Set High To Avoid breaking due to  some Noisy measurements
  dispBreaks <- seq(0,preyCntRange[2]+5,5) 
  
  # Truncate Prey Count To Maximum Of Range
  datDispersion[!is.na(datDispersion$PreyCount) & 
                  datDispersion$PreyCount >  preyCntRange[2],"PreyCount"] = preyCntRange[2]
  
 
  
  strRG <- list(LF=c("LE","LL"),NF=c("NE","NL"),DF=c("DE","DL")) #LF=c("LE","LL")
  for (strG in names(strRG) )
  {
    message(strG) 
   
    
    
    ##Sample From Each Prey Density Bin Equally
    datDispersionG <- datDispersion[!is.na(datDispersion$PreyCount) &
                                          !is.na(datDispersion$MSD) &
                                          datDispersion$groupID %in% strRG[[strG]], ]

    ## Correct Prey Count Noise in Empty Condition
    datDispersionG[datDispersionG$groupID == strRG[[strG]][1],"PreyCount" ] = 0
    
    ##Sample Equally From Each bin so we obtain more uniform MSD data samples across densities
    datDispersion_Subset <- getUnifPreyCountSample(datDispersionG,dataSamples,dispBreaks)
    
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
  
  message("Save JAGS results to file:",paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD_N",dataSamples,modelFileName,".RData"))
  
  save(draw,modelData,m,steps,thin,
       file=paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD_N",dataSamples,modelFileName,".RData"))
  
  return (list(draw=draw,data=modelData))
}


## Model Init Randomizer
inits_func_fixRho <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      tau0 = rgamma(1, 50, rate=1),
      tau = rgamma(1, 50, rate=1),
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


## Model Init Randomizer
inits_func_VarRho <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      tau0 = rgamma(1, 50, rate=1),
      tau = rgamma(1, 50, rate=1),
      rho = rgamma(1, 10, rate=1),
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





modelFileName <- vector()
modelFileName[1] <-modelFixedRho(10,1,0.035)
modelFileName[2] <-modelFixedRho(50,1,0.035) #**********
modelFileName[3] <-modelFixedRho(50,1,0.075) #**********
#modelFileName[3] <-modelFixedRho(150,1,0.035)
modelFileName[4] <-modelFixedRho(250,1,0.035)
modelFileName[5] <-modelFixedRho(10,1,0.025)
modelFileName[6] <-modelFixedRho(50,1,0.025) #******** Best One 
modelFileName[7] <-modelFixedRho(150,1,0.025)
modelFileName[8] <-modelFixedRho(250,1,0.025)
modelFileName[9] <-modelFixedRho(10,1,0.015)
modelFileName[10] <-modelFixedRho(50,1,0.015)
modelFileName[11] <-modelFixedRho(150,1,0.015)
modelFileName[12] <-modelFixedRho(250,1,0.015)
modelFileName[13] <-modelFixedRho(10,1,0.035) ## Works Nicely
modelFileName[14] <-modelFixedRho(5,20,0.1) ## Try Weak Correlation - Wider Band
modelFileName[15] <-modelFixedRho(10,20,0.05) ## *****
modelFileName[16] <-modelFixedRho(15,40,0.05) ## Try Weak Correlation - Wider Band
modelFileName[16] <-modelFixedRho(5,10,0.025) ## Try Weak Correlation - Wider Band
modelFileName[17] <-modelFixedRho(430,10,0.025) # * Follow up from 6: Centre Tau more tight around 45 - Covariance Wide
modelFileName[18] <-modelFixedRho(250,6,0.025) # *Follow up from 17:  Make Tau Narrow So confidence Intervals remain wide
modelFileName[19] <-modelFixedRho(850,20,0.025) # *Follow up from 17: Not Tested
modelFileName[20] <-modelFixedRho(1250,30,0.025) # *Follow up from 17:   Not Tested
modelFileName[21] <-modelFixedRho(2250,55,0.025) # *Follow up from 17: - This Showed to be narrow band
modelFileName[22] <-modelFixedRho(20,1/2,0.025) # *** Looks Good  / Follow up from 21: - Make Prior Broader Around Working Region 
modelFileName[23] <-modelFixedRho(21,20,0.05) ## *** Looks Ok (Followup from 15)

##Variable Rho prior
#modelFileName[23] <-modelVarRho(20,1/2,1,10) # ***Follow up from 21: - Make Prior Broader Around Working Region 
modelFileName[24] <-modelVarRho(20,1/2,1/2,10) # ***Follow up from 21: - Make Prior Broader Around Working Region
modelFileName[25] <-modelVarRho(50,1,1/2,10) # ***Follow up from 21: - Make Prior Broader Around Working Region
modelFileName[26] <-modelVarRho(50,1,1,1) # ***Follow up from 21: - Make Prior Broader Around Working Region
modelFileName[27] <-modelVarRho(50,1,10,1) # ***Follow up from 21: - Make Prior Broader Around Working Region

modelFileName[28] <-modelVarRho(15,40,10,1) # XX Straight line *Follow up from 21: 
modelFileName[29] <-modelVarRho(15,40,1,1) # XX Straight line ***Follow up from 21: - Make Prior Broader Around Working Region
modelFileName[30] <-modelVarRho(15,40,1,1/2) # XX Straight line***Follow up from 21: - Make Prior Broader Around Working Region
modelFileName[31] <-modelVarRho(15,40,1,1/4) # XX Straight line ***Follow up from 21: - Make Prior Broader Around Working Region 
#mean(dgamma(-1000:10000,shape=1,rate=1))

#plot(dgamma(1:80,shape=100,rate=1),main="tau")


t = 2
## Prepare Data - 
preyCntRange <- c(0,60) ## Prey Density Range to Include in Model
message(paste(" Loading Dispersion Dat List to Analyse... "))
datDispersion <- loadDispersionData(FALSE,t)  
datHuntLabelledEventsSBMerged <- getLabelledHuntEventsSet()

# 
# #setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-2) #not to overload your computer
# registerDoParallel(cl)
# 
# ## Run Across Conditions
# vretM = foreach (t_model= modelFileName,.combine=c) %dopar%
#   {
#     
#     inferGPModel_MSDVsPreyDensity(burn_in=150,steps=2000,dataSamples=350,thin=2,t_model)
#     
#   }

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
       xlab="Prey Density (Rotifers/10ml)",
       cex=1.4,
       cex.axis = 1.7,
       cex.lab = 1.7,
       ylim = c(0,31),##preyCntRange,
       xlim = c(1,60),##preyCntRange,
       #asp=1,
       #log="y",
       pch=pointTypeScheme$LL,
       #sub=paste("GP tau:",format(mean(draw[["LF"]]$tau),digits=4 ),
       #           "tau0:",format(mean(draw[["LF"]]$tau0),digits=4 ) ,
       #           "rho:",format(mean(draw[["LF"]]$rho),digits=4 ) )  
  )
  
  
  plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=plot_Chain)
  #plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=2)
  #plot_res(ind,draw[["LF"]],modelData$LF$prey,modelData$LF$MSD, colourH[1],0.05,pointTypeScheme$LL,chain=3)
  
  plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL)
  #plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL,chain=2)
  #plot_res(ind,draw[["NF"]],modelData$NF$prey,modelData$NF$MSD,colourH[2],0.05,pointTypeScheme$NL,chain=3)
  
  plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL)
  #plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL,chain=2)
  #plot_res(ind,draw[["DF"]],modelData$DF$prey,modelData$DF$MSD,colourH[3],0.05,pointTypeScheme$DL,chain=3)
  
  legend("topright",legend = c(paste("LF #",modelData$LF$N),paste("NF #",modelData$NF$N ),paste("DF #",modelData$DF$N)),
         col=c(colourDataScheme[["LF"]]$Evoked,colourDataScheme[["NF"]]$Evoked,colourDataScheme[["DF"]]$Evoked),
         pch=c(pointTypeScheme$LL,pointTypeScheme$NL,pointTypeScheme$DL ),cex=1.5 )
  
  
  dev.off()
}

colourH <- c(rgb(0.01,0.7,0.01,0.5),rgb(0.9,0.01,0.01,0.5),rgb(0.01,0.01,0.9,0.5),rgb(0.00,0.00,0.0,1.0))
ind = 10

### RUN MOdel Sequence
for (i in c(6))
{
  retM <- inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=200,thin=2, modelFileName[i] ,inits_func = inits_func_fixRho)
  draw <- retM[[1]]
  modelData <- retM[[2]]
  plotPDFOutput(modelData,draw,modelFileName[i])
}

# 
strSuffix <- "model-tauS50R1-rho0.025.tmp" #modelFileName[2]
load(file=paste0(strDataExportDir,"/jags_GPPreyDensityVsMSD",strSuffix,".RData"))
#load(file=paste0(strDataExportDir,"/jags_GPPreyDensityVsMSDmodel-tauS10R1-rho0.025.tmp.RData"))
plotPDFOutput(modelData,draw,strSuffix)


for (t_model in sample(modelFileName) )
{
  retM <- inferGPModel_MSDVsPreyDensity(burn_in=150,steps=1000,dataSamples=390,thin=2,t_model)
  draw <- retM[[1]]
  modelData <- retM[[2]]
  plotPDFOutput(modelData,draw,t_model)
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




### TEST BSAR model 
# bsamGP: An R Package for Bayesian Spectral Analysis Models Using Gaussian Process Priors Seongil

library(bsamGP)

##Sample Equally From Each bin so we obtain more uniform MSD data samples across densities
datDispersion_Subset <- NA
strRG <- list(LF=c("LE","LL"),NF=c("NE","NL"),DF=c("DE","DL")) #LF=c("LE","LL")
strG <- "LF"
  message(strG) 
  ##Sample From Each Prey Density Bin Equally
  datDispersionG <- datDispersion[!is.na(datDispersion$PreyCount) &
                                    !is.na(datDispersion$MSD) &
                                    datDispersion$groupID %in% strRG[[strG]], ]
  
  ## Correct Prey Count Noise in Empty Condition
  datDispersionG[datDispersionG$groupID == strRG[[strG]][1],"PreyCount" ] = 0
  ##Sample Equally From Each bin so we obtain more uniform MSD data samples across densities
  datDispersion_Subset <- getUnifPreyCountSample(datDispersionG,dataSamples,dispBreaks)



## Synthetic Data - Example
n <- NROW(datDispersion_Subset)

f <- function(x) {
   5-10*x+8*exp(-100*(x-0.3)^2)-
     8 * exp(- 100 * (x - 0.7) ^ 2)
}

w <- rnorm(n, sd = 0.5)
x <- runif(n)
y <- 2 * w + f(x) + rnorm(n, sd = 1)


y <- datDispersion_Subset$MSD
x <- datDispersion_Subset$PreyCount

plot(x,y)

nbasis <- 10
prior <- list(beta_m0 = numeric(20), beta_v0 = diag(100, 2), w0 = 20,
              tau2_m0 = 10, tau2_v0 = 10, sigma2_m0 = 5, sigma2_v0 = 10)

mcmc <- list(nblow = 10000, nskip = 10, smcmc = 5000, ndisp = 5000)
fout <- bsar(y ~ w + fs(x), nbasis = nbasis, mcmc = mcmc, prior = prior,
             shape = "Free", marginal.likelihood = TRUE, spm.adequacy = TRUE)

library("coda")
post <- as.mcmc(data.frame(beta = fout$mcmc.draws$beta,                 sigma2 = fout$mcmc.draws$sigma2))
plot(post)
plot(fout)

##Now Predict 
fit <- fitted(fout, HPD = FALSE)
plot(fit, ggplot2 = FALSE)
