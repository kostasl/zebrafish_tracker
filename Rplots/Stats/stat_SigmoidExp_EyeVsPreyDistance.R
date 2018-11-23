##  24-10-2018 - Estimates Vergence OnSet Distance - 
### Fitting a sigmoid and Exp to the eye Vergence Data of Retracked Hunt Events (The same ones that we used to show the underhooting)
## Assumes A slow and Fast Muscle Action
### Model fits Eye Vergence / Detecting Onset And Rate Of Converge In the Near Prey Region/After Vergence 
## Produces a plot comparing onset distance (τ) , with no striking distance shown actually LL seems to be a proader density
### Note : 20 points Padding is added before the furthest point, making Phi Vergence angle 0, such that lowest V angle Of Sigmoid sits low.
##


##Model Each Hunt Event Individually / And obtain Group Statistic and Regresion of eye vergence vs Distance To Prey
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")
library(rjags)
library(runjags)

##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.3),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
pchL <- c(16,2,4)


#
#These RC params Work Well to Smooth LF And NF
burn_in=1000;
steps=12000;
thin=3;
nchains <-4

dataFrac <- 1.0 ##Fraction Of Hunt Episodes to Include in DataSet
sampleFraction  <- 0.65 ##Fraction of Points to Use from Each Hunt Episode's data
fitseqNo <- 5
npad <- 1

##THe Growth Model : Carlin and Gelfand (1991) present a nonconjugate Bayesian analysis of the following data set from Ratkowsky (1983):
modelGCSigmoidInd  <- "model
{
  
  for( i in 1 : N ) {
   phi_hat[ hidx[i],i] <-  phi_0[ hidx[i] ] +   (phi_max[hidx[i]] - phi_0[ hidx[i] ])/( 1 + exp( -lambda[ hidx[i] ]*( ( tau[ hidx[i] ] - distP[i]    ) ) ) )

   ###OUT Set Region Of Exp Growth Model ##
   # s[hidx[i],i] <- step( distP[i] - u1[hidx[i]] )*step( tau[ hidx[i] ] -distP[i] )     # step( phi_max[hidx[i]] - phi_0[hidx[i]] ) #step(u0[ hidx[i] ] - distP[i]  )  
   
   ## Define Exp Growth Model 
   phi_exp[ hidx[i],i] <- alpha[hidx[i]]*exp( gamma[ hidx[i] ]* ( tau[ hidx[i] ] -  distP[i]))
   
   ### Conditionally Include the exp Model
   phi[i] ~ dnorm(phi_exp[ hidx[i],i]  + phi_hat[ hidx[i],i]  , var_inv[hidx[i]] ) #s[hidx[i],i]+1 
   

  }
  
  
  ## Priors
  limDist <- max(distMax)

  
  for(i in 1:max(hidx) ) { 
    phi_max[i] ~ dnorm(65,1e-3)T(0,150) ##I(0,100) # Max Eye Vergence Angle
    phi_0[i] ~ dnorm(0.01, 1e-3)T(0,60)  # Idle Eye Position
    lambda[i] ~ dgamma(1, 1) #dnorm(100.0, 1e-3)T(0,) # RiseRate of Eye Vs Prey Distance Sigmoid
    gamma[i] ~ dgamma(1, 1) #dnorm(0.5, 1e-3)I(0,)  # RiseRate of Eye Vs Prey Distance After Sig Rise dunif(0.5, 0.000001)
    alpha[i] ~ dunif(1,3)
    tau[i] ~ dnorm(distMax[i], 1e-2) ##inflexion point, sample from where furthest point of Hunt event is found
    var_inv[i] ~ dgamma(0.001, 0.001) ##Draw   ##Precision
    
    sigma[i] <- 1 / sqrt(var_inv[ i])    

  }
  
  
  
  
  
}"
  
  ## Plot the Eye Vs Distance data points and the regression variations ##
  plotEyeGCFit <- function(pp,strGroup,dataSubset,drawS)
  {
    vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
    
    vX  <- seq(0,5,by=0.01)
    vPP <- which (dataSubset$hidx == pp)
    
    etau    <- (tail(drawS$tau[pp,,],n=100))
    ephimax <- (tail(drawS$phi_max[pp,,],n=100))
    ephi0   <- (tail(drawS$phi_0[pp,,],n=100))
    elambda <- (tail(drawS$lambda[pp,,],n=100))
    egamma  <- (tail(drawS$gamma[pp,,],n=100))
    ealpha  <- (tail(drawS$alpha[pp,,],n=100))
    
    plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(0,85),
         main=paste(strGroup,pp," (",vRegIdx[pp],")"), 
         bg=colourP[2],col=colourP[1],cex=0.5)
    ## Draw The 100 Variotons before the fit converged      
    for (k in 1:NROW(etau) )
    {
      vY  <-  ealpha[k]*exp(egamma[k]*(etau[k]-vX) )+ ephi0[k]   +  (ephimax[k] -ephi0[k]  )/(1+exp( -(elambda[k]   *(etau[k] -vX )   ) ) ) 
      
      #vY_l  <- quantile(drawS$phi_0[pp,,])[1]   - ( quantile (drawS$lambda[pp])[1] )*((( quantile(drawS$gamma[pp,,])[1] )^( quantile(drawS$u0[pp])[1] - (vX) ) ) ) #
      #vY_u  <- quantile(drawS$phi_0[pp,,])[5]   - (quantile (drawS$lambda[pp,,])[5])*((( quantile(drawS$gamma[pp,,])[5] )^( quantile(drawS$u0[pp,,])[5] - (vX) ) ) ) #
      #      #points(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(-85,85),main=paste("L",pp), bg=colourP[2],col=colourP[1],cex=0.5)
      lines( vX ,vY,type="l",col=colourR[3],lwd=1)
      #        lines( vX ,vY_u,type="l",col=colourR[4],lwd=1)
    }

  }
  
  
  
  ## Plot average regressed function ##
  ## plot( exp(0.1*(-vx+80))+  10 + (90-10)/(1+exp(-100*(60-vx) ))   ,ylim=c(0,400))
  plotGCSig <- function (drawS,dataSubset,n=NA,groupID){
    
    bPlotIndividualEvents <- FALSE
    ## compute 2D kernel density, see MASS book, pp. 130-131
    max_x <- 7
    nlevels <- 12
    
    if (is.na(n))
      n <- NROW(unique(dataSubset$hidx))
    
    vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
    
    vsampleP <- sample(unique(dataSubset$hidx),n)
    vsub <- which (dataSubset$hidx %in% vsampleP)
    
    z <- kde2d(dataSubset$distP, dataSubset$phi, n=80)
    
    # X11()
    
    plot(dataSubset$distP[vsub],dataSubset$phi[vsub],pch=21,xlim=c(0,max_x),ylim=c(0,80),
         main=paste("Model Fit : Eye Vergence Vs Distance Data ",strGroupID[groupID]),
                    ylab=expression(paste("Eye Vergence ",Phi," (degrees)") ),
                    xlab=expression(paste("Distance from Prey (mm)") ),
                     bg=colourP[groupID],col="#FFFFFFAA",cex=0.5)
    #points(dataSubset$distToPrey[vsub],dataSubset$vAngle[vsub],pch=21,xlim=c(0,5),ylim=c(0,80),main="LL", bg=colourP[4],col=colourP[1],cex=0.5)
    contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
    ## Plot The Mean Curve of the selected subset of curves
    vX  <- seq(0,max_x,by=0.01)##max(drawS$u0[vsampleP])
    ##  (phi_max[hidx[i]] - phi_0[hidx[i]] )/(1-exp(-lambda[ hidx[i] ]*(u0[ hidx[i] ]  - distP[i] ) ))
    etau    <- mean(tail(drawS$tau[vsampleP],n=100))
    ephimax <- mean(tail(drawS$phi_max[vsampleP],n=100))
    ephi0   <- mean(tail(drawS$phi_0[vsampleP],n=100))
    elambda <- mean(tail(drawS$lambda[vsampleP],n=100))
    egamma  <- mean(tail(drawS$gamma[vsampleP],n=100))
    ealpha  <- mean(tail(drawS$alpha[vsampleP],n=100))
    
    vY  <-    ealpha*exp(egamma*(etau-vX) ) + ephi0 + (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    
    etau <- quantile((drawS$tau[vsampleP]))[2]
    vY_l  <-  ealpha*exp(egamma*(etau-vX) )+  ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    
    etau <- quantile((drawS$tau[vsampleP]))[4]
    vY_u  <-  ealpha*exp(egamma*(etau-vX) )+  ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    
    
    #vY_u <-  quantile(drawS$phi_0[vsampleP])[4]-(quantile(drawS$lambda[vsampleP])[4])*((quantile(drawS$gamma[vsampleP])[4]^( quantile(drawS$u0[vsampleP])[4] - (vX) ) ) )
    #vY_l <-  quantile(drawS$phi_0[vsampleP])[2]-(quantile(drawS$lambda[vsampleP])[2])*((quantile(drawS$gamma[vsampleP])[2]^( quantile(drawS$u0[vsampleP])[2] - (vX) ) ) )
    lines( vX ,vY,xlim=c(0,max_x),ylim=c(0,80),type="l",col="black",lwd=3)
    lines( vX ,vY_u,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    lines( vX ,vY_l,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    
    dev.off()
    
    if (!bPlotIndividualEvents)
    return(NA)
    ##plot individual Curve Fits
    for (pp in vsampleP) ##Go through Each hidx 
    {
      #phi_0[hidx[i]] - lambda[ hidx[i] ] * pow(gamma[hidx[i]],distMax[i] - distP[i] )   
      #      #X11()

      pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",vRegIdx[pp],".pdf",sep="")) 

      plotEyeGCFit(pp,strGroupID[groupID],dataSubset,drawS)
    
      dev.off()
    
  } ##For Each Sampled Hunt Event 
    
}##END oF Function 
  
plotConvergenceDiagnostics <- function(strGroupID,drawS,dataS)
{
  
  vRegIdx <- unique(dataS$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
  N <- NROW(drawS$tau[,1,1])
  for (idxH in 1:N)
  {
    
    ## Plot multipage for param convergence ##
    pdf(onefile=TRUE,file= paste(strPlotExportPath,"/stat/diag/stat_SigExpFit_",strGroupID,"_",vRegIdx[idxH],".pdf",sep="")) 
    ## plot the regression lines and the data
    plotEyeGCFit(idxH,strGroupID,dataS,drawS) 
    
    plot(drawS$tau[idxH,,1],type='l',ylim=c(0,4),main=paste("tau",idxH," (",vRegIdx[idxH],")") )
    lines(drawS$tau[idxH,,2],type='l',col="red")
    lines(drawS$tau[idxH,,3],type='l',col="blue")
    #dev.off()
    
    ##gelmal rubin diag ##
    print(paste(idxH," :--" )  )
    chains_tau <- mcmc(drawS$tau[idxH,,],thin=thin)
    lmcmc_tau <- mcmc.list(chains_tau[,1],chains_tau[,2],chains_tau[,3])
    gdiag_psrf <- gelman.diag(lmcmc_tau,autoburnin=TRUE )$psrf
    #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_gelman_SigExpFit_gamma",idxH,".pdf",sep="")) 
    gelman.plot( lmcmc_tau,autoburnin=TRUE,max.bins=100, ylim=c(0.99,1.5),
                 main=paste("tau psrf:", round(gdiag_psrf[1]*100)/100 ) )
    
    
    
    #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_SigExpFit_gamma",idxH,".pdf",sep="")) 
    plot(drawS$gamma[idxH,,1],type='l',ylim=c(0,4),main=paste("V rise rate gamma ",idxH," (",vRegIdx[idxH],")") )
    lines(drawS$gamma[idxH,,2],type='l',col="red")
    lines(drawS$gamma[idxH,,3],type='l',col="blue")
    #dev.off()
    
    chains_gamma <- mcmc(drawS$gamma[idxH,,],thin=thin)
    lmcmc_gamma <- mcmc.list(chains_gamma[,1],chains_gamma[,2],chains_gamma[,3])
    gdiag_psrf <- gelman.diag(lmcmc_gamma,autoburnin=TRUE )$psrf
    #pdf(file= paste(strPlotExportPath,"/stat/diag/stat_gelman_SigExpFit_gamma",idxH,".pdf",sep="")) 
    gelman.plot( lmcmc_gamma,autoburnin=TRUE,max.bins=100, ylim=c(0.99,1.5),
                 main=paste("gamma psrf:", round(gdiag_psrf[1]*100)/100 ) )
    
    
    dev.off()
  }
}
  
  
  
  
  ####Select Subset Of Data To Analyse
  
  strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
  message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
  datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
  
  lEyeMotionDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData.rds",sep="") ) #Processed Registry on which we add )
  
  datEyeVsPreyCombinedAll <-  data.frame( do.call(rbind,lEyeMotionDat ) )
  
  datEyeVsPreyCombinedAll <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$doesCaptureStrike == 1,] # Select The Ones With a strike Capture
  
  strGroupID <- levels(datTrackedEventsRegister$groupID)
  
  
  ##Add The Empty Test Conditions
  #strProcDataFileName <-paste("setn14-D5-18-HuntEvents-Merged",sep="") ##To Which To Save After Loading
  #datHuntLabelledEventsKL <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))
  #datHuntStatE <- makeHuntStat(datHuntLabelledEventsKL)
  #datHuntLabelledEventsKLEmpty <- datHuntLabelledEventsKL[datHuntLabelledEventsKL$groupID %in% c("DE","LE","NE"),]
  lRegIdx <- list()
  ldatsubSet <-list()
  
  ## Get Event Counts Within Range ##
  ldatREyePoints <- list()
  ldatLEyePoints <- list()
  ldatVEyePoints <- list()
  lnDat          <- list()
  lnMaxDistanceToPrey <- list()
  
  ##Do all this processing to add a sequence index To The hunt Event + make vergence angle INdex 
  for (g in strGroupID) {
    lRegIdx[[g]] <- unique(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == g),"RegistarIdx"])
    ldatLEyePoints[[g]] <- list()
    
    ##Take Each Hunting Event of Group And Select a Subset of the records. add to dataStruct List which will be fed to Model
    for (h in 1:NROW(lRegIdx[[g]]) )
    {
      ldatsubSet[[g]] <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == g) &
                                                   datEyeVsPreyCombinedAll$RegistarIdx %in% lRegIdx[[g]][h] ,]  ##Select THe Ones With Capture
      
      ldatsubSet[[g]] <- ldatsubSet[[g]][sample(NROW(ldatsubSet[[g]]),sampleFraction*NROW(ldatsubSet[[g]] ) ) ,] ##Sample Points 
      
      ldatLEyePoints[[g]][[h]] <- cbind(ldatsubSet[[g]]$LEyeAngle,
                                        as.numeric(ldatsubSet[[g]]$DistToPrey),
                                        as.numeric(ldatsubSet[[g]]$DistToPreyInit ),
                                        ldatsubSet[[g]]$RegistarIdx,
                                        h)
      
      ldatREyePoints[[g]][[h]] <- cbind(ldatsubSet[[g]]$REyeAngle,
                                        as.numeric(ldatsubSet[[g]]$DistToPrey),
                                        as.numeric(ldatsubSet[[g]]$DistToPreyInit ),
                                        ldatsubSet[[g]]$RegistarIdx,
                                        h)
      
      ##Make Data      
      ldatVEyePoints[[g]][[h]] <- cbind(
        vAngle=ldatsubSet[[g]]$LEyeAngle-ldatsubSet[[g]]$REyeAngle,
        distToPrey=as.numeric(ldatsubSet[[g]]$DistToPrey),
        initDistToPrey=as.numeric(ldatsubSet[[g]]$DistToPreyInit ),
        RegistarIdx=ldatsubSet[[g]]$RegistarIdx,
        seqIdx=h)
      
      
      ##Augment with Idle phi entries for this hunt Event- to go up to 6 mm - 0 
      missingRegion <- 6 - head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1) 
     
      missingRegion <- 0 ##Do Not Pad The data
            
      if (missingRegion > 0)
      {
        ### Max Angle When Info Is Missing is 10
        datpadding <- cbind(vAngle=rep( max(0.1,min( c(10,ldatVEyePoints[[g]][[h]][,"vAngle"] ) ) )  ,npad),
                            distToPrey = seq(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),6,length=npad),
                            initDistToPrey = rep(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),npad),
                            RegistarIdx = rep(head(as.numeric(ldatsubSet[[g]]$RegistarIdx ),n=1),npad),
                            seqIdx = h)
        
         ldatVEyePoints[[g]][[h]] <- rbind(ldatVEyePoints[[g]][[h]],datpadding) ##Add The Zero Phi Data
      }
      
      
      #ldatLEyePoints[[g]][[h]] <- ldatLEyePoints[[g]][[h]][!is.na(ldatLEyePoints[[g]][[h]][,2]),]
      #ldatREyePoints[[g]][[h]] <- ldatREyePoints[[g]][[h]][!is.na(ldatREyePoints[[g]][[h]][,2]),]
      #ldatVEyePoints[[g]][[h]] <- ldatVEyePoints[[g]][[h]][!is.na(ldatVEyePoints[[g]][[h]][,2]),]
      
      lnDat[[g]][[h]] <- NROW(ldatLEyePoints[[g]][[h]]) ##Not Used Anymore
      lnMaxDistanceToPrey[[g]][[h]] <- as.numeric(head(ldatsubSet[[g]]$DistToPreyInit,1)  ) ##Hold Unique Value Of Max Distance To Prey
    } ##For Each Hunt Event

  } ##For Each Group[]
  datVEyePointsLL <- data.frame( do.call(rbind,ldatVEyePoints[["LL"]] ) ) 
  datVEyePointsNL <- data.frame( do.call(rbind,ldatVEyePoints[["NL"]] ) ) 
  datVEyePointsDL <- data.frame( do.call(rbind,ldatVEyePoints[["DL"]] ) ) 
  
  
  
  
  ##Larva Event Counts Slice
  nDatLL <- NROW(datVEyePointsLL)
  nDatNL <- NROW(datVEyePointsNL)
  nDatDL <- NROW(datVEyePointsDL)
  
  ##Test limit data
  ## Subset Dat For Speed
  vsubIdx <- sample(NROW(lRegIdx[["LL"]]),NROW(lRegIdx[["LL"]])*dataFrac)
  datVEyePointsLL_Sub <- datVEyePointsLL[datVEyePointsLL$seqIdx %in% vsubIdx ,] #
  dataLL=list(phi=datVEyePointsLL_Sub$vAngle,
              distP=datVEyePointsLL_Sub$distToPrey ,
              N=NROW(datVEyePointsLL_Sub),
              distMax=lnMaxDistanceToPrey[["LL"]], #Put All distances in So We can Ref By Index #datVEyePointsLL_Sub$initDistToPrey,
              hidx=datVEyePointsLL_Sub$seqIdx, ##Define an idx Key To link to the 
              RegistrarIdx=datVEyePointsLL_Sub$RegistarIdx);
  
  
  ##Test limit data
  ## Subset Dat For Speed
  
  vsubIdx <-sample(NROW(lRegIdx[["NL"]]),NROW(lRegIdx[["NL"]])*dataFrac)
  datVEyePointsNL_Sub <- datVEyePointsNL[datVEyePointsNL$seqIdx %in% vsubIdx,] 
  dataNL=list(phi=datVEyePointsNL_Sub$vAngle,
              distP=datVEyePointsNL_Sub$distToPrey ,
              N=NROW(datVEyePointsNL_Sub),
              distMax=lnMaxDistanceToPrey[["NL"]],#datVEyePointsNL_Sub$initDistToPrey,
              hidx=datVEyePointsNL_Sub$seqIdx,
              RegistrarIdx=datVEyePointsNL_Sub$RegistarIdx);
  
  ##Test limit data
  ## Subset Dat For Speed
  vsubIdx <-sample(NROW(lRegIdx[["DL"]]),NROW(lRegIdx[["DL"]])*dataFrac)
  datVEyePointsDL_Sub <- datVEyePointsDL[datVEyePointsDL$seqIdx %in% vsubIdx ,] 
  dataDL=list(phi=datVEyePointsDL_Sub$vAngle,
              distP=datVEyePointsDL_Sub$distToPrey ,
              N=NROW(datVEyePointsDL_Sub),
              distMax=lnMaxDistanceToPrey[["DL"]],
              hidx=datVEyePointsDL_Sub$seqIdx,
              RegistrarIdx=datVEyePointsDL_Sub$RegistarIdx);
  
  
  
  varnames=c("phi_0","phi_max","lambda","gamma","sigma","alpha","tau") #"gamma"
  
  
  
  fileConn=file("modelSig.tmp")
  #writeLines(modelGPV1,fileConn);
  writeLines(modelGCSigmoidInd,fileConn);
  close(fileConn)
  
  mLL=jags.model(file="modelSig.tmp",n.chains=nchains,data=dataLL);
  update(mLL,burn_in);
  #drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
   resultsLL <- run.jags(mLL,method = "parallel",monitor = varnames,n.chains = nchains,data=dataLL,thin = thin, sample = steps )
  drawLL
  #sampLL <- coda.samples(mLL,                      variable.names=varnames,                      n.iter=steps, progress.bar="none")
  
  #X11()

  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_LL_C.pdf",sep="")) 

  plotGCSig(drawLL,dataLL,n=NA,groupID=2)
  #dev.off()
  #plotExpRes(drawLL,dataLL)
  
  
  
  
  #dev.off()
  ########################
  ## NL ###
  mNL=jags.model(file="modelSig.tmp",n.chains=nchains,data=dataNL);
  update(mNL,burn_in)
  drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_NL_C.pdf",sep="")) 
  plotGCSig(drawNL,dataNL,n=NA,groupID=3)
  #dev.off()
  ############
  ### DL ###
  mDL=jags.model(file="modelSig.tmp",n.chains=nchains,data=dataDL);
  update(mDL,burn_in)
  drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_DL_C.pdf",sep="")) 
  plotGCSig(drawDL,dataDL,n=NA,groupID=1)
  #dev.off()
  
  
  ## SHOW Hunt OnSet Densities ###
  ## Plot Vergence Angle Distributions/Densities
  dLLphi<-density(drawLL$tau)
  dDLphi<-density(drawDL$tau)
  dNLphi<-density(drawNL$tau)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_SigmoidFit_CompareOnset_tau2.pdf",sep="")) 
  #X11()
  plot(dLLphi,col=colourH[2],type="l",lwd=2,ylim=c(0,1.0),main="Vergence Onset Vs Distance",xlab=expression(paste(" ",tau, " (mm)" ) ) )
  lines(dNLphi,col=colourH[3],lwd=2)
  lines(dDLphi,col=colourH[1],lwd=2)
  #legend("topright",legend =  strGroupID,fill=colourH)
  legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lnDat[["DL"]]),NROW(lnDat[["LL"]]) ,NROW(lnDat[["NL"]]) ) )
         ,fill=colourH,lty=c(1,2,3))
  dev.off()
  
  save(dataLL,dataDL,dataNL,drawLL,drawDL,drawNL,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_",fitseqNo,".RData",sep=""))      
  
  ##Save All  
  save.image(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit",fitseqNo,".RData",sep="") )
       
####################
  
  ############ Checking Convergence ####
  ## Compare convergence between the 3 chains for each trace 
  
  #### CalcInformation ##
  load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_",fitseqNo,".RData",sep=""))
  
  idxH <- 1
  
  
  ## The gelman.diag gives you the scale reduction factors for each parameter.
  ## A factor of 1 means that between variance and within chain variance are equal, larger 
  # values mean that there is still a notable difference between chains. 
  
  dataS <- dataDL
  drawS <- drawDL
  strGroupID <- "DL"
  plotConvergenceDiagnostics(strGroupID,drawS, dataS)

  dataS <- dataNL
  drawS <- drawNL
  strGroupID <- "NL"
  plotConvergenceDiagnostics(strGroupID,drawS, dataS)
  
  dataS <- dataLL
  drawS <- drawLL
  strGroupID <- "LL"
  plotConvergenceDiagnostics(strGroupID,drawS, dataS)
  
  