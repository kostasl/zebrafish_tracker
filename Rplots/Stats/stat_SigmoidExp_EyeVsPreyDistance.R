##  24-10-2018 - Estimates Vergence OnSet Distance - 
### Fitting a sigmoid to the eye Vergence Data of Retracked Hunt Events (The same ones that we used to show the underhooting)
### Model To Detect Onset Of Vergence And Compare HuntOnset To Distance From Prey Among Groups
### 20 points Padding is added before the furthest point, making Phi Vergence angle 0, such that lowest V angle Of Sigmoid sits low.
## Produces a plot comparing onset distance (τ) , with no striking distance shown actually LL seems to be a proader density
##


##Model Each Hunt Event Individually / And obtain Group Statistic and Regresion of eye vergence vs Distance To Prey
source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

##THe Growth Model : Carlin and Gelfand (1991) present a nonconjugate Bayesian analysis of the following data set from Ratkowsky (1983):
modelGCSigmoidInd  <- "model
{

  for( i in 1 : N ) {  
   phi_hat_sigmoid[hidx[i],i] <- phi_0[hidx[i]]  
                                +(phi_max[hidx[i]] - phi_0[ hidx[i] ])/( 1 + exp( -gamma[hidx[i]]*( ( tau[ hidx[i] ] - distP[i]) ) ) )
                                +exp(lambda[hidx[i]]*( tau[ hidx[i] ] - distP[i]))
   phi[i] ~ dnorm( phi_hat[ hidx[i],i], sigma_inv[hidx[i]] ) 
  }
  
  
  ## Priors
  limDist <- max(distMax)

  for(i in 1:max(hidx) ) {

    phi_max[i] ~ dnorm(65,1e-3) ##I(0,100) # Max Eye Vergence Angle
    phi_0[i] ~ dnorm(1.0, 1e-3)I(0,phi_max[i]) # Idle Eye Position
  
    gamma[i] ~ dnorm(100.0,10 )I(10,) #dgamma(1, 1) # RiseRate of Eye Vs Prey Distance
    lambda[i] ~ dgamma(1, 1)
    tau[i] ~ dnorm(distMax[i], 1e-2) ##inflexion point, sample from where furthest point of Hunt event is found

  # Sigma On Eye Angle when  In Or Out of hunt region 
  
    sigma_inv[i] ~ dgamma(0.001, 0.001) ##Draw 
  
  }
  
  
}"
  
  
  plotGCSig <- function (drawS,dataSubset,n=NA,groupID){
    
    ## compute 2D kernel density, see MASS book, pp. 130-131
    max_x <- 7
    nlevels <- 12
    
    if (is.na(n))
      n <- NROW(unique(dataSubset$hidx))
    
    vsampleP <- sample(unique(dataSubset$hidx),n)
    vsub <- which (dataSubset$hidx %in% vsampleP)
    
    z <- kde2d(dataSubset$distP, dataSubset$phi, n=80)
    
    # X11()
    plot(dataSubset$distP[vsub],dataSubset$phi[vsub],pch=21,xlim=c(0,max_x),ylim=c(0,80),main=paste("Sigmoid fit to eye vergence onset",
                                                                                                    strGroupID[groupID] ), bg=colourP[groupID],col="#FFFFFFAA",cex=0.5)
    #points(dataSubset$distToPrey[vsub],dataSubset$vAngle[vsub],pch=21,xlim=c(0,5),ylim=c(0,80),main="LL", bg=colourP[4],col=colourP[1],cex=0.5)
    contour(z, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
    ## Plot The Mean Curve of the selected subset of curves
    vX  <- seq(0,max_x,by=0.01)##max(drawS$u0[vsampleP])
    ##  (phi_max[hidx[i]] - phi_0[hidx[i]] )/(1-exp(-lambda[ hidx[i] ]*(u0[ hidx[i] ]  - distP[i] ) ))
    etau    <- mean(tail(drawS$tau[vsampleP],n=100))
    ephimax <- mean(tail(drawS$phi_max[vsampleP],n=100))
    ephi0   <- mean(tail(drawS$phi_0[vsampleP],n=100))
    elambda <- mean(tail(drawS$lambda[vsampleP],n=100))
    egamma <- mean(tail(drawS$gamma[vsampleP],n=100))
    
    vY  <-    ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) +  
    
    etau <- quantile((drawS$tau[vsampleP]))[2]
    vY_l  <-   ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    
    etau <- quantile((drawS$tau[vsampleP]))[4]
    vY_u  <-   ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
    
    
    #vY_u <-  quantile(drawS$phi_0[vsampleP])[4]-(quantile(drawS$lambda[vsampleP])[4])*((quantile(drawS$gamma[vsampleP])[4]^( quantile(drawS$u0[vsampleP])[4] - (vX) ) ) )
    #vY_l <-  quantile(drawS$phi_0[vsampleP])[2]-(quantile(drawS$lambda[vsampleP])[2])*((quantile(drawS$gamma[vsampleP])[2]^( quantile(drawS$u0[vsampleP])[2] - (vX) ) ) )
    lines( vX ,vY,xlim=c(0,max_x),ylim=c(0,80),type="l",col="black",lwd=3)
    lines( vX ,vY_u,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    lines( vX ,vY_l,xlim=c(0,max_x),ylim=c(0,80),type="l",col="red",lwd=0.5)
    
    dev.off()
    
    ##plot individual Curve Fits
    for (pp in vsampleP)
    {
      #phi_0[hidx[i]] - lambda[ hidx[i] ] * pow(gamma[hidx[i]],distMax[i] - distP[i] )   
      vX  <- seq(0,5,by=0.01)
      etau <- drawS$tau[pp]
      ephimax <- drawS$phi_max[pp] 
      ephi0 <- drawS$phi_0[pp]
      elambda <- drawS$lambda[pp]
      vY  <-    ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
      
      #elambda <- quantile(tail(drawS$lambda[pp,,],n=300))[2]
      etau <- quantile(tail(drawS$tau[pp,,],n=300))[2]
      #ephimax <- quantile(tail(drawS$phi_max[pp,,],n=300))[2]
      #ephi0 <- quantile(tail(drawS$phi_0[pp,,],n=300))[2]
      #vY    <- (drawS$phi_0[pp] ) - ( (drawS$lambda[pp]))*(((drawS$gamma[pp])^( drawS$u0[pp] - (vX) ) ) ) #
      vY_l  <-   ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
      
      #elambda <- quantile(tail(drawS$lambda[pp,,],n=300))[4]
      etau <- quantile(tail(drawS$tau[pp,,],n=300))[4]
      #ephimax <- quantile(tail(drawS$phi_max[pp,,],n=300))[4]
      #ephi0 <- quantile(tail(drawS$phi_0[pp,,],n=300))[4]
      
      vY_u  <-   ephi0   +  (ephimax -ephi0  )/(1+exp( -(elambda )  *(etau -(vX)   ) ) ) 
      
      #vY_l  <- quantile(drawS$phi_0[pp,,])[1]   - ( quantile (drawS$lambda[pp])[1] )*((( quantile(drawS$gamma[pp,,])[1] )^( quantile(drawS$u0[pp])[1] - (vX) ) ) ) #
      #vY_u  <- quantile(drawS$phi_0[pp,,])[5]   - (quantile (drawS$lambda[pp,,])[5])*((( quantile(drawS$gamma[pp,,])[5] )^( quantile(drawS$u0[pp,,])[5] - (vX) ) ) ) #
      
      vPP <- which (dataSubset$hidx == pp)
      
      pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",pp,".pdf",sep="")) 
      #      #X11()
      plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(0,85),main=paste(strGroupID[groupID],pp), bg=colourP[2],col=colourP[1],cex=0.5)
      #      #points(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,5),ylim=c(-85,85),main=paste("L",pp), bg=colourP[2],col=colourP[1],cex=0.5)
      lines( vX ,vY,type="l",col=colourR[3],lwd=2)
      lines( vX ,vY_l,type="l",col=colourR[4],lwd=1)
      lines( vX ,vY_u,type="l",col=colourR[4],lwd=1)
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
  sampleFraction  <- 0.25
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
      npad <- 20

            
      if (missingRegion > 0)
      {
        datpadding <- cbind(vAngle=rep(0,npad),
                            distToPrey = seq(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),6,length=npad),
                            initDistToPrey = rep(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),npad),
                            RegistarIdx = rep(head(as.numeric(ldatsubSet[[g]]$RegistarIdx ),n=1),npad),
                            seqIdx = h)
        
 #        ldatVEyePoints[[g]][[h]] <- rbind(ldatVEyePoints[[g]][[h]],datpadding) ##Add The Zero Phi Data
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
  
  
  ##For the 3 Groups 
  colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
  colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points DL,LL,NL
  colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
  pchL <- c(16,2,4)
  #
  #These RC params Work Well to Smooth LF And NF
  burn_in=100;
  steps=10000;
  thin=10;
  
  dataFrac <- 1.0 ##Fraction Of Hunt Episodes to Include in DataSet
  
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
              hidx=datVEyePointsLL_Sub$seqIdx );
  
  
  ##Test limit data
  ## Subset Dat For Speed
  
  vsubIdx <-sample(NROW(lRegIdx[["NL"]]),NROW(lRegIdx[["NL"]])*dataFrac)
  datVEyePointsNL_Sub <- datVEyePointsNL[datVEyePointsNL$seqIdx %in% vsubIdx,] 
  dataNL=list(phi=datVEyePointsNL_Sub$vAngle,
              distP=datVEyePointsNL_Sub$distToPrey ,
              N=NROW(datVEyePointsNL_Sub),
              distMax=lnMaxDistanceToPrey[["NL"]],#datVEyePointsNL_Sub$initDistToPrey,
              hidx=datVEyePointsNL_Sub$seqIdx );
  
  ##Test limit data
  ## Subset Dat For Speed
  vsubIdx <-sample(NROW(lRegIdx[["DL"]]),NROW(lRegIdx[["DL"]])*dataFrac)
  datVEyePointsDL_Sub <- datVEyePointsDL[datVEyePointsDL$seqIdx %in% vsubIdx ,] 
  dataDL=list(phi=datVEyePointsDL_Sub$vAngle,
              distP=datVEyePointsDL_Sub$distToPrey ,
              N=NROW(datVEyePointsDL_Sub),
              distMax=lnMaxDistanceToPrey[["DL"]],
              hidx=datVEyePointsDL_Sub$seqIdx );
  
  
  
  varnames=c("phi_0","phi_max","lambda","gamma","sigma_inv","tau") #"gamma"
  
  
  library(rjags)
  fileConn=file("modelSig.tmp")
  #writeLines(modelGPV1,fileConn);
  writeLines(modelGCSigmoidInd,fileConn);
  close(fileConn)
  
  mLL=jags.model(file="modelSig.tmp",data=dataLL);
  update(mLL,burn_in);
  drawLL=jags.samples(mLL,steps,thin=thin,variable.names=varnames)
  #sampLL <- coda.samples(mLL,                      variable.names=varnames,                      n.iter=steps, progress.bar="none")
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_LL.pdf",sep="")) 
  plotGCSig(drawLL,dataLL,n=NA,groupID=2)
  #dev.off()
  #plotExpRes(drawLL,dataLL)
  
  
  
  
  #dev.off()
  ########################
  ## NL ###
  mNL=jags.model(file="modelSig.tmp",data=dataNL);
  update(mNL,burn_in)
  drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_NL.pdf",sep="")) 
  plotGCSig(drawNL,dataNL,n=NA,groupID=3)
  #dev.off()
  ############
  ### DL ###
  mDL=jags.model(file="modelSig.tmp",data=dataDL);
  update(mDL,burn_in)
  drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_DL.pdf",sep="")) 
  plotGCSig(drawDL,dataDL,n=NA,groupID=1)
  #dev.off()
  
  
  ## SHOW Hunt OnSet Densities ###
  ## Plot Vergence Angle Distributions/Densities
  dLLphi<-density(drawLL$tau)
  dDLphi<-density(drawDL$tau)
  dNLphi<-density(drawNL$tau)
  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_SigmoidFit_CompareOnset_tau.pdf",sep="")) 
  #X11()
  plot(dLLphi,col=colourH[2],type="l",lwd=2,ylim=c(0,1.0),main="Vergence Onset Vs Distance",xlab=expression(paste(" ",tau, " (mm)" ) ) )
  lines(dNLphi,col=colourH[3],lwd=2)
  lines(dDLphi,col=colourH[1],lwd=2)
  #legend("topright",legend =  strGroupID,fill=colourH)
  legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lnDat[["DL"]]),NROW(lnDat[["LL"]]) ,NROW(lnDat[["NL"]]) ) )
         ,fill=colourH,lty=c(1,2,3))
  dev.off()
  
  save.image(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit.RData",sep="") )
  
  #### CalcInformation ##
  load(paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit.RData",sep=""))
  
  
  #library("entropy")
  
  phi_hat <- function(x,Ulist){
    return(Ulist$phi_0  
          +(Ulist$phi_max - Ulist$phi_0)/( 1 + exp( -Ulist$gamma*( Ulist$tau - x)  ) )
          +Ulist$C*exp(Ulist$lambda*( Ulist$tau - x)) )
  }
  
  #Returns the FrequenciesAround mean (phihat)
  phiDens <- function(phi,x,Ulist)
  {
  
    return( dnorm(phi,mean=phi_hat(x,Ulist),sd=Ulist$sigma ) ) 
  }
  
  ##
  ## Calc Info In single Sample & Hunt Event
  InfoCalc <- function(DistMin,DistMax,Ulist)
  {
    PhiRange <- seq(0,90,1)
    DistRange <- seq(DistMin,DistMax,0.1)
    
    Grid <- expand.grid(PhiRange,DistRange)
    PVec=rep(0,NROW(Grid))
    
    for (i in 1:NROW(Grid) )
    {
      PVec[i] <- phiDens(Grid[i,1],Grid[i,2],Ulist)
    }
    
    PVec=PVec/sum(PVec)
    
    # Convert Pvec to a matrix 
    PMatrix=matrix(PVec,length(PhiRange),length(DistRange))
    
    image(PMatrix)
    MargVec=rowSums(PMatrix) ### Marginalize to obtain P(Response/Phi)
    Iloc=PMatrix/MargVec*length(DistRange) ##Information On Local x/For Each X - Assume X is unif. and so Prob[X]=1/Length(X)
    ###row sum
    sel=PMatrix>0
    #INFO=sums(PMatrix[sel]*log2(Iloc[sel]) )
    INFO=colSums(PMatrix*log2(Iloc) )
    return(INFO)
  }
    
  calcInfoOfHuntEvent <- function(drawLL)
  {
    Ulist <- list(phi_0=2,phi_max=35,gamma=100,tau=2,lambda=2,C=1,sigma=3)
    plot(InfoCalc(DistMin = 0.5,DistMax = 4,Ulist = Ulist))
    
    ##Next Is integrate Over sampled points, and Make Ii hidx matrix 
  }
  