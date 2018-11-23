## 23-10-2018
### Consider Information content in eye vergence on the distance to prey
## We use Model and the sampled parameter values to obtain an estimate of the mean 
##  information content in each hunt episode 
## We calculate the mutual information using P(Phi | X), assuming X is uniform

source("DataLabelling/labelHuntEvents_lib.r") ##for convertToScoreLabel
source("TrackerDataFilesImport_lib.r")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")

fitseqNo <- 5

#### CalcInformation ##
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_",fitseqNo,".RData",sep="") )

## Can Load Last Inf Matrix 
#load(file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit.RData",sep=""))
       

#### Load the Tracked Hunts Register ###
strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register",".rds",sep="") #Processed Registry on which we add 
message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 

strGroupID <- levels(datTrackedEventsRegister$groupID)

NSamples <-10 ## Number of Fit Lines (rows) to add in info matrrix


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.20,0.40,0.5,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.4),rgb(0.00,0.00,0.0,0.3)) ##Region (Transparency)
pchL <- c(16,2,4)
ltL <-  c(1,2,3) ##line Types

#library("entropy")
## Estimate Response Based On Model's mean response ## 
phi_hat <- function(x,Ulist){
  return(Ulist$phi_0  
         +(Ulist$phi_max - Ulist$phi_0)/( 1 + exp( -Ulist$lambda*( Ulist$tau - x)  ) )
         +Ulist$C*exp(Ulist$gamma*( Ulist$tau - x)) )
}

#Returns the FrequenciesAround mean (phihat)
phiDens <- function(phi,x,Ulist)
{
### Gaussian around mean response  ###
  return( dnorm(phi,mean=phi_hat(x,Ulist),sd=Ulist$sigma   ) ) #sqrt()
}

##
## Calc Info In single Sample & Hunt Event
InfoCalc <- function(DistRange,minDist,maxDist,Ulist)
{
  PhiRange <- seq(0,100,1) ##We limit The information Obtained To Reasonable Ranges Of Phi (Vergence Angle)
  
  
  #Generates a list with all possible pairs of in the range of Phi And Distance X 
  Grid <- expand.grid(PhiRange,DistRange)
  PVec=rep(0,NROW(Grid)) 
  ##Calc Density for all input Space
  for (i in 1:NROW(Grid) )
  {
    PVec[i] <- phiDens(Grid[i,1],Grid[i,2],Ulist)
  }
  ##Normalize to Probability
  PVec=PVec/sum(PVec)
  
  # Convert Pvec to a matrix / 
  PMatrix=matrix(PVec,nrow=length(PhiRange),ncol=length(DistRange),byrow = FALSE)
  
  ##Image shows a transpose Of the PMatrix, - In reality Rows contain Phi, and Cols are X 
  #image(t(PMatrix),y=PhiRange,x=DistRange)
  MargVec=rowSums(PMatrix) ### Marginalize across X to obtain P(Response/Phi)
  
  Iloc=PMatrix/MargVec*length(DistRange) ##Information On Local x/For Each X - Assume X is unif. and so Prob[X]=1/Length(X)
  
  ###row sum
  sel=PMatrix>0
  #INFO=sums(PMatrix[sel]*log2(Iloc[sel]) )
  ### Return Marginals I_xPhi
  mInfo <- PMatrix*log2(Iloc)
  
  ##Calc Marginals only on relevant Region - up to min distance from Target
  ##Make Binary Array indicating OutOf Bound regions
  vIntRegion<-as.numeric( DistRange >= minDist)
  mIntRegion <- t(matrix(vIntRegion,dim(mInfo)[2],dim(mInfo)[1]))
  
  INFO=colSums(mInfo*mIntRegion,na.rm=TRUE )
  
  return(INFO)
}


## Calculates Mutual Information X to Phi 
##  Fitted Function instances for each Hunt event
##  Assume distance is encoded in Dx steps, in N bits
## Returns a matrix Giving N Inf. Measures for each Hunt event
calcInfoOfHuntEvent <- function(drawS,dataSubset,n=NA,groupID)
{
  bPlot <- TRUE
  DistMin = 0.5
  DistMax =  5
  ##Assume X distance is encoded using 5 bits
  DistRange <- seq(DistMin,DistMax, (DistMax-DistMin )/(2^5-1)  )
  
 
  NHuntEvents <- NROW(unique(dataSubset$hidx) )
  
  vRegIdx <- unique(dataSubset$RegistrarIdx) ##Get Vector Of RegIdx That Associate with the sample Sequence
  vsampleP <- unique(dataSubset$hidx) # Obtain the seq id, so we can link to drawSamples
  
  if (is.na(n))
  {
    n <- NROW(unique(dataSubset$hidx))
    ##Note:Suffled and so the content of vsampleP maintains the hidx from where we can recover the RegIdx
    vsampleP <- sample(unique(dataSubset$hidx),n) 
  }
  
  
  vsub <- which (dataSubset$hidx %in% vsampleP)
  
  
  
  mInfMatrix <- matrix(nrow=NSamples,ncol=NROW(vsampleP) )
  vsampleRegisterIdx <- vector()
  
  hCnt <- 0
  idxChain <- 1
  for (h in vsampleP)
  {
    print(h)
    hCnt <- hCnt + 1
    ##Hold RegIDx For This Hunt Event 
    vsampleRegisterIdx[hCnt] <- vRegIdx[h] ##Save the Registry Idx So we can refer back to which Fish This belongs to 
    
    vphi_0_sub <- tail(drawS$phi_0[h,,idxChain],n=NSamples)
    vphi_max_sub <- tail(drawS$phi_max[h,,idxChain],n=NSamples)
    vgamma_sub <- tail(drawS$gamma[h,,idxChain],n=NSamples)
    vlambda_sub <- tail(drawS$lambda[h,,idxChain],n=NSamples)
    vtau_sub   <- tail(drawS$tau[h,,idxChain],n=NSamples)
    vsigma_sub <- tail(drawS$sigma[h,,idxChain],n=NSamples)
    valpha_sub <- tail(drawS$alpha[h,,idxChain],n=NSamples)
    
    ##Use the minimum Distance to prey, to obtain the effective integration Zone That information should be measured from 
    DistRangeEffective_min <- min(dataSubset$distP[dataSubset$hidx == h])
    DistRangeEffective_max <- max(dataSubset$distP[dataSubset$hidx == h])
    
    vPP <- which (dataSubset$hidx == h)
    #pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_",strGroupID[groupID],"_Sigmoid_",pp,".pdf",sep="")) 
    if (bPlot)
    {
      pdf(file= paste(strPlotExportPath,"/stat/stat_InfMeasure_EyeVsDistance5mm_",strGroupID[groupID],"_SigExp_",vRegIdx[h],".pdf",sep="")) 
      par(mar = c(5,5,2,5))
      plot(dataSubset$distP[vPP],dataSubset$phi[vPP],pch=19,xlim=c(0,max(DistRange)),ylim=c(0,90),
           main=paste(strGroupID[groupID],h," (",vRegIdx[h],")") , 
           bg=colourP[2],col=colourP[1],
           cex=0.5,
           xlab="Distance From Prey",ylab=expression(paste(Phi)) )  
      #Plot The Fitted Function
      vY  <-  mean(valpha_sub)*exp(mean(vgamma_sub)*( mean(vtau_sub)-DistRange) )+ mean(vphi_0_sub)   +  ( mean(vphi_max_sub) - mean(vphi_0_sub)  )/(1+exp( -( mean(vlambda_sub)   *( mean(vtau_sub) -DistRange )   ) ) ) 
      lines( DistRange ,vY,type="l",col=colourR[3],lwd=3)
    }
    
    lPlist <- list()
    for (i in 1:NSamples)
    {
      
      lPlist[[i]] <- list(phi_0=vphi_0_sub[i],phi_max=vphi_max_sub[i], gamma= vgamma_sub[i],
                    tau=vtau_sub[i],
                    lambda=vlambda_sub[i], 
                    C=valpha_sub[i] #vC_sub[i]
                    , sigma=vsigma_sub[i] )
      
      ##Information Is in The Sum / Marginal Across X ##
      ##Keep Distance Fixed accross so we are comparing against the same amount of information baseline, 
      ## from which mutual inf gives us the reduction in uncertainty
      ## Ex. We need log2(NROW(DistRange)) of bits to encode distance at the res. defined in DistRange
      
      vPPhiPerX <- InfoCalc(DistRange,DistRangeEffective_min,DistRangeEffective_max,Ulist = lPlist[[i]]) ##Get MutInf Per X

      ##Plot The Information Content Of The fitted Function ##
      if (bPlot)
      {
        par(new=T) 
          plot(DistRange,vPPhiPerX,ylim=c(0,2.5),xlim=c(0,max(DistRange) ),axes=F,type="p",pch=19,col=colourP[4], xlab=NA,
               ylab=NA,sub=paste("x requires ", round(100*log2(NROW(DistRange)) )/100,"bits"  ) )
         lines(DistRange,rev(cumsum(rev(vPPhiPerX))),ylim=c(0,2.5),xlim=c(0,max(DistRange) ),type="l",col=colourR[4], xlab=NA, ylab=NA )
      }
        
      mInfMatrix[i,hCnt] <- sum(vPPhiPerX) ## Mutual Inf Of X and Phi
    } 
    
    if (bPlot)
    {
      axis(side=4)
      mtext(side = 4, line = 3, 'Mutual information (bits)')
      dev.off()
    }
  }##For Each Hunt Event   
  
  
  ##Next Is integrate Over sampled points, and Make Ii hidx matrix 
  lInfoMatStruct <- list(infoMatrix=mInfMatrix,vsampleRegisterIdx=vsampleRegisterIdx,vsamplePSeqIdx=vsampleP)
  return(lInfoMatStruct)
}


##Returns Data frame With Mean Information in EyeVergence Along With Undershoot Ratio For each Tracked Event
mergeFirstTurnWithInformation <- function(datFirstBouts,lInfStruct)
{
  
  ##Subset the Event We have a first Turn Bout and Have Measured Eye INf (which are filtered for doing a strike attack)
  datBoutSubset <- datFirstBouts[datFirstBouts$RegistarIdx %in%  lInfStruct$vsampleRegisterIdx,]
  
  ## Retrieve Turn Bouts in right order to match the information Matric Column Order
  lboutSet <- list()
  idx <- 0
  for (idxR in lInfStruct$vsampleRegisterIdx)
  {
    idx <- idx +1
    if (NROW(datFirstBouts[datFirstBouts$RegistarIdx == idxR, ]) == 0 )
      next()
    
    rec <- datFirstBouts[datFirstBouts$RegistarIdx == idxR, ]
    vInfMeasure <- mean(lInfStruct$infoMatrix[,idx])
    
    lboutSet[[idx]] <- list(RegistarIdx=idxR,hidx= lInfStruct$vsamplePSeqIdx[idx], 
                            UnderShootAngle=abs( (rec$OnSetAngleToPrey)-abs(rec$Turn) ),
                            UnderShootRatio = (abs(rec$Turn)/abs(rec$OnSetAngleToPrey) ) ,
                            MInf=vInfMeasure) 
  }
  
  ##Combine Into A dataframe
  datFirstBoutVsInf<- data.frame(do.call(rbind,lboutSet ))
  
  
  return(datFirstBoutVsInf)
}

################# End of Library Functions ##### 

### Begin Script ###

## Sample Matrices Of Information / Retursn Struct Containing Mat and Id vectors ##
lInfStructLL <- calcInfoOfHuntEvent(drawLL,dataLL,groupID=2)
lInfStructNL <- calcInfoOfHuntEvent(drawNL,dataNL,groupID=3)
lInfStructDL <- calcInfoOfHuntEvent(drawDL,dataDL,groupID=1)

mInfMatrixLL <- lInfStructLL$infoMatrix
mInfMatrixNL <- lInfStructNL$infoMatrix
mInfMatrixDL <- lInfStructDL$infoMatrix


X11()
hist(mInfMatrixLL,col=colourH[2],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(mInfMatrixDL,col=colourH[1],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(mInfMatrixNL,col=colourH[3],xlim=c(0,3),breaks = seq(0,3,1/20))


X11()
hist(colMeans(mInfMatrixLL),col=colourH[2],xlim=c(0,3),breaks = seq(0,3,1/20),main="mean Inf Per Event ")
X11()
hist(colMeans(mInfMatrixDL),col=colourH[1],xlim=c(0,3),breaks = seq(0,3,1/20))
X11()
hist(colMeans(mInfMatrixNL),col=colourH[3],xlim=c(0,3),breaks = seq(0,3,1/20))


### Plot CDF ###
## Match the N 
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_CDF_4.pdf",sep=""))
subset_mInfMatrixLL <- mInfMatrixLL[,sample(1:58,58)]
plot(ecdf(mInfMatrixDL),col=colourH[1],main="Information In Eye Vergence CDF",
     xlab="Information (bits) ",lty=1,lwd=2,xlim=c(0,2.5))
plot(ecdf(mInfMatrixLL),col=colourH[2],add=T,lty=2,lwd=2)
plot(ecdf(mInfMatrixNL),col=colourH[3],add=T,lty=3,lwd=2)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()


##Plot Density 
dLLphi<-density(mInfMatrixLL)
dDLphi<-density(mInfMatrixDL)
dNLphi<-density(mInfMatrixNL)

#X11()
pdf(file= paste(strPlotExportPath,"/stat/stat_InfSigmoidExp_EyeVsDistance_Density_4.pdf",sep=""))
plot(dLLphi,col=colourH[2],type="l",lwd=3,lty=2,
     ylim=c(0,2),main="Mutual Information Distance to Eye Vergence ",
     xlab=expression(paste(" (bits)" ) ) )
lines(dNLphi,col=colourH[3],lwd=3,lty=3)
lines(dDLphi,col=colourH[1],lwd=3,lty=1)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NCOL(mInfMatrixDL),NCOL(subset_mInfMatrixLL) ,NCOL(mInfMatrixNL) ) ) 
       ,col=colourH,lty=c(1,2,3),lwd=2)
dev.off()

save(lInfStructLL,lInfStructDL,lInfStructNL,drawLL,drawDL,drawNL,file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit_4.RData",sep=""))      

##########
####  Examine Correation Of Information Vs Undershooting ###

##Load List the information measured ##
load(file=paste(strDataExportDir,"/stat_infoMat_EyeVergenceVsDistance_sigmoidFit5mm-5bit_4.RData",sep=""))

### Load Regression Data ###
load(paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJAgsOUt_",fitseqNo,".RData",sep=""))
## Load the First Bout turn data --
lFirstBoutPoints <- readRDS(file=paste(strDataExportDir,"/huntEpisodeAnalysis_FirstBoutData",".rds",sep="") ) #Processed Registry on which we add )

##Convert to data frame 
datFirstBouts <- data.frame(do.call(rbind,lFirstBoutPoints )) ##data.frame(lFirstBoutPoints[["LL"]] )
##Select The First Bout Relevant Records - 
#Note we are subsetting the events at the intersection of having a 1st bout turn to prey AND a strike to prey

##Make Data Frame For Debug of Events #
datTrackedEventsRegisterWithFirstBout <- datTrackedEventsRegister[datFirstBouts$RegistarIdx,]

datFirstBoutVsInfLL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructLL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfLL$RegistarIdx),"groupID"]) == "LL" )  ##Check for Errors in Reg idx - Group should match registry

datFirstBoutVsInfNL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructNL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfNL$RegistarIdx),"groupID"]) == "NL" )  ##Check for Errors in Reg idx - Group should match registry

datFirstBoutVsInfDL <- mergeFirstTurnWithInformation(datFirstBouts,lInfStructDL  )
## Error Check That subset of data using RegIdx Indeed belongs to the intended group
stopifnot(unique(datTrackedEventsRegister[unlist(datFirstBoutVsInfDL$RegistarIdx),"groupID"]) == "DL" )  ##Check for Errors in Reg idx - Group should match registry



## plot Undershot Ratio ###
pdf(file= paste(strPlotExportPath,"/stat/stat_InfVsTurnRatio.pdf",sep=""))
plot(unlist(datFirstBoutVsInfLL$UnderShootRatio),unlist(datFirstBoutVsInfLL$MInf),
     ylim=c(0,2),xlim=c(0,2),
     xlab=( expression(paste(Phi,"/",theta," Turn Ratio  ") )  ),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot ",col=colourH[2],pch=pchL[2])
segments(1,-10,1,20);

text(unlist(datFirstBoutVsInfLL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfLL$MInf)*1.01,unlist(datFirstBoutVsInfLL$RegistarIdx),cex=0.7)

points(unlist(datFirstBoutVsInfNL$UnderShootRatio),unlist(datFirstBoutVsInfNL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[3],pch=pchL[3])
text(unlist(datFirstBoutVsInfNL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfNL$MInf)*1.01,unlist(datFirstBoutVsInfNL$RegistarIdx),cex=0.7,col=colourP[3])

points(unlist(datFirstBoutVsInfDL$UnderShootRatio),unlist(datFirstBoutVsInfDL$MInf),ylim=c(0,2),xlim=c(0,2),
       col=colourH[1],pch=pchL[1])
text(unlist(datFirstBoutVsInfDL$UnderShootRatio)*1.01,unlist(datFirstBoutVsInfDL$MInf)*1.01,datFirstBoutVsInfDL$RegistarIdx,cex=0.7,col=colourP[1])

legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(datFirstBoutVsInfDL),NROW(datFirstBoutVsInfLL) ,NROW(datFirstBoutVsInfNL) ) ) 
       ,col=colourH,pch=pchL,lty=c(1,2,3),lwd=2)
dev.off()

##Plot Vs Undershoot Angle ###
plot(datFirstBoutVsInfLL$UnderShootAngle,datFirstBoutVsInfLL$MInf,
     ylim=c(0,2),xlim=c(-60,60),xlab=("OnSetAngleToPrey - Turn Angle"),ylab="mutual Inf in Eye V",
     main="Information Vs Undershoot ",pch=pchL[2])
points(datFirstBoutVsInfNL$UnderShootAngle,datFirstBoutVsInfNL$MInf,ylim=c(0,2),
       col="red",pch=pchL[3])
points(datFirstBoutVsInfDL$UnderShootAngle,datFirstBoutVsInfDL$MInf,ylim=c(0,2),
       col="blue",pch=pchL[1])



########








## Entropy ##
#if you use this package please cite: Jean Hausser and Korbinian Strimmer. 2009. Entropy inference
#and the James-Stein estimator, with application to nonlinear gene association networks.  J. Mach.
#Learn.  Res.10:  1469-1484.  Available online fromhttp://jmlr.csail.mit.edu/papers/v10/  hausser09a.html.
library("entropy")

KL.empirical(binNLphi,binLLphi,unit="log2")


mp <- barplot(c(entropy(binDLphi,unit="log2"),entropy(binLLphi,unit="log2"),entropy(binNLphi,unit="log2")),col=colourR,xlab="Group ",
              main=expression(paste("Entropy of eye vergence angle during hunt" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg" ) ) ,
              ylab="(bits)", ylim=c(0,8) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")


## Mutual Information ##
binLLphiVsDist <- discretize2d(dataLL$phi, dataLL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binNLphiVsDist <- discretize2d(dataNL$phi, dataNL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))
binDLphiVsDist <- discretize2d(dataDL$phi, dataDL$distP, numBins1=80, numBins2=10,r1=c(0,80),r2=c(0.5,2))

## Joint Entropy
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H12_NL = entropy(binNLphiVsDist ,unit="log2")
H12_DL = entropy(binDLphiVsDist,unit="log2" )
##Plot Mutual Information
mp <- barplot(c(H12_DL,H12_LL,H12_NL),col=colourR,xlab="Group ",
              main=expression(paste("Joint entropy prey distance and eye vergence" ) ),sub=expression(paste("Resolution: ", Phi,": 1 deg, Distance : 0.5mm" ) ) ,
              ylab="(bits)", ylim=c(0,10) ) #legend.text = strGroupID
text(mp,y=c(-0.3,-0.3,-0.3),labels=strGroupID , xpd = TRUE, col = "black")

# mutual information
mi.empirical(binLLphiVsDist,unit="log2" ) # //Mutual Informations
mi.empirical(binNLphiVsDist,unit="log2" ) # approximately zero
mi.empirical(binDLphiVsDist,unit="log2" ) # approximately zero

# another way to compute mutual information
# compute marginal entropies
H12_LL = entropy(binLLphiVsDist,unit="log2" )
H1_LL = entropy(rowSums(binLLphiVsDist),unit="log2")
H2_LL = entropy(colSums(binLLphiVsDist),unit="log2")
H1_LL+H2_LL-H12_LL # mutual entropy

H12_NL = entropy(binNLphiVsDist,unit="log2" )
H1_NL = entropy(rowSums(binNLphiVsDist),unit="log2")
H2_NL = entropy(colSums(binNLphiVsDist),unit="log2")
H1_NL+H2_NL-H12_NL # mutual entropy


H12_DL = entropy(binDLphiVsDist,unit="log2" )
H1_DL = entropy(rowSums(binDLphiVsDist),unit="log2")
H2_DL = entropy(colSums(binDLphiVsDist),unit="log2")
H1_DL+H2_DL-H12_DL # mutual entropy



# sample from continuous uniform distribution
x1 = runif(10000)
X11()
hist(x1, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories
y1 = discretize(x1, numBins=10, r=c(0,1))
# compute entropy from counts
entropy(y1) # empirical estimate near theoretical maximu
log(10)

#discretize2d(x1, x2, numBins1=10, numBins2=10)

# sample from a non-uniform distribution
X11()
hist(x2, xlim=c(0,1), freq=FALSE)
# discretize into 10 categories and estimate entropy

#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_LL_E.pdf",sep=""))
X11()
hist(drawLL$lambda[,,],col=colourR[2],main="Lambda")
hist(drawNL$lambda[,,],col=colourR[3],add=TRUE)
hist(drawDL$lambda[,,],col=colourR[1],add=TRUE)

X11()
hist(drawLL$gamma[,,],col=colourR[2],main="gamma")
hist(drawNL$gamma[,,],col=colourR[3],add=TRUE)
hist(drawDL$gamma[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_0[,,],col=colourR[2],main="Phi 0")
hist(drawNL$phi_0[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_0[,,],col=colourR[1],add=TRUE)


X11()
hist(drawLL$phi_max[,,],col=colourR[2],main="Phi Max")
hist(drawNL$phi_max[,,],col=colourR[3],add=TRUE)
hist(drawDL$phi_max[,,],col=colourR[1],add=TRUE)
#plot(drawLL$phi_max[3,,])

X11()
hist(drawLL$sigma[,1,,],main="LL")



X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawLL$u1[,,],breaks=50,xlim=c(0,7),col=colourH[2])
hist(drawLL$u0[,,],breaks=50,xlim=c(0,7),add=TRUE,col=colourH[2])

#dev.off()
########################
## NL ###
mNL=jags.model(file="model.tmp",data=dataNL);
drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)


X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_NL_F.pdf",sep=""))
plotGCRes(drawNL,dataNL,groupID=3)
#dev.off()

## Plot the infered function NL

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_Rate_lambda_NL_E.pdf",sep=""))
hist(drawNL$lambda[1,,1],main="NL")
#dev.off()

X11()
#pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_StartEnd_u0_NL_E.pdf",sep=""))
hist(drawNL$u1[1,,1],breaks=50,xlim=c(0,7),col="red")
hist(drawNL$u0[1,,1],breaks=50,xlim=c(0,7),add=TRUE,col="red")
#dev.off()





pdf(file= paste(strPlotExportPath,"/stat/stat_densityolinregressionslope.pdf",sep=""))
plot(dDLb,col=colourH[1],xlim=c(0.5,1.2),lwd=3,lty=1,ylim=c(0,20),main="Density Inference of Turn-To-Prey Slope ")
lines(dLLb,col=colourH[2],xlim=c(0.5,1.2),lwd=3,lty=2)

lines(dNLb,col=colourH[3],xlim=c(0.5,1.2),lwd=3,lty=3)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       ,fill=colourL,lty=c(1,2,3))
dev.off()

### PLot Scatter with regression lines with Conf intervals##
#X11()

pdf(file= paste(strPlotExportPath,"/stat/stat_TurnToPrey_LinearRegression.pdf",sep=""))
plot(lFirstBoutPoints[["DL"]][,1], lFirstBoutPoints[["DL"]][,2],
     main=paste("Turn Size Vs Bearing To Prey ", sep=""),
     xlab="Bearing To Prey prior to Bout",ylab="Bearing Change After Bout",xlim=c(-100,100),
     ylim=c(-100,100),
     col=colourP[1] ,pch=pchL[1]) ##boutSeq The order In Which The Occurred Coloured from Dark To Lighter
##Draw 0 Vertical Line
segments(0,-90,0,90); segments(-90,0,90,0); segments(-90,-90,90,90,lwd=1,lty=2);
#text(lFirstBoutPoints[["DL"]][,1]+2,lFirstBoutPoints[["DL"]][,2]+5,labels=lFirstBoutPoints[["DL"]][,3],cex=0.8,col="darkblue")
abline(lm(lFirstBoutPoints[["DL"]][,2] ~ lFirstBoutPoints[["DL"]][,1]),col=colourH[4],lwd=1.0) ##Fit Line / Regression
abline(a=muDLa,b=muDLb,col=colourH[1],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[1],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawDL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawDL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[1],lwd=4.0) ##Fit Line / Regression

#abline( lsfit(lFirstBoutPoints[["DL"]][,2], lFirstBoutPoints[["DL"]][,1] ) ,col=colourH[1],lwd=2.0)
##LL
points(lFirstBoutPoints[["LL"]][,1], lFirstBoutPoints[["LL"]][,2],pch=pchL[2],col=colourP[2])
#text(lFirstBoutPoints[["LL"]][,1]+2,lFirstBoutPoints[["LL"]][,2]+5,labels=lFirstBoutPoints[["LL"]][,3],cex=0.8,col="darkgreen")
abline(lm(lFirstBoutPoints[["LL"]][,2] ~ lFirstBoutPoints[["LL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muLLa,b=muLLb,col=colourH[2],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[2],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawLL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawLL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[2],lwd=4.0) ##Fit Line / Regression

#abline(lsfit(lFirstBoutPoints[["LL"]][,2], lFirstBoutPoints[["LL"]][,1] ) ,col=colourH[2],lwd=2.0)
##NL
points(lFirstBoutPoints[["NL"]][,1], lFirstBoutPoints[["NL"]][,2],pch=pchL[3],col=colourP[3])
#text(lFirstBoutPoints[["NL"]][,1]+2,lFirstBoutPoints[["NL"]][,2]+5,labels=lFirstBoutPoints[["NL"]][,3],cex=0.8,col="darkred")
abline(lm(lFirstBoutPoints[["NL"]][,2] ~ lFirstBoutPoints[["NL"]][,1]),col=colourH[4],lwd=1.0)
abline(a=muNLa,b=muNLb,col=colourH[3],lwd=1.5) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[2],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[2],col=colourR[3],lwd=4.0) ##Fit Line / Regression
abline(a=quantile(drawNL$beta[,(steps-ind):steps,1][1,])[3],b=quantile(drawNL$beta[,(steps-ind):steps,1][2,])[3],col=colourR[3],lwd=4.0) ##Fit Line / Regression
#abline( lsfit(lFirstBoutPoints[["NL"]][,2], lFirstBoutPoints[["NL"]][,1] ) ,col=colourH[3],lwd=2.0)
legend("topleft",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(lFirstBoutPoints[["DL"]][,1]),NROW(lFirstBoutPoints[["LL"]][,1]) ,NROW(lFirstBoutPoints[["NL"]][,1] ) ) )
       , pch=pchL,col=colourL)

dev.off()





##Plot Densities Summary

sampNL <- coda.samples(mNL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
sampDL <- coda.samples(mDL,                      variable.names=c("beta","sigma"),                      n.iter=20000, progress.bar="none")
X11()
plot(sampLL)
X11()
plot(sampNL)
X11()
plot(sampDL,main="DL")


