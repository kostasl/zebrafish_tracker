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
source("Stats/stat_SigmoidExp_EyeVsPreyDistance_lib.R")
### Hunting Episode Analysis ####
source("HuntingEventAnalysis_lib.r")
library(rjags)
library(runjags)


#
#These RC params Work Well to Smooth LF And NF
burn_in=1000;
steps=10000;
thin=3;
nchains <-3
n.cores <- 6
timings <- vector('numeric', 3)

dataFrac <- 1.0 ##Fraction Of Hunt Episodes to Include in DataSet
sampleFraction  <- 0.75 ##Fraction of Points to Use from Each Hunt Episode's data
fitseqNo <- 11
npad <- 1


##For the 3 Groups 
colourH <- c(rgb(0.01,0.01,0.9,0.8),rgb(0.01,0.7,0.01,0.8),rgb(0.9,0.01,0.01,0.8),rgb(0.00,0.00,0.0,1.0)) ##Legend
colourP <- c(rgb(0.01,0.01,0.8,0.5),rgb(0.01,0.6,0.01,0.5),rgb(0.8,0.01,0.01,0.5),rgb(0.00,0.00,0.0,1.0)) ##points DL,LL,NL
colourR <- c(rgb(0.01,0.01,0.9,0.4),rgb(0.01,0.7,0.01,0.4),rgb(0.9,0.01,0.01,0.3),rgb(0.00,0.00,0.0,1.0)) ##Region (Transparency)
pchL <- c(16,2,4)

  
  ####Select Subset Of Data To Analyse
  
  strRegisterDataFileName <- paste(strDataExportDir,"/setn_huntEventsTrackAnalysis_Register_SetB",".rds",sep="") #Processed Registry on which we add 
  message(paste(" Importing Retracked HuntEvents from:",strRegisterDataFileName))
  datTrackedEventsRegister <- readRDS(strRegisterDataFileName) ## THis is the Processed Register File On 
  
  lEyeMotionDat <- readRDS(paste(strDataExportDir,"/huntEpisodeAnalysis_EyeMotionData.rds",sep="") ) #Processed Registry on which we add )
  
  datEyeVsPreyCombinedAll <-  data.frame( do.call(rbind,lEyeMotionDat ) )
  
  ##datEyeVsPreyCombinedAll <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$doesCaptureStrike == 1,] # Select The Ones With a strike Capture
  
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
  lnInitPhi <- list() ##For Passing Prior for phi0
  ##Do all this processing to add a sequence index To The hunt Event + make vergence angle INdex 
  for (g in strGroupID) {
    lRegIdx[[g]] <- unique(datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == g),"RegistarIdx"])
   
    ldatLEyePoints[[g]] <- list()
    
    ##Take Each Hunting Event of Group And Select a Subset of the records. add to dataStruct List which will be fed to Model
    for (h in 1:NROW(lRegIdx[[g]]) )
    {
      ldatsubSet[[g]] <- datEyeVsPreyCombinedAll[datEyeVsPreyCombinedAll$groupID == which(strGroupID == g) &
                                                   datEyeVsPreyCombinedAll$RegistarIdx %in% lRegIdx[[g]][h] ,]  ##Select THe Ones With Capture
      ## Get Initial Eye V - Phi 0 
      recPhi0 <- head(ldatsubSet[[g]][ldatsubSet[[g]]$DistToPrey == ldatsubSet[[g]]$DistToPreyInit, ],1)
      
      
      ldatsubSet[[g]] <- ldatsubSet[[g]][sample(NROW(ldatsubSet[[g]]),sampleFraction*NROW(ldatsubSet[[g]] ) ) ,] ##SubSample Points 
      
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
        initVAngle=min(recPhi0$LEyeAngle-recPhi0$REyeAngle) ,
        RegistarIdx=ldatsubSet[[g]]$RegistarIdx,
        seqIdx=h)
      
      
      ##Augment with Idle phi entries for this hunt Event- to go up to 6 mm - 0 
      missingRegion <- 6 - head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1) 
      
      missingRegion <- 0 ##Do Not Pad The data
      
      if (missingRegion > 0)
      {
        ### Max Angle When Info Is Missing is 20
        datpadding <- cbind(vAngle=rep( max(0.1,min( c(20,ldatVEyePoints[[g]][[h]][,"vAngle"] ) ) )  ,npad),
                            distToPrey = seq(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),6,length=npad),
                            initDistToPrey = rep(head(as.numeric(ldatsubSet[[g]]$DistToPreyInit ),n=1),npad),
                            initVAngle = rep( ldatVEyePoints[[g]][[h]][,"initVAngle"] ,npad),
                            RegistarIdx = rep(head(as.numeric(ldatsubSet[[g]]$RegistarIdx ),n=1),npad),
                            seqIdx = h)
        
        ldatVEyePoints[[g]][[h]] <- rbind(ldatVEyePoints[[g]][[h]],datpadding) ##Add The Zero Phi Data
      }
      
      
      #ldatLEyePoints[[g]][[h]] <- ldatLEyePoints[[g]][[h]][!is.na(ldatLEyePoints[[g]][[h]][,2]),]
      #ldatREyePoints[[g]][[h]] <- ldatREyePoints[[g]][[h]][!is.na(ldatREyePoints[[g]][[h]][,2]),]
      #ldatVEyePoints[[g]][[h]] <- ldatVEyePoints[[g]][[h]][!is.na(ldatVEyePoints[[g]][[h]][,2]),]
      
      lnDat[[g]][[h]] <- NROW(ldatLEyePoints[[g]][[h]]) ##Not Used Anymore
      lnMaxDistanceToPrey[[g]][[h]] <- as.numeric(head(ldatsubSet[[g]]$DistToPreyInit,1)  ) ##Hold Unique Value Of Max Distance To Prey
      #lnInitPhi[[g]][[h]] <- as.numeric(head(ldatVEyePoints[[g]]$initVAngle,1)  ) ##Hold Unique Value Of Max Distance To Prey
    } ##For Each Hunt Event
    
  } ##For Each Group[]
  datVEyePointsLL <- data.frame( do.call(rbind,ldatVEyePoints[["LL"]] ) ) 
  datVEyePointsNL <- data.frame( do.call(rbind,ldatVEyePoints[["NL"]] ) ) 
  datVEyePointsDL <- data.frame( do.call(rbind,ldatVEyePoints[["DL"]] ) ) 
  
  
  ##SUBSET LL Filter Out Relevant Trajectories for Regression  #
  ## Best to focus on those that can fit the regressor ##
  vIdxExcludeAllExcept <-  unique(datVEyePointsLL[!(datVEyePointsLL$RegistarIdx %in% c(169,165)),]$RegistarIdx)  ##Ones that do not fit
  ##Exclude the ones Marked as unsuitable for the regressor
  vExcludesubIdx <- c(134,95,85,57,42,168,61) # 141
  ## OR Simply Subset Dat For Speed
  #vsubIdx <- sample(NROW(lRegIdx[["LL"]]),NROW(lRegIdx[["LL"]])*dataFrac)

  datVEyePointsLL_Sub <- datVEyePointsLL[!(datVEyePointsLL$RegistarIdx %in% vExcludesubIdx ),] #
  dataLL=list(phi=datVEyePointsLL_Sub$vAngle,
              distP=datVEyePointsLL_Sub$distToPrey ,
              N=NROW(datVEyePointsLL_Sub),
              distMax=lnMaxDistanceToPrey[["LL"]][unique(datVEyePointsLL_Sub$seqIdx)],##Re arrange the MaxDistances vector to the subset  ##lnMaxDistanceToPrey[["LL"]], #Put All distances in So We can Ref By Index #datVEyePointsLL_Sub$initDistToPrey,
              initPhi=unique(datVEyePointsLL_Sub$initVAngle),##Assumes Unique preserves order of rec appearance (seqIdx) ##Re arrange the MaxDistances vector to the subset  ##lnMaxDistanceToPrey[["LL"]], #Put All distances in So We can Ref By Index #datVEyePointsLL_Sub$initDistToPrey,
              ##Define an idx Key To link to the priors for each Eye Trajectory 
              hidx=as.numeric(factor(datVEyePointsLL_Sub$seqIdx)), ##Fix Seq to be 1...N over subset of data, 
              RegistrarIdx=datVEyePointsLL_Sub$RegistarIdx);
  
  
  
  ## NL  Subset Data ## 
  #vsubIdx <-sample(NROW(lRegIdx[["NL"]]),NROW(lRegIdx[["NL"]])*dataFrac)
  vIdxExcludeAllExcept <-  unique(datVEyePointsLL[!(datVEyePointsLL$RegistarIdx %in% c(154,152,130,104)),]$RegistarIdx)  ##Ones that do not fit
  ##Exclude the ones Marked as unsuitable for the regressor
  vExcludesubIdx <- c(147,112,111,109,100,98,44,26,22,20,19,13,12,11,10,152,99) 
  ## OR Simply Subset Dat For Speed
  #vsubIdx <- sample(NROW(lRegIdx[["LL"]]),NROW(lRegIdx[["LL"]])*dataFrac)
  datVEyePointsNL_Sub <- datVEyePointsNL[!(datVEyePointsNL$RegistarIdx %in% vExcludesubIdx ),] #
  dataNL=list(phi=datVEyePointsNL_Sub$vAngle,
              distP=datVEyePointsNL_Sub$distToPrey ,
              N=NROW(datVEyePointsNL_Sub),
              distMax=lnMaxDistanceToPrey[["NL"]][unique(datVEyePointsNL_Sub$seqIdx)],#datVEyePointsNL_Sub$initDistToPrey,
              initPhi=unique(datVEyePointsNL_Sub$initVAngle),##Assumes Unique preserves order of rec appearance (seqIdx) ##Re arrange the MaxDistances vector to the subset  ##lnMaxDistanceToPrey[["LL"]], #Put All distances in So We can Ref By Index #datVEyePointsLL_Sub$initDistToPrey,
              hidx=as.numeric(factor(datVEyePointsNL_Sub$seqIdx)),
              RegistrarIdx=datVEyePointsNL_Sub$RegistarIdx);
  
  ## DL SubSet Data
  #vsubIdx <-sample(NROW(lRegIdx[["DL"]]),NROW(lRegIdx[["DL"]])*dataFrac)
  vIdxExcludeAllExcept <-  unique(datVEyePointsLL[!(datVEyePointsLL$RegistarIdx %in% c(164,162,132,128)),]$RegistarIdx)  ##Ones that do not fit
  ##Exclude the ones Marked as unsuitable for the regressor
  vExcludesubIdx <- c(159,144,129,121,120,56,48,41,40,25,7)
  ## OR Simply Subset Dat For Speed
  #vsubIdx <- sample(NROW(lRegIdx[["LL"]]),NROW(lRegIdx[["LL"]])*dataFrac)
  datVEyePointsDL_Sub <- datVEyePointsDL[!(datVEyePointsDL$RegistarIdx %in% vExcludesubIdx ),] #
  dataDL=list(phi=datVEyePointsDL_Sub$vAngle,
              distP=datVEyePointsDL_Sub$distToPrey ,
              N=NROW(datVEyePointsDL_Sub),
              distMax=lnMaxDistanceToPrey[["DL"]][unique(datVEyePointsDL_Sub$seqIdx)],
              initPhi=unique(datVEyePointsDL_Sub$initVAngle),##Assumes Unique preserves order of rec appearance (seqIdx) ##Re arrange the MaxDistances vector to the subset  ##lnMaxDistanceToPrey[["LL"]], #Put All distances in So We can Ref By Index #datVEyePointsLL_Sub$initDistToPrey,
              hidx=as.numeric(factor(datVEyePointsDL_Sub$seqIdx)), ## Trick to reassign seq numbers of data subset
              RegistrarIdx=datVEyePointsDL_Sub$RegistarIdx);
  
  
  varnames=c("phi_0","phi_max","lambda","gamma","sigma","alpha","tau") #"gamma"
  fileConn=file("modelSig.tmp")
  #writeLines(modelGPV1,fileConn);
  writeLines(modelGCSigmoidInd,fileConn);
  close(fileConn)
  
  ### RUN METHOD 1 - CLASSIC RJAGS SAMPLES ###
  ### SETUP And Run THE LL Model Fit ###
  timer <- proc.time()
  mLL=jags.model(file="modelSig.tmp",
                 n.chains=nchains,
                 inits= initfunct(nchains, max(unique(dataLL$hidx) )),
                 data=dataLL);
  
  #update(mLL,burn_in);
  
  drawLL=jags.samples(mLL,
                      n.iter=steps,
                      thin=thin,
                      variable.names=varnames
  )
  time.taken <- proc.time() - timer
  timings[1] <- time.taken[3]/60
  
  
  ######### SAVE ##
  save(dataLL,drawLL,mLL,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_LL",fitseqNo,".RData",sep=""))      
  
  
  ## Plot The LL Fit ##
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_LL_C-",fitseqNo,".pdf",sep="")) 
  plotGCSig(drawLL,dataLL,n=NA,groupID=2)
  dev.off()
  
  
  ########################  N L  #########
  ## NL ### 
  timer <- proc.time()
  mNL=jags.model(file="modelSig.tmp",
                 n.chains=nchains,
                 inits= initfunct(nchains, NROW(unique(dataNL$hidx) )),
                 data=dataNL);
  #update(mNL,burn_in)
  drawNL=jags.samples(mNL,steps,thin=thin,variable.names=varnames)
  time.taken <- proc.time() - timer
  timings[2] <- time.taken[3]
  
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_NL_C-",fitseqNo,".pdf",sep="")) 
  plotGCSig(drawNL,dataNL,n=NA,groupID=3)
  dev.off()
  
  
  ######### SAVE ##
  save(dataNL,drawNL,mNL,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_NL",fitseqNo,".RData",sep=""))      
  
  
  ############ D L ######  ##
  ### DL ###
  timer <- proc.time()
  mDL=jags.model(file="modelSig.tmp",
                 n.chains=nchains,
                 inits= initfunct(nchains, NROW(unique(dataDL$hidx) )),
                 data=dataDL);
  #update(mDL,burn_in)
  drawDL=jags.samples(mDL,steps,thin=thin,variable.names=varnames)
  
  time.taken <- proc.time() - timer
  timings[3] <- time.taken[3]
  
  ######### SAVE ##
  save(dataDL,drawDL,mDL,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_DL",fitseqNo,".RData",sep=""))      
  
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_GroupSigmoidFit_DL_C-",fitseqNo,".pdf",sep="")) 
  plotGCSig(drawDL,dataDL,n=NA,groupID=1)
  dev.off()
  
  
  
  ### METHOD 2 - Using RUNJAGS ########
  ### #### Do it The Parallel Way RUN JAGS ##### ##
  ## Set up a distributed computing cluster with 2 nodes:
  library(parallel)
  library('coda')
  
  nchains <- 5
  n.cores <- nchains
  cl <- makeCluster(n.cores)
  timer <- proc.time()
  resultsLL <- run.jags(modelGCSigmoidInd,method = "rjparallel",
                        monitor = varnames,n.chains = nchains,
                        data=dataLL,thin = thin,
                        sample = steps,
                        inits = initfunct(nchains, NROW(unique(dataLL$hidx) )),
                        cl=cl)
  
  time.taken <- proc.time() - timer
  timings[1] <- time.taken[3]
  # Write the current model representation to file:
  write.jagsfile(resultsLL, file=paste(strDataExportDir,"/models/model_SigExpLL.txt",sep="" ) )
  mcmcLL.object <- as.mcmc.list(resultsLL)
  save(dataLL,resultsLL,timings,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RunJags_LL",fitseqNo,"-Parallel.RData",sep=""))      
  
  stopCluster(cl)
  
  
  ## NL    RUN JAGS ##
  cl <- makeCluster(n.cores)
  timer <- proc.time() 
  resultsNL <- run.jags(modelGCSigmoidInd,method = "rjparallel",
                         monitor = varnames,n.chains = nchains,
                         data=dataNL,
                         thin = thin,
                         sample = steps,
                         inits = initfunct(nchains, NROW(unique(dataNL$hidx) )),
                         cl=cl)
  time.taken <- proc.time() - timer
  timings[2] <- time.taken[3]
  # 
  # # Write the current model representation to file:
   write.jagsfile(resultsNL, file=paste(strDataExportDir,"/models/model_SigExpNL.txt",sep="" ) )
   save(dataNL,resultsNL,timings,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RunJags_NL",fitseqNo,"-Parallel.RData",sep=""))      
   mcmcNL.object <- as.mcmc.list(resultsNL)
   stopCluster(cl)

   ### DL RUN JAGS ##
   cl <- makeCluster(n.cores)
   timer <- proc.time() 
   resultsDL <- run.jags(modelGCSigmoidInd,method = "rjparallel",
                         monitor = varnames,n.chains = nchains,
                         data=dataDL,
                         thin = thin,
                         sample = steps,
                         inits = initfunct(nchains, NROW(unique(dataDL$hidx) )),
                         cl=cl)
   time.taken <- proc.time() - timer
   timings[3] <- time.taken[3]
   
  stopCluster(cl)
  
  # Write the current model representation to file:
  write.jagsfile(resultsDL, file=paste(strDataExportDir,"/models/model_SigExpDL.txt",sep="" ) )
  save(dataDL,resultsDL,timings,file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RunJags_DL",fitseqNo,"-Parallel.RData",sep=""))      
  mcmcDL.object <- as.mcmc.list(resultsDL)
  #######################
  ########################
  
  
  
  
  ############ Checking Convergence ############
  ## Compare convergence between the 3 chains for each trace 
  
  #### LOAD MCMC Data CalcInformation ##
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_LL",fitseqNo,".RData",sep=""))
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_NL",fitseqNo,".RData",sep=""))
load(file=paste(strDataExportDir,"/stat_EyeVergenceVsDistance_sigmoidFit_RJags_DL",fitseqNo,".RData",sep=""))
  
  
  idxH <- 1
  
  ## The gelman.diag gives you the scale reduction factors for each parameter.
  ## A factor of 1 means that between variance and within chain variance are equal, larger 
  # values mean that there is still a notable difference between chains. 
  
  dataS <- dataDL
  drawS <- drawDL
  strGroupID <- "DL"
  lFitScores_DL <- plotConvergenceDiagnostics(strGroupID,drawS, dataS)
  datFitScores_DL <- data.frame(do.call(rbind,lFitScores_DL))
  
  dataS <- dataNL
  drawS <- drawNL
  strGroupID <- "NL"
  lFitScores_NL <- plotConvergenceDiagnostics(strGroupID,drawS, dataS)
  datFitScores_NL <- data.frame(do.call(rbind,lFitScores_NL))
  
  dataS <- dataLL
  drawS <- drawLL
  strGroupID <- "LL"
  lFitScores_LL <- plotConvergenceDiagnostics(strGroupID,drawS, dataS)
  datFitScores_LL <- data.frame(do.call(rbind,lFitScores_LL))
  
  mErrorThres <- median(unlist(datFitScores_LL$meansqFitError) )
  datFitScores_LL[datFitScores_LL$meansqFitError < mErrorThres,]
  ####################
  
  
  ### List of Good Fit RegIdxs 
  vRegIdxGoodFit_LL <- c(3,4,5,23,18,29,30,31,37,38,49,50,53,54,58,62,66,68,69,70,74,75,76,77,79,83,84,94,113,114,133,135,137,138,139,140,155,165,166,167,169)
  vRegIdxGoodFit_DL <- c(1,2,8,14,16,17,18,24,27,55,65,72,73,78,80,86,118,119,127,128,131,145,157,158)
  vRegIdxGoodFit_NL <- c(15,21,32,33,34,35,36,39,45,46,60,81,87,88,89,93,101,102,103,104,105,107,108,110,115,116,123,124,125,126,130,143,148,149,150,151)
  ## Now Translate to dataIdx which We use to locate the  samples in the drawStructure
  vRegIdx <- unique(dataLL$RegistrarIdx)
  vIdxGoodFit_LL <- which(vRegIdx %in% vRegIdxGoodFit_LL)
  
  vRegIdx <- unique(dataDL$RegistrarIdx)
  vIdxGoodFit_DL <- which(vRegIdx %in% vRegIdxGoodFit_DL)
  
  vRegIdx <- unique(dataNL$RegistrarIdx)
  vIdxGoodFit_NL <- which(vRegIdx %in% vRegIdxGoodFit_NL)
  
  ### Examine Differences In Statistics of Curves
  ## SHOW Hunt OnSet Densities ###
  ## Plot Vergence Angle Distributions/Densities
  
   
  chain <-1 ##Which Chain To Draw The stats from 

  ### Compare Mean Onset Distance Tau 
  cdf_NL <- ecdf(unlist(drawNL$tau[vIdxGoodFit_NL,,chain]))
  cdf_DL <- ecdf(unlist(drawDL$tau[vIdxGoodFit_DL,,chain]))
  cdf_LL <-  ecdf(unlist(drawLL$tau[vIdxGoodFit_LL,,chain]))
  
  
  dLLphi<-density(unlist(drawLL$tau[vIdxGoodFit_LL,,chain]))
  dDLphi<-density(unlist(drawDL$tau[vIdxGoodFit_DL,,chain]))
  dNLphi<-density(unlist(drawNL$tau[vIdxGoodFit_NL,,chain]))

  
  #X11()
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_SigmoidFit_CompareOnset_tau2.pdf",sep=""),onefile=TRUE ) 
  #X11()
  plot(dLLphi,col=colourH[2],
       type="l",
       lwd=2,ylim=c(0,2.0),
       main="Vergence Onset densities",xlab=expression(paste(" ",tau, " (mm)" ) ) )
  lines(dNLphi,col=colourH[3],lwd=2)
  lines(dDLphi,col=colourH[1],lwd=2)
  #legend("topright",legend =  strGroupID,fill=colourH)
  legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(vIdxGoodFit_DL),NROW(vIdxGoodFit_LL) ,NROW(vIdxGoodFit_NL) ) )
         ,fill=colourH,lty=c(1,2,3))
  ###CDF Plot
  plot(cdf_NL,col=colourR[3],main="Onset Distance CDF")  
  lines(cdf_DL,col=colourR[1])
  lines(cdf_LL,col=colourR[2])
  
  
  dev.off()
  
  ############################################### ##
  ## Compare V Rise Rate after Vergences  Gamma  ###
  cdf_NL <- ecdf(unlist(drawNL$gamma[vIdxGoodFit_NL,,chain]))
  cdf_DL <- ecdf(unlist(drawDL$gamma[vIdxGoodFit_DL,,chain]))
  cdf_LL <-  ecdf(unlist(drawLL$gamma[vIdxGoodFit_LL,,chain]))
  
    
  dLLphi<-density(unlist(drawLL$gamma[vIdxGoodFit_LL,,chain]))
  dDLphi<-density(unlist(drawDL$gamma[vIdxGoodFit_DL,,chain]))
  dNLphi<-density(unlist(drawNL$gamma[vIdxGoodFit_NL,,chain]))
  
  ### Plot Gamma ###
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_SigmoidFit_CompareVRiseRate_gamma2.pdf",sep=""),onefile=TRUE ) 
  plot(dLLphi,col=colourH[2],type="l",lwd=2,
       ylim=c(0,1.5),main="Vergence Rise (Gamma) ",xlab=expression(paste(" ",gamma, " ()" ) ) )
  lines(dNLphi,col=colourH[3],lwd=2)
  lines(dDLphi,col=colourH[1],lwd=2)
  
    #legend("topright",legend =  strGroupID,fill=colourH)
  legend("topright",legend=paste(c("DL n=","LL n=","NL n="),c(NROW(vIdxGoodFit_DL),NROW(vIdxGoodFit_LL) ,NROW(vIdxGoodFit_NL) ) )
         ,fill=colourH,lty=c(1,2,3))
  
  ## Show Shift in Higher Rates Using CDF 
  ## Here NF seems to Lead the race, LF only wins by small margin
  plot(cdf_NL,col=colourR[3],main="gamma  Exp Rise Rate")  
  lines(cdf_DL,col=colourR[1])
  lines(cdf_LL,col=colourR[2])
  
  dev.off()
  
  ### PLot Gamma Vs Tau ###
  pdf(file= paste(strPlotExportPath,"/stat/stat_EyeVsDistance_SigmoidFit_gammaVStau.pdf",sep=""),onefile=TRUE ) 
  
  plot(unlist(drawLL$tau[vIdxGoodFit_LL,,chain]) , unlist(drawLL$gamma[vIdxGoodFit_LL,,chain]),col=colourR[2],
       main="Vergence rate (gamma) vs Onset distance tau ",
       xlab="distance (mm)",ylab="Rise rate (gamma)")  
  points(unlist(drawDL$tau[vIdxGoodFit_DL,,chain]), unlist(drawDL$gamma[vIdxGoodFit_DL,,chain]),col=colourR[1])  
  points(unlist(drawNL$tau[vIdxGoodFit_NL,,chain]), unlist(drawNL$gamma[vIdxGoodFit_NL,,chain]),col=colourR[3])  
  dev.off()
    
 ###
  ## Compare V Rise Rate after Vergences  Gamma  ###
  cdf_NL <- ecdf(unlist(drawNL$lambda[vIdxGoodFit_NL,,chain]))
  cdf_DL <- ecdf(unlist(drawDL$lambda[vIdxGoodFit_DL,,chain]))
  cdf_LL <-  ecdf(unlist(drawLL$lambda[vIdxGoodFit_LL,,chain]))
  
  ## Show Shift in Higher Rates Using CDF 
  ## Here NF seems to Lead the race, LF only wins by small margin
  plot(cdf_NL,col=colourR[3],main="lambda Hunting Onset Rise Rate")  
  lines(cdf_DL,col=colourR[1])
  lines(cdf_LL,col=colourR[2])
  
  #############
  ## Compare Max V Angle  ###
  cdf_NL <- ecdf(unlist(drawNL$phi_max[vIdxGoodFit_NL,,chain]))
  cdf_DL <- ecdf(unlist(drawDL$phi_max[vIdxGoodFit_DL,,chain]))
  cdf_LL <-  ecdf(unlist(drawLL$phi_max[vIdxGoodFit_LL,,chain]))
  
  ## Show Shift in Angles Rates Using CDF 
  ## Here NF seems to Lead the race, LF only wins by small margin
  plot(cdf_NL,col=colourR[3],main="Phi max Angle")  
  lines(cdf_DL,col=colourR[1])
  lines(cdf_LL,col=colourR[2])
  
  #############
  ## Compare Max V Angle  ###
  cdf_NL <- ecdf(unlist(drawNL$phi_0 [vIdxGoodFit_NL,,chain]))
  cdf_DL <- ecdf(unlist(drawDL$phi_0[vIdxGoodFit_DL,,chain]))
  cdf_LL <-  ecdf(unlist(drawLL$phi_0[vIdxGoodFit_LL,,chain]))
  
  ## Show Shift in Angles Rates Using CDF 
  ## Here NF seems to Lead the race, LF only wins by small margin
  plot(cdf_NL,col=colourR[3],main="Phi 0 Angle")  
  lines(cdf_DL,col=colourR[1])
  lines(cdf_LL,col=colourR[2])
  
  
    
  
  
  ## Load And RE-run from File ##
  # ByPass Silly path error for jags file
  # strWorkD <- getwd()
  # setwd("/")
  # #modelLL <- read.jagsfile(file=)
  # resultsLL <- run.jags(paste(strDataExportDir,"models/model_SigExpLL.txt",sep="" ),
  #                       data = dataLL,sample=100,adapt = 0,
  #                       keep.jags.files=paste(strDataExportDir,"/models/model_SigExpDL_results.txt",sep="" ) )
  # 
  # resultsNL <- run.jags(paste(strDataExportDir,"/models/model_SigExpNL.txt",sep="" ),sample=100,adapt = 0)
  # resultsDL <- run.jags(paste(strDataExportDir,"/models/model_SigExpDL.txt",sep="" ),sample=100)
  # 
  # 
  #X11()
  