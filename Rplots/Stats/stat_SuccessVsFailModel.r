## stat Model Success Vs Failure - Plots Distribution For HuntRate and Success Rate
# Makes A model Assuming Poisson Rate lambda and Hunt Success Probability q
# Note: We can exclude the fish that produced no events when excluding eventID == 0- 

source("TrackerDataFilesImport_lib.r")
source("DataLabelling/labelHuntEvents_lib.r")


## The mixture Of Poisson Drawing from Gammaa, gives a negative binomial
modelNBinom="model {
         
         for(i in 1:3) {
             #lambda[i] ~ dgamma(1,1) ## 
             r[i] ~ dgamma(1,1) ##
             q[i] ~ dunif(0.0,1)

             f[i] ~ dbeta(1,1)
             p[i] ~ dbeta(1,1)
             t[i] ~ dbeta(1,1) ##Prob Of Enganging With Prey Given HuntMode Is On
         }

         for(i in 1:NTOT){
             Events[i] ~  dnegbin(q[ID[i]],r[ID[i]] )
             TrackPrey[i] ~ dbinom(t[ID[i]],Events[i])
             Success[i] ~ dbinom(p[ID[i]],TrackPrey[i])
             Fail[i] ~ dbinom(f[ID[i]],TrackPrey[i])
             
         }
}"




## Model (Wrongly) Assumes Single Distribution Prior For Rates - Yet, fitting the group event rates shows
## that it is best to assume a mixture of hunt rates exist in each group, such that the group hunt rates appears neg Binom.
## wiki suggest gamma prior for lambda
modelPoisson="model {
         
         for(i in 1:3) {
             lambda[i] ~ dgamma(1,1) ##Suggested in wiki 
             #lambda[i] ~ dexp(1) 
             q[i] ~ dbeta(1,1)
             p[i] ~ dbeta(1,1)
             t[i] ~ dbeta(1,1) ##Prob Of Enganging With Prey Given HuntMode Is On
         }

         for(i in 1:NTOT){
             Events[i] ~ dpois(lambda[ID[i]])
             TrackPrey[i] ~ dbinom(t[ID[i]],Events[i])
             Success[i] ~ dbinom(q[ID[i]],TrackPrey[i])
             Fail[i] ~ dbinom(p[ID[i]],TrackPrey[i])
             
         }
}"


#strProcDataFileName <- "setn14-HuntEventsFixExpID-SB-Updated"
#strProcDataFileName <-paste("setn-12-HuntEvents-SB-ALL_19-07-18",sep="") ## Latest Updated HuntEvent Labelled data
strProcDataFileName <- "setn15-HuntEvents-SB-Updated-Merged"
message(paste(" Loading Hunt Event List to Analyse... "))
#load(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".RData",sep="" )) ##Save With Dataset Idx Identifier
datHuntLabelledEventsSB <- readRDS(file=paste(strDatDir,"/LabelledSet/",strProcDataFileName,".rds",sep="" ))

##We Can Choose To Exclude The Fish That Produced No Hunting Events
datHuntLabelledEventsSB <- datHuntLabelledEventsSB[datHuntLabelledEventsSB$eventID != 0 &
                                                     datHuntLabelledEventsSB$groupID %in% c("LL","NL","DL") ,]
datFishSuccessRate <- getHuntSuccessPerFish(datHuntLabelledEventsSB)
datFishSuccessRate$groupID <- factor(datFishSuccessRate$groupID)
strGroups <-levels(datFishSuccessRate$groupID)

NRecCount_DL <- table(datFishSuccessRate$groupID)["DL"]
NRecCount_NL <- table(datFishSuccessRate$groupID)["NL"]
NRecCount_LL <- table(datFishSuccessRate$groupID)["LL"]

datatest=list(Success=datFishSuccessRate$Success,
              Fail=datFishSuccessRate$Fails,
              TrackPrey=datFishSuccessRate$Fails+datFishSuccessRate$Success, ##Number of Prey Engangements
              Events=datFishSuccessRate$HuntEvents, ##Includes No_Target - ie cases where stimuli Trigger Fish HuntMode But no Prey Tracking Seems to take place
              ID=as.numeric(datFishSuccessRate$groupID),
              NTOT=nrow(datFishSuccessRate));

varnames1=c("q","p","t","f","r")
burn_in=1000;
steps=100000;
thin=10;
chains=3

library(rjags)
strModelName = "model1.tmp"
fileConn=file(strModelName)
writeLines(modelNBinom,fileConn);
close(fileConn)

m=jags.model(file=strModelName,data=datatest,n.chains=chains);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames1)


#######PLOT RESULTS
#X11()
colourH <- c("#0303E663","#03B30363","#E6030363")
colourD <- c("#0303E623","#03B30323","#E6030323")
colourL <- c("#0303E6AF","#03B303AF","#E60303AF")
#for(i in 1:3) hist(draw$q[i,,1],breaks=seq(0,0.5,0.01),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

#X11()
#for(i in 1:3) hist(draw$lambda[i,,1],breaks=seq(5,20,0.2),col=colourH[i],add=!(i==1))
#legend("topright", legend=strGroups,fill=colourL)

##Density in 2D of Success Vs Fail
#X11()

strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation_Success.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of success, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))

nlevels <- 5
zLL <- kde2d(draw$q[2,,1], draw$lambda[2,,1], n=80)
zNL <- kde2d(draw$q[3,,1], draw$lambda[3,,1], n=80)
zDL <- kde2d(draw$q[1,,1], draw$lambda[1,,1], n=80)


Range_ylim <- c(7,17)
  plot(draw$q[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=Range_ylim,xlim=c(0,0.5),pch=19,
       main="Bayesian Estimation for Hunt Rate and Efficiency",
       xlab="Probability of Success q",
       ylab=(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
  contour(zDL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  points(draw$q[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=Range_ylim,xlim=c(0.1,0.5),pch=19)
  contour(zLL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  points(draw$q[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=Range_ylim,xlim=c(0.1,0.5),pch=19)
  contour(zNL, drawlabels=FALSE, nlevels=nlevels,add=TRUE)
  legend("topright", legend=paste(strGroups," n=",c(NRecCount_DL,NRecCount_LL,NRecCount_NL)),fill=colourL)
dev.off()


strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation_Fails.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of Engaging with Prey and Failing, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
plot(draw$p[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=Range_ylim,xlim=c(0,1),pch=19,
     main="Bayesian Estimation for Hunt Rate and Efficiency (Fails)",
     xlab="Probability of Failing p",
     ylab=(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
points(draw$p[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=Range_ylim,xlim=c(0.1,1),pch=19)
points(draw$p[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=Range_ylim,xlim=c(0.1,1),pch=19)
legend("topright", legend=strGroups,fill=colourL)
dev.off()

strPlotName = paste(strPlotExportPath,"/stat/stat_HuntRateAndEfficiencyEstimation_Track.pdf",sep="")
pdf(strPlotName,width=8,height=8,title="Bayesian Inference on distribution of hunt rate parameter and probability of Engaging with Prey, based on labelled data set",onefile = TRUE) #col=(as.integer(filtereddatAllFrames$expID))
plot(draw$t[1,,1], draw$lambda[1,,1],col=colourD[1],ylim=c(5,26),xlim=c(0,1),pch=19,
     main="Bayesian Estimation for Hunt Rate and Tracking Efficiency / Engaging with prey ",
     xlab="Probability of Tracking Prey t",
     ylab=(expression(paste("Hunt Rate ",lambda ) ) )  ) #paste("Hunt Rate", )
points(draw$t[2,,1], draw$lambda[2,,1],col=colourD[2],ylim=c(5,26),xlim=c(0,1),pch=19)
points(draw$t[3,,1], draw$lambda[3,,1],col=colourD[3],ylim=c(5,26),xlim=c(0,1),pch=19)
legend("topright", legend=strGroups,fill=colourL)
dev.off()

#hist(draw$q[1,,1],breaks=seq(0,30,length=100),col=colourH[1])
#hist(drawNL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[2])
#hist(drawDL$q[1,,1],breaks=seq(0,30,length=100),add=T,col=colourH[3])
#legend(2,500,legend = c("LL","NL","DL"),fill=c(colourH[1],colourH[2],colourH[3]))

#hist(drawLL$qq[1,,1],breaks=seq(0,50,length=100))
#hist(drawNL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(1,0,0,.4))
#hist(drawDL$qq[1,,1],breaks=seq(0,50,length=100),add=T,col=rgb(0,1,0,.4))

