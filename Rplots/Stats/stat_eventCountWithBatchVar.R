##Done With Giovanni on 14/2/2018 to analyse why significance was lost after adding the 08-02-2018 |Exp Set ##
##### Jags Baysian Statistics ###
## Testing How Variance Increases With Each Exp. Set ###
## PLoting With Increasing Datapoints reveals Effect on NF Live group coming from 08-02-2018 data ###
## This may suggest that either A mistake occured on the 08-02-2018 NF fish, 
## or that our Poisson Single rate prior is wrong   ## 


modelG="model { 
q ~ dnorm(15,0.0001)T(0,400)

for(j in 1:NTOT){
n[j] ~ dpois(q)
}


}"

modelI="model {
  qmean ~ dnorm(10,0.001)T(0,400)
  tt ~ dnorm(0,0.001)T(0,100)
  for(j in 1:NDATASET){
    q[j] ~ dnorm(qmean,tt)
  }

  for(j in 1:NTOT){
    n[j] ~ dpois(q[datasetID[j]])
  }
}"

fromchain=1000

datHuntVsPreyNL <- cbind(datHuntStat[,"vHInitialPreyCount"]$NL , as.numeric(datHuntStat[,"vHLarvaEventCount"]$NL) )

nDat <- nrow(datHuntVsPreyNL)
nDatList <- (nDat-10):nDat
datSliceNL <- (datHuntVsPreyNL[,2])

library(rjags)
varnames1=c("q","qmean")
burn_in=1000;
steps=30000;
thin=10;

strModelName = "modelI.tmp"
fileConn=file(strModelName)
writeLines(modelI,fileConn);
close(fileConn)

nNLlist<-vector("list",length(nDatList))
dataNLlist<-vector("list",length(nDatList))
mNLjm=vector("list",length(nDatList))
draw=vector("list",length(nDatList))
datasetIDlist=vector("list",length(nDatList))

for(i in 1:length(nDatList)){
  datasetIDlist[[i]]=sort(rep(1:20,4))[1:nDatList[i]]
  nNLlist[[i]]=as.numeric(datSliceNL[1:nDatList[i]])
  dataNLlist[[i]]=list(n=nNLlist[[i]],NTOT=nDatList[i],datasetID=datasetIDlist[[i]],NDATASET=max(datasetIDlist[[i]]));
  mNLjm[[i]]=jags.model(file=strModelName,data=dataNLlist[[i]]);
  update(mNLjm[[i]],burn_in)
  draw[[i]]=jags.samples(mNLjm[[i]],steps,thin=thin,variable.names=varnames1)
}

##q[idx,sampleID,chainID]


layout(matrix(1:(2*ceiling(length(nDatList)/2)),ncol=2))
for(i in 1:length(nDatList)) hist(draw[[i]]$qmean[1,(steps/10-fromchain):(steps/10),1],
                   breaks=seq(0,50,length=100),
                   col=gray.colors(4)[i],
                   xlab="Rate Parameter",
                   add=F,main=i)



