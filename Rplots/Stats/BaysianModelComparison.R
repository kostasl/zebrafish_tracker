### Model Comparison ###

### Model Evidence ### 
## Compare Likelyhoods between models for undershooting (linear slope fit) 
## Establish if DryFed Data Belong to NF rather then LF
## 
## This was given by Giovanni , based on model evidence formula (see wikipedia Bayesian linear regression)

### Plot Undershoot Raw Data
hist(dataLL$turn/dataLL$bearing)
hist(dataNL$turn/dataNL$bearing)
hist(dataDL$turn/dataDL$bearing)

##Synthetic - TEst Data
x=seq(-10,10,0.1)
c11=0.5
c12=0.75
c21=0.5
c22=2.95

test_data1=cbind(x,c11+c12*x+rnorm(length(x)))
test_data2=cbind(x,c21+c22*x+rnorm(length(x)))
test_data3=rbind(test_data1,test_data2)

plot(test_data1)
points(test_data2,col="red")
###

getParams <- function(data,a0=1,b0=1,sigma0=1){
  n=nrow(data)
  y=data[,2]
  X=cbind(1,data[,1])
  Lambda0 = diag(sigma0,2)
  Lambda  = t(X)%*%X+Lambda0
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
  mu0 = c(1,0)
  mu  = solve(t(X)%*%X+Lambda0)%*%(t(X)%*%X%*%beta_hat+Lambda0%*%mu0)
  a=a0+n/2
  b=b0+0.5*(t(y)%*%y+t(mu0)%*%Lambda0%*%mu0-t(mu)%*%Lambda%*%mu)
  
  return(list(n=n,a0=a0,b0=b0,sigma0=sigma0,a=a,b=b,mu=mu,lambda=Lambda))
}

logML <- function(par){
  res=-par$n/2*log(2*pi)+0.5*log(det(diag(par$sigma0,2))/det( par$lambda)) + par$a0*log(par$b0) - par$a*log(par$b) + lgamma(par$a) - lgamma(par$a0)
  return(res)
}

## Example ##
p1=getParams(test_data1); lML1=logML(p1)
p2=getParams(test_data2); lML2=logML(p2)
p3=getParams(test_data3); lML3=logML(p3)
## +ve means they are from separate Sources
## -ve is interpreted as common model best describes these
logR=(lML1+lML2)-lML3


b0=1
a0=1
MLparamsLL <- getParams( cbind(dataLL$turn,dataLL$bearing),a0,b0 )
MLparamsDL <- getParams( cbind(dataDL$turn,dataDL$bearing),a0,b0 )
MLparamsNL <- getParams( cbind(dataNL$turn,dataNL$bearing),a0,b0 )

dataNLDL <- rbind(cbind(dataNL$turn,dataNL$bearing),cbind(dataDL$turn,dataDL$bearing))
MLparamsNLDL <- getParams( dataNLDL,a0,b0 )

dataDLLL <- rbind(cbind(dataDL$turn,dataDL$bearing),cbind(dataLL$turn,dataLL$bearing))
MLparamsDLLL <- getParams( dataDLLL,a0,b0 )

dataLLNL <- rbind(cbind(dataNL$turn,dataNL$bearing),cbind(dataLL$turn,dataLL$bearing))
MLparamsLLNL <- getParams( dataLLNL,a0,b0 )


## Calcilate Probability of Model Given Data
logML_LL <- logML(MLparamsLL)
logML_DL <- logML(MLparamsDL)
logML_NL <- logML(MLparamsNL)
logML_NLDL <- logML(MLparamsNLDL)
logML_DLLL <- logML(MLparamsDLLL)
logML_LLNL <- logML(MLparamsLLNL)

## Compare Models On Log Likehoods
## +ve means they are from separate Sources
## -ve is interpreted as common model best describes these
logR_DLNL=(logML_DL+logML_NL)-logML_NLDL
logR_LLDL=(logML_DL+logML_LL)-logML_DLLL
logR_LLNL=(logML_NL+logML_LL)-logML_LLNL



####### OLD ###
getParams_old <- function(data,a0=1,b0=1,sigma0=1){
  n=nrow(data)
  y=data[,2]
  X=cbind(1,data[,1])
  Lambda0 = diag(sigma0,2)
  Lambda  = t(X)%*%X+Lambda0 
  beta_hat = solve(t(X)%*%X) %*%t(X)%*%y
  mu0 = c(1,0)               
  mu  = solve(t(X)%*%X+Lambda0)%*%(t(X)%*%X%*%beta_hat+Lambda0%*%mu0)
  a=a0+n/2                           
  b=b0+0.5*(t(y)%*%y+t(mu0)%*%Lambda0%*%mu0-t(mu)%*%Lambda%*%mu)
  
  return(list(n=n,a=a,b=b,mu=mu,lambda=Lambda))  
}

##Marginal Likelyhood 
MarginalLikelihood <- function(MLParams,a0,b0)
{
  return (1/(2*pi)^(MLParams$n/2))* sqrt( det(diag(sigma0,2))/det( MLParams$lambda))*((b0^a0)/(MLParams$b^MLParams$a)) *(gamma(MLParams$a)/gamma(a
}



b0=1
a0=1
MLparamsLL <- getParams( cbind(dataLL$turn,dataLL$bearing),a0,b0 )
MLparamsDL <- getParams( cbind(dataDL$turn,dataDL$bearing),a0,b0 )
MLparamsNL <- getParams( cbind(dataNL$turn,dataNL$bearing),a0,b0 )

dataNLDL <- rbind(cbind(dataNL$turn,dataNL$bearing),cbind(dataDL$turn,dataDL$bearing))
MLparamsNLDL <- getParams( dataNLDL,a0,b0 )

dataDLLL <- rbind(cbind(dataDL$turn,dataDL$bearing),cbind(dataLL$turn,dataLL$bearing))
MLparamsDLLL <- getParams( dataDLLL,a0,b0 )

## Calcilate Probability of Model Given Data
ML_LL <- MarginalLikelihood(MLparamsLL,a0,b0)
ML_DL <- MarginalLikelihood(MLparamsDL,a0,b0)
ML_NL <- MarginalLikelihood(MLparamsNL,a0,b0)
ML_NLDL <- MarginalLikelihood(MLparamsNLDL,a0,b0)
ML_DLLL <- MarginalLikelihood(MLparamsDLLL,a0,b0)

##Now Compare ##
# A value of K > 1 means that M1 is more strongly supported by the data under consideration than M2.
ML_DL*ML_NL/(ML_NLDL) ## This is equal to 1 
ML_NL/(ML_NLDL)
ML_DL/(ML_NLDL)

## Check For COmparing DL LL 
ML_DL*ML_LL/(ML_DLLL)
ML_LL/(ML_DLLL)
ML_DL/(ML_DLLL)

mean(dataDL$turn/dataDL$bearing)
mean(dataNL$turn/dataNL$bearing)
