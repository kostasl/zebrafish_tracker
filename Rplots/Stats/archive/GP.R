model="model {
  # Likelihood
  
  n ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:N) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(tau, 2)+pow(tau0,2)
  
    for(j in (i+1):N) {
      Sigma[i,j] <- pow(tau,2) * exp( - 0.5*rho^2 * pow(food[i] - food[j], 2) )
      Sigma[j,i] <- Sigma[i,j]
    }
  }
 
  alpha=0 
  tau0 ~ dgamma(2,0.2) 
  tau  ~ dgamma(2,.2) 
  rho = 0.1
  
}"

#load("data/setn-12-D-5-16-datHuntStat.RData")
foodlevels=datHuntStat[,"vHInitialPreyCount"]$LL
ord=order(foodlevels)
foodlevels=foodlevels[ord]
counts=as.numeric(datHuntStat[,"vHLarvaEventCount"]$LL)[ord]
data=list(n=counts,food=foodlevels,N=length(counts));
varnames=c("tau","rho","alpha","tau0")
burn_in=100;
steps=1000;
thin=1;

library(rjags)
fileConn=file("model.tmp")
writeLines(model,fileConn);
close(fileConn)

m=jags.model(file="model.tmp",data=data);
update(m,burn_in)
draw=jags.samples(m,steps,thin=thin,variable.names=varnames)

SE <- function(Xi,Xj, rho,tau) tau^2*exp(-0.5*(Xi - Xj) ^ 2 * rho^2)
cov <- function(X, Y, rho,tau) outer(X, Y, SE, rho,tau)

plot_res<- function(ind,qq=0.05){
	plot(foodlevels,counts,col="red")
	x_predict=seq(1,40,0.1)
	Ef=matrix(NA,ncol=length(x_predict),nrow=ind)
	for(j in 1:ind){
		i=steps-j+1
		cov_xx_inv=solve(cov(foodlevels,foodlevels,draw$rho[1,i,1],draw$tau[1,i,1])+diag(draw$tau0[1,i,1],length(foodlevels)))
		Ef[j,] <- cov(x_predict, foodlevels,draw$rho[1,i,1],draw$tau[1,i,1]) %*% cov_xx_inv %*% counts 
	}
	mu=apply(Ef,2,mean)
	band=apply(Ef,2,quantile,probs=c(qq,1-qq))
	lines(x_predict,mu,lwd=2)
	polygon(c(x_predict,rev(x_predict)),c(band[1,],rev(band[2,])),col="green")
}

