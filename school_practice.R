#load libraries
library(rstan)
library(RColorBrewer)

#simulate some data
set.seed(20161110)
N<-100                 #sample size
J<-10                  #number of plant species
id<-rep(1:J,each=10)   #index of plant species
K<-3                   #number of regression coefficients


gamma <- c(2,-1,3)  #population-level regression coefficient
tau   <- c(0.3,2,1) #standard deviation of the group-level coefficient
sigma <- 1          #standard deviation of individual observations
beta<-mapply(function(g,t) rnorm(J,g,t),g=gamma,t=tau) #group-level regression coefficients

#the model matrix
X<-model.matrix(~x+y,data=data.frame(x=runif(N,-2,2),y=runif(N,-2,2)))
y<-vector(length = N)
for(n in 1:N){
  #simulate response data
  y[n]<-rnorm(1,X[n,]%*%beta[id[n],],sigma)
}

#run the model
m_hier<-stan(file="school_practice.stan",data=list(N=N,J=J,K=K,id=id,X=X,y=y))

print(m_hier)
