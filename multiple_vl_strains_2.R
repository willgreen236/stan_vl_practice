library(rstan)

N <- 3000
K <- 30
S <- 2
t <- rep(seq(0,99,1),K)
set.seed(1)
a <- c(rnorm(n=K, mean=0.55, sd=0.05))
b <- rnorm(n=K, mean=0.1, sd=0.02)
t_onset <- rnorm(n=K, mean=20, sd=5)
tmax <- rnorm(n=K, mean=30, sd=4)
strain <- c(rep(1, N/2), rep(2, N/2))

individual <- vector(length=N)
for (i in 1:length(individual)) individual[i] = floor((i-1)/(N/K) + 1)
  
vl <- vector(length = N)
for (i in 1:N) vl[i] = max(rnorm(n=1,(a[floor((i-1)/(N/K))+1]+b[floor((i-1)/(N/K))+1])/(b[floor((i-1)/(N/K))+1]*exp(-a[floor((i-1)/(N/K))+1]*(t[i]+t_onset[floor((i-1)/(N/K))+1]-tmax[floor((i-1)/(N/K))+1])) + a[floor((i-1)/(N/K))+1] * exp(b[floor((i-1)/(N/K))+1]*(t[i]+t_onset[floor((i-1)/(N/K))+1]-tmax[floor((i-1)/(N/K))+1]))), 0.01),0)

model_multiple_vl <- stan_model("multiple_vl_strains_2.stan")

initfun <- function(...) {
  list(a_bar=0.1, a_sigma = 0.1, a=rep(0.1,30),a_d_bar=0.1,a_d_sigma=0.1, b_bar=0.05, b_sigma = 0.01, b=rep(0.05, 30), tmax = rep(21,30), t_bar = 25, t_sigma = 4, sigma=0.0001)
}

fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), init=initfun, iter=1000, chains=1, seed=1)

print(fit)



vl_func <- function(a, b, tmax, t_onset, t){
  return((a+b)/(b*exp(-a*(t+t_onset-tmax))+a*exp(b*(t+t_onset-tmax))))
}

vl_test <- data.frame(t = seq(-20, 40), 
                      vl = vl_func(a=0.3, b=0.1, tmax = 30, t_onset=20, t = seq(-20:40)))

plot(vl_test$t, vl_test$vl)