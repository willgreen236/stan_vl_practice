library(rstan)

N <- 1000
K <- 10
t <- rep(seq(0,99,1),K)
set.seed(1)
a <- c(rnorm(n=K/2, mean=0.55, sd=0.05), rnorm(n=K/2, mean=0.45, sd=0.05))
b <- rnorm(n=K, mean=0.1, sd=0.02)
tmax <- rnorm(n=K, mean=30, sd=4)
strain <- c(rep(1, N/2), rep(2, N/2))

individual <- vector(length=N)
for (i in 1:length(individual)) individual[i] = floor((i-1)/(N/K) + 1)
  
vl <- vector(length = N)
for (i in 1:N) vl[i] = max(rnorm(n=1,(a[floor((i-1)/(N/K))+1]+b[floor((i-1)/(N/K))+1])/(b[floor((i-1)/(N/K))+1]*exp(-a[floor((i-1)/(N/K))+1]*(t[i]-tmax[floor((i-1)/(N/K))+1])) + a[floor((i-1)/(N/K))+1] * exp(b[floor((i-1)/(N/K))+1]*(t[i]-tmax[floor((i-1)/(N/K))+1]))), 0.01),0)

plot(vl[1:300])

model_multiple_vl <- stan_model("multiple_vl.stan")

initfun <- function(...) {
  list(a_bar=0.1, a_sigma = 0.1, a=rep(0.1,10), b_bar=0.05, b_sigma = 0.01, b=rep(0.05, 10), tmax = rep(21,10), t_bar = 25, t_sigma = 4, sigma=0.0001)
}

fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), iter=10000, chains=1, seed=1)
fit2 <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), init=initfun, iter=1000, chains=1, seed=1)

extract(model_multiple_vl)

print(fit)
print(fit2)
