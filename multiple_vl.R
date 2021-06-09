library(rstan)

N <- 300
K <- 3
t <- rep(seq(0,99,1),3)
a <- c(0.6,0.5,0.4)
b <- c(0.1,0.1, 0.1)
tmax <- c(30,30,30)
individual <- c(rep(1,100),rep(2,100),rep(3,100))

vl <- vector(length = N)
for (i in 1:N) vl[i] = max(rnorm(n=1,(a[floor((i-1)/100)+1]+b[floor((i-1)/100)+1])/(b[floor((i-1)/100)+1]*exp(-a[floor((i-1)/100)+1]*(t[i]-tmax[floor((i-1)/100)+1])) + a[floor((i-1)/100)+1] * exp(b[floor((i-1)/100)+1]*(t[i]-tmax[floor((i-1)/100)+1]))), 0.01),0)

plot(vl)

model_multiple_vl <- stan_model("multiple_vl.stan")

initfun <- function(...) {
  list(a_bar=0.1, a_sigma = 0.1, a=c(0.1,0.1,0.1), b=0.05, tmax = 21, sigma=0.0001)
}

fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), iter=10000, init = initfun, chains=1, seed=1)

print(fit)

