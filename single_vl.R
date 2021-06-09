library(rstan)

N <- 300
t <- rep(seq(0,99,1),3)
a <- c(0.6,0.5,0.4)
b <- c(0.08,0.1, 0.12)
tmax <- c(29,30,32)

vl <- vector(length = N)
for (i in 1:N) vl[i] = max(rnorm(n=1,(a[floor((i-1)/100)+1]+b[floor((i-1)/100)+1])/(b[floor((i-1)/100)+1]*exp(-a[floor((i-1)/100)+1]*(t[i]-tmax[floor((i-1)/100)+1])) + a[floor((i-1)/100)+1] * exp(b[floor((i-1)/100)+1]*(t[i]-tmax[floor((i-1)/100)+1]))), 0.01),0)

plot(vl)

model_single_vl <- stan_model("single_vl.stan")

initfun <- function(...) {
  list(a=0.1, b=0.05, tmax = 21, sigma=0.0001)
}

fit <- sampling(model_single_vl, list(N=N, t=t, vl = vl), iter=1000, init = initfun, chains=4, cores=4, seed=1)

print(fit)

