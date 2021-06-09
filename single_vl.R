library(rstan)

N <- 100
t <- seq(0,99,1)
a <- 0.5
b <- 0.1
tmax <- 30
vl <- vector(length = N)
for (i in 1:N) vl[i] = max(rnorm(n=1,(a+b)/(b*exp(-a*(t[i]-tmax)) + a * exp(b*(t[i]-tmax))), 0.01),0)

plot(vl)

model_single_vl <- stan_model("single_vl.stan")

initfun <- function(...) {
  list(a=0.3, b=0.05, tmax = 25, sigma=0.01)
}

fit <- sampling(model_single_vl, list(N=N, t=t, vl = vl), iter=1000, init = initfun, chains=1, seed=1)

print(fit)

