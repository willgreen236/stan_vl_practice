library(rstan)

N <- 100
t <- seq(0, N-1, 1)

vl_func <- function(a, b, tmax, t, log_vlmax) return(log10(10^log_vlmax*(a+b)/(b*exp(-a*(t-tmax))+a*exp(b*(t-tmax)))))

vl <- rnorm(n=N,mean=vl_func(0.5, 0.1, 30, t, 8), sd=0.01) * rbinom(n=100,size=1,p=0.7)
length(which(vl==0))
plot(vl)

single_case_vl <- stan_model("single_case_fn.stan")

initfun <- function(...) {
  list(a=0.3, b=0.05, tmax = 21, log_vlmax = 5, sigma=0.0001, l_fn = 0.99)
}

fit <- sampling(single_case_vl, list(N=N, t=t, vl = vl), init=initfun, iter=10000, chains=1, seed=1)

print(fit)
