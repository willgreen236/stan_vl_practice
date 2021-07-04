library(rstan)

K <- 30
t_points <- 100
N <- K*t_points
t <- rep(seq(0,t_points-1),K)
Z <- vector(length=N)
individual = vector(length = N)
for(i in 1:N) individual[i] = floor((i-1)/t_points+1)

set.seed(1)
a <- c(rnorm(n=K, mean=0.55, sd=0.05))
b <- rnorm(n=K, mean=0.1, sd=0.02)
t_onset <- rnorm(n=K, mean=20, sd=5)
tmax <- rnorm(n=K, mean=30, sd=4)
vlmax <- rnorm(n=K, mean=8, sd=1)

vl_func <- function(a, b, tmax, t, log_vlmax) return(log10(10^log_vlmax*(a+b)/(b*exp(-a*(t-tmax))+a*exp(b*(t-tmax)))))

for (i in 1:N){
  Z[i] = pmax(rnorm(n=1,vl_func(a[individual[i]], b[individual[i]], tmax[individual[i]], t[i], vlmax[individual[i]]), 0.01) * rbinom(n=1, size=1, prob=0.9),4)
} 

plot(Z)  

model_multiple_vl <- stan_model("hierarchical_with_fn_tn.stan")

initfun <- function(...) {
  list(a_bar=0.5, a_sigma = 0.1, a=rep(0.5,K), b_bar=0.2, b_sigma = 0.01, b=rep(0.2, K), tmax = 21, t_bar = 20, log_vlmax_bar = 5, log_vlmax_sigma = 1, log_vlmax = rep(5,K), t_sigma = 4, sigma=0.01, l_fn = 0.9)
}

initfun_nh <- function(...) {
  list(a=rep(0.35,K), b=rep(0.1,K), tmax = rep(21,K), log_vlmax = rep(5,K), sigma=0.1, l_fn = 0.9)
}


p <- proc.time()
fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, individual=individual, vl = Z), init=initfun, iter=1000, chains=1, seed=1)
(proc.time()-p)["elapsed"]

q <- proc.time()
fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, individual=individual, vl = Z), init=initfun_nh, iter=200000, chains=1, cores=1, seed=1)
(proc.time()-q)["elapsed"]

print(fit)

post_fit <- function(stan_fit, param) return(as.numeric(get_posterior_mean(stan_fit,par=param)))

get_posterior_mean(fit)[1:100,]

checking_df <- function(individual){
    data.frame(t = seq(0,99),
               input = Z[((individual-1)*t_points + 1):(individual*100)],
               fit = vl_func(post_fit(fit, paste0("a[", individual,"]")), post_fit(fit, paste0("b[", individual,"]")), post_fit(fit, paste0("tmax[", individual,"]")), seq(0,99), post_fit(fit, paste0("log_vlmax[", individual,"]")))
  )
}

ggplot(checking_df(30), aes(x=t)) +
  geom_line(aes(y=fit)) +
  geom_point(aes(y=input))
