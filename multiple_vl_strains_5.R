library(rstan)

N <- 1000
K <- 10
S <- 2
t <- rep(seq(0,99,1),K)
set.seed(1)
a <- c(rnorm(n=K, mean=0.55, sd=0.05))
b <- rnorm(n=K, mean=0.1, sd=0.02)
t_onset <- rnorm(n=K, mean=20, sd=5)
tmax <- rnorm(n=K, mean=30, sd=4)
vlmax <- rnorm(n=K, mean=8, sd=1)
strain <- c(rep(1, N/2), rep(2, N/2))

individual <- vector(length=N)
for (i in 1:length(individual)) individual[i] = floor((i-1)/(N/K) + 1)
  
vl <- vector(length = N)
for (i in 1:N){
  index <- floor((i-1)/(N/K))+1
  vl[i] = max(rnorm(n=1,log10(10^vlmax[index]*(a[index]+b[index])/(b[index]*exp(-a[index]*(t[i]-tmax[index])) + a[index] * exp(b[index]*(t[i]-tmax[index])))), 0.01),0)
} 

plot(vl[1:100])  

model_multiple_vl <- stan_model("multiple_vl_strains_5.stan")

initfun <- function(...) {
  list(a_bar=0.3, a_sigma = 0.1, a=rep(0.3,K), b_bar=0.05, b_sigma = 0.01, b=rep(0.05, K), tmax = rep(21,K), t_bar = 25, log_vlmax_bar = 7, log_vlmax_sigma = 1, log_vlmax = rep(8,K), t_sigma = 4, sigma=0.0001)
}

fit <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), init=initfun, iter=1000, chains=1, seed=1)

post_fit <- function(stan_fit, param) return(as.numeric(get_posterior_mean(stan_fit,par=param)))

vl_func <- function(a, b, tmax, t, log_vlmax) return(log10(10^log_vlmax*(a+b)/(b*exp(-a*(t-tmax))+a*exp(b*(t-tmax)))))

checking_df <- data.frame(
  t = seq(0,99),
  input = vl[1:100],
  fit = vl_func(post_fit(fit, "a[1]"), post_fit(fit, "b[1]"), post_fit(fit, "tmax[1]"), seq(0,99), post_fit(fit, "log_vlmax[1]"))
)

ggplot(checking_df, aes(x=t)) +
  geom_line(aes(y=fit)) +
  geom_point(aes(y=input))
