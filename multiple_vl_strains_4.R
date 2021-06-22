library(rstan)

N <- 3000
K <- 30
S <- 2

vector(length = N)

set.seed(1)
a <- c(rnorm(n=K, mean=0.55, sd=0.05))
b <- rnorm(n=K, mean=0.1, sd=0.02)
t_onset <- round(rnorm(n=K, mean=20, sd=5),0)
tmax <- round(rnorm(n=K, mean=30, sd=4),0)
log_vlmax <- rnorm(K, 7, 1)
strain <- c(rep(1, N/2), rep(2, N/2))

t <- vector(length = N)
for (i in 1:K) t[(i*100-99):(i*100)] = seq(-t_onset[i], 100-t_onset[i]-1,1)

individual <- vector(length=N)
for (i in 1:length(individual)) individual[i] = floor((i-1)/(N/K) + 1)
  
vl <- vector(length = N)
for (i in 1:N){
  index <- floor((i-1)/(N/K))+1
  vl[i] = max(log10(rnorm(n=1,10^log_vlmax[index]*(a[index]+b[index])/(b[index]*exp(-a[index]*(t[i]+t_onset[index]-tmax[index])) + a[index] * exp(b[index]*(t[i]+t_onset[index]-tmax[index]))), 0.01)),0)
}

which(is.na(vl)==TRUE)

plot(vl[1:100])
  
model_multiple_vl <- stan_model("multiple_vl_strains_4.stan")

#initfun <- function(...) {
#  list(a_bar=0.1, a_sigma = 0.1, a=rep(0.1,K), b_bar=0.05, b_sigma = 0.01, b=rep(0.05, K), tmax = rep(21,K), t_bar = 25, t_sigma = 4, log_vlmax = rep(6,K), log_vlmax_bar = 6, log_vlmax_sigma = 1, sigma=0.0001)
#}

fit3 <- sampling(model_multiple_vl, list(N=N, K=K, t=t, vl = vl, individual = individual), iter=2000, chains=1, seed=1)

#print(fit)
#print(fit2)
print(fit3)

vl_func <- function(a, b, tmax, t_onset, t, log_vlmax){
  return(log10(10^log_vlmax*(a+b)/(b*exp(-a*(t+t_onset-tmax))+a*exp(b*(t+t_onset-tmax)))))
}

tmax - t_onset

vl_test <- data.frame(t = seq(1,100), 
                      vl = vl_func(a=0.42, b=2.27, tmax = 3.6, t_onset=1.38, log_vlmax = 6.34, t = t[1:100]),
                      vl_init = vl[1:100])

ggplot(vl_test, aes(x=t)) +
  geom_line(aes(y=vl)) +
  geom_point(aes(y=vl_init))
