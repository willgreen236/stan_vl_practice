library(rstan)

set.seed(1)
nStudy <- 10
tumour_status <- rbinom(n=nStudy, size = 1, p=0.2)
score <- vector(length = nStudy)

for(i in 1:nStudy) {
  if(tumour_status[i] == 0) score[i] = rbinom(n=1, size=20, p=0.3)
  if(tumour_status[i] == 1) score[i] = rbinom(n=1, size=20, p=0.7)
}

model_tumour <- stan_model("tumour_practice.stan")

#initfun <- function(...) {
#  list(a_bar=0.1, a_sigma = 0.1, a=c(0.1,0.1,0.1), b_bar=0.05, b_sigma = 0.01, b=c(0.05, 0.05, 0.05), tmax = c(20,21,22), t_bar = 25, t_sigma = 4, sigma=0.0001)
#}

fit <- sampling(model_tumour, list(nStudy = nStudy, N = N, X=score), iter=1000, chains=1, seed=1)

print(fit)

