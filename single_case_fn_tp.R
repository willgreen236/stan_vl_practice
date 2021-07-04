library(rstan)

N <- 100
t <- seq(0, N-1, 1)

vl_func <- function(a, b, tmax, t, log_vlmax) return(log10(10^log_vlmax*(a+b)/(b*exp(-a*(t-tmax))+a*exp(b*(t-tmax)))))

vl <- pmax(rnorm(n=N,mean=vl_func(0.5, 0.2, 30, t, 8), sd=0.01) * rbinom(n=100,size=1,p=0.91),4)
length(which(vl==4))
plot(vl)

#single_case_vl <- stan_model("single_case_fn_tp.stan")
single_case_vl <- stan_model("single_case_fn_tp_2.stan")

initfun_simple <- function(...) {
  list(a=0.3, b=0.1, tmax = 25, log_vlmax = 5, sigma=0.1, l_fn = 0.9)
}

#expose_stan_functions("single_case_fn_tp.stan")

fit_simple <- sampling(single_case_vl, list(N=N, t=t, vl = vl), init=initfun_simple, iter=1000, chains=1, seed=1)

print(fit_simple)

post_fit <- function(stan_fit, param) return(as.numeric(get_posterior_mean(stan_fit,par=param)))

checking_df <- data.frame(
  t = seq(0,99),
  input = vl[1:100],
  fit = vl_func(post_fit(fit_simple, "a"), post_fit(fit_simple, "b"), post_fit(fit_simple, "tmax"), seq(0,99), post_fit(fit_simple, "log_vlmax"))
)

ggplot(checking_df, aes(x=t)) +
  geom_line(aes(y=fit)) +
  geom_point(aes(y=input))
