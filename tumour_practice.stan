data{
  int<lower=1> nStudy; // number of studies
  int N; // number of samples per study
  int<lower=0, upper=N> X[nStudy]; // number successes
}

parameters{
  ordered[2] alpha;
}

transformed parameters{
  real<lower=0, upper=1> theta[2];
  matrix[N,2] lp;
  
  for(i in 1:2){
    theta[i] = inv_logit(alpha[i]);
  }
  
  for(n in 1:N){
    for(s in 1:2){
      lp[n,s] = log(0.5) + binomial_logit_lpmf(X[N]|N, alpha[s]);
    }
  }
  
}

model{
  for(n in 1:N){
    target += log_sum_exp(lp[n]);
  }
}