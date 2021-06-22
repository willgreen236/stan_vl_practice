functions {
  real vlfunc(real a1, real b1, real t1, real t_onset, real tmax1){
    return((a1 + b1)/(b1 * exp(-a1*(t1+t_onset-tmax1)) + a1 * exp(b1*(t1+t_onset-tmax1))));
  } 
}

data { 
  int N;             // number of level-1 observations
  int<lower = 0> K;  // number of level-2 clusters (individuals)

  real t[N];         // time points of observations 
  int individual[N]; // identity of individuals

  real vl[N];        // viral load reads
}

parameters{
  real<lower=0> a[K];
  real<lower=0> b[K];
  real<lower=0> tmax[K];
  real<lower=0> t_onset[K];
  real<lower=0> sigma;
  real<lower=0> a_bar;
  real<lower=0> a_sigma;
  real<lower=0> b_bar;
  real<lower=0> b_sigma;
  real<lower=0> t_bar;
  real<lower=0> t_sigma;
  real<lower=0> t_onset_bar;
  real<lower=0> t_onset_sigma_bar;
  
}

model{
  for(i in 1:N){
    vl[i] ~ normal(vlfunc(a[individual[i]], b[individual[i]], t[i], t_onset[individual[i]], tmax[individual[i]]), sigma);
  }
  // priors
  a ~ normal(a_bar,a_sigma);
  b ~ normal(b_bar,b_sigma);
  tmax ~ normal(t_bar,t_sigma);
  t_onset ~ normal(t_onset_bar, t_onset_sigma_bar);
  sigma ~ uniform(0,1);
  
  // hyper-priors
  a_bar ~ normal(1,1);
  a_sigma ~ normal(1,1);
  b_bar ~ normal(1,1);
  b_sigma ~ normal(1,1);
  t_bar ~ normal(30,10);
  t_sigma ~ normal(10,3);
  t_onset_bar ~ normal(20,20);
  t_onset_sigma_bar ~ normal(0,10);
  
}