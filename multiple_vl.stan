functions {
  real vlfunc(real a1, real b1, real t1, real tmax1){
    return((a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1))));
  }
}

data { 
  int N;
  int K;
  real t[N];
  real vl[N];
  int individual[N];
}

parameters{
  real<lower=0> a[K];
  real<lower=0> b;
  real<lower=0> tmax;
  real<lower=0> sigma;
  real<lower=0> a_bar;
  real<lower=0> a_sigma;
}

model{
  for(i in 1:N){
    vl[i] ~ normal(vlfunc(a[individual[i]], b, t[i], tmax), sigma);
  }
  // priors
  a ~ normal(a_bar,a_sigma);
  b ~ uniform(0,1);
  tmax ~ uniform(20,40);
  sigma ~ uniform(0,1);
  
  // hyper-priors
  a_bar ~ normal(1,1);
  a_sigma ~ normal(1,1);
}