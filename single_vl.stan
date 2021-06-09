functions {
  real vlfunc(real a1, real b1, real t1, real tmax1){
    return((a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1))));
  }
}

data { 
  int N;
  real t[N];
  real vl[N];
}

parameters{
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> tmax;
  real<lower=0> sigma;
}

model{
  for(i in 1:N){
    vl[i] ~ normal(vlfunc(a, b, t[i], tmax), sigma);
  }
  a ~ uniform(0,1);
  b ~ uniform(0,1);
  tmax ~ uniform(20,40);
  sigma ~ uniform(0,1);
}