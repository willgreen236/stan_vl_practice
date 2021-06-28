functions {
  real vlfunc(real a1, real b1, real t1, real tmax1, real log_vlmax1){
    return(log10(10^log_vlmax1*(a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1)))));
  } 
  real example_lpdf(real aX, real a1, real b1, real t1, real tmax1, real log_vlmax1, real aSigma, real l_fn1){
    real output;
    
    if(aX > 0)  output = log(1-l_fn1) + normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
    if(aX == 0) output = log(l_fn1); //+ normal_lpdf(0 | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
    return(output);
  }
}

data { 
  int N;             // number of time points
  int<lower = 0> K;  // number of individuals

  real t[N];         // time points of observations 
  //int individual[N]; // identity of individuals

  matrix[N,K] vl;        // viral load reads
}

parameters{
  real<lower=0> a[K];
  real<lower=0> b[K];
  real<lower=0> tmax[K];
  real<lower=0> log_vlmax[K];
  real<lower=0> sigma;
  real<lower=0> a_bar;
  real<lower=0> a_sigma;
  real<lower=0> b_bar;
  real<lower=0> b_sigma;
  real<lower=0> t_bar;
  real<lower=0> t_sigma;
  real<lower=0> log_vlmax_bar;
  real<lower=0> log_vlmax_sigma;
  
  //real<lower=0> l_fp; // false positive probability
  real<lower=0> l_fn; // false negative probability
}

model{
  for(i in 1:N){
    for(j in 1:K){
      vl[i,j] ~ example(a[j], b[j], t[i], tmax[j], log_vlmax[j], sigma, l_fn);
    }
  }
  // priors
  a ~ normal(a_bar,a_sigma);
  b ~ normal(b_bar,b_sigma);
  tmax ~ normal(t_bar,t_sigma);
  log_vlmax ~ normal(log_vlmax_bar, log_vlmax_sigma);
  sigma ~ uniform(0,1);
  
  // hyper-priors
  a_bar ~ normal(1,1);
  a_sigma ~ normal(1,1);
  b_bar ~ normal(1,1);
  b_sigma ~ normal(1,1);
  t_bar ~ normal(30,10);
  t_sigma ~ normal(10,3);
  log_vlmax_bar ~ normal(10,4);
  log_vlmax_sigma ~ normal(2,2);
  
}