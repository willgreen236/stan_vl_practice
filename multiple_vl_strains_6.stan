functions {
  real vlfunc(real a1, real b1, real t1, real tmax1, real log_vlmax1){
    return(log10(10^log_vlmax1*(a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1)))));
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
  
//  real<lower=0> l_fp; // false positive probability
//  real<lower=0> l_fn; // false negative probability
}

//transformed parameters{
//   real l_fpm;
//   real l_fnm;
//   real fp;
//   real fn;
   
//   fp = exp(-l_fp);
//   fn = exp(-l_fn);
//   l_fpm = log(1.0-fp);
//   l_fnm = log(1.0-fn);
//}

model{
  for(i in 1:N){
    vl[i] ~ normal(vlfunc(a[individual[i]], b[individual[i]], t[i], tmax[individual[i]], log_vlmax[individual[i]]), sigma);
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