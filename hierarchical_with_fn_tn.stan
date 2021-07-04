functions {
  real vlfunc(real a1, real b1, real t1, real tmax1, real log_vlmax1){
    return(log10(10^log_vlmax1*(a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1)))));
  } 
}

data { 
  int N;             // number of data points
  int<lower = 0> K;  // number of individuals

  real t[N];         // time points of observations 
  int individual[N]; // identity of individuals
  vector[N] vl;        // viral load reads
}

parameters{
  //non-hierarchical parmaeters
  //real<lower=0> a;
  //real<lower=0> b;
  //real<lower=0> tmax;
  //real<lower=0> log_vlmax;
  
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
  real<lower=0, upper=1> l_fn; // false negative probability
}

transformed parameters{
  real ll_total;
  vector[N] log_lik;
  real fn;
  {
    real l_fnm;
    fn = exp(-l_fn);
    l_fnm = log(1.0-fn);
    
    //for(n in 1:N){
    //  if(vl[n] > 4.0)       log_lik[n] = log(1-l_fn) + normal_lpdf(vl[individual[n]] | vlfunc(a[individual[n]], b[individual[n]], t[individual[n]], tmax[individual[n]], log_vlmax[individual[n]]), sigma);
    //  else if(vl[n] == 4.0) log_lik[n] = log(l_fn + (1-l_fn) * exp(normal_lcdf(vl[individual[n]] | vlfunc(a[individual[n]], b[individual[n]], t[individual[n]], tmax[individual[n]], log_vlmax[individual[n]]), sigma)));
    //}
    
    //non-hierarchical likelihood
    for(n in 1:N){
      if(vl[n] > 4.0)       log_lik[n] = log(1-l_fn) + normal_lpdf(vl[n] | vlfunc(a[individual[n]], b[individual[n]], t[n], tmax[individual[n]], log_vlmax[individual[n]]), sigma);
      else if(vl[n] == 4.0) log_lik[n] = log(l_fn + (1-l_fn) * exp(normal_lcdf(vl[n] | vlfunc(a[individual[n]], b[individual[n]], t[n], tmax[individual[n]], log_vlmax[individual[n]]), sigma)));
    }
    
    ll_total = sum(log_lik);
  }
}

model{
  target += ll_total;
  
  //non-hierarchical priors
  //a ~ normal(1, 1);
  //b ~ normal(1, 1);
  //tmax ~ normal(30,10);
  //log_vlmax ~ normal(10, 4);
  sigma ~ uniform(0, 1);
  
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