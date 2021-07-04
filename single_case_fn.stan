functions {
  real vlfunc(real a1, real b1, real t1, real tmax1, real log_vlmax1){
    return(log10(10^log_vlmax1*(a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1)))));
  } 
  real example_lpdf(real aX, real a1, real b1, real t1, real tmax1, real log_vlmax1, real aSigma, real l_fn1){
    real output;
    if(aX > 4)  output = log(1-l_fn1) + normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
    if(aX == 4) output = log(  l_fn1 +  (1-l_fn1) * exp(normal_lcdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma)));
    
    //print(output);
    //if(aX > 4)  output = log(1-l_fn1) + normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
    //if(aX == 4) output = log(  l_fn1 +  (1-l_fn1) * exp(normal_lcdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma)/normal_lcdf(100 | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma)));
    
    return(output);
  }
}

data { 
  int N;             // number of time points
  real t[N];         // time points of observations 
  real vl[N];        // viral load reads
}

parameters{
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> tmax;
  real<lower=0> log_vlmax;
  real<lower=0> sigma;

  real<lower=0, upper=1> l_fn; // false negative probability
}

model{
  for(i in 1:N){
    vl[i] ~ example(a, b, t[i], tmax, log_vlmax, sigma, l_fn);
  }
  
  // priors
  a ~ normal(1,1);
  b ~ normal(1,1);
  tmax ~ normal(30,10);
  log_vlmax ~ normal(10, 4);
  sigma ~ uniform(0,1);
  
}