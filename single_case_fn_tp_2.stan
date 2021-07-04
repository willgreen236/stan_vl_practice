functions {
  real vlfunc(real a1, real b1, real t1, real tmax1, real log_vlmax1){
    //print(a1);
    return(log10(10^log_vlmax1*(a1 + b1)/(b1 * exp(-a1*(t1-tmax1)) + a1 * exp(b1*(t1-tmax1)))));
  } 
  //real example_lpdf(real aX, real a1, real b1, real t1, real tmax1, real log_vlmax1, real aSigma, real l_fn1){
  //  real output;
  //  if(aX > 4)  output = log(1-l_fn1) + normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
  //  if(aX == 4) output = log(  l_fn1 +  (1-l_fn1) * exp(normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma)));
    
    //if(aX == 4) output = log_sum_exp(log(1-exp(-l_fn1))+normal_lcdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma), -l_fn1);
    //if(aX > 4) output = log(1-exp(-l_fn1))+normal_lpdf(aX | vlfunc(a1, b1, t1, tmax1, log_vlmax1), aSigma);
    
  //return(output);
  //}
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

//transformed parameters{
// 
//  matrix[N, 2] lp;
  
//  for(i in 1:N){
//    for(s in 1:2){
//      if(vl[i] > 4){
//        lp[i,1] = log(1-l_fn) + normal_lpdf(vl[i] | vlfunc(a, b, t[i], tmax, log_vlmax), sigma);
//        lp[i,2] = 0;
//      }  
//      if(vl[i] == 4){
//        lp[i, 1] = log(1-l_fn) + normal_lcdf(vl[i] | vlfunc(a, b, t[i], tmax, log_vlmax), sigma);
//        lp[i, 2] = log(l_fn);
//      }    
//    }
//  }
//}

transformed parameters{
  real ll_total;
  vector[N] log_lik;
  //real fp;
  real fn;
  
  {
    real l_fnm;
    real cum_dens_lik;
    fn = exp(-l_fn);
    l_fnm = log(1.0-fn);
    
    for(n in 1:N){
      //if(vl[n] == 4.0) log_lik[n] =  log_sum_exp(l_fnm+normal_lcdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma),-l_fn);
      //else             log_lik[n] =              l_fnm+normal_lpdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma);
      
      //if(vl[n] > 4.0)       log_lik[n] = log(1-l_fn) + normal_lpdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma);
      //else if(vl[n] == 4.0) log_lik[n] = log(l_fn)
      
      if(vl[n] > 4.0)       log_lik[n] = log(1-l_fn) + normal_lpdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma);
      else if(vl[n] == 4.0) log_lik[n] = log(l_fn + (1-l_fn) * exp(normal_lcdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma)));
      
      //if(vl[n] > 4.0) log_lik[n] = log(1-l_fn) + normal_lpdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma);
      //else if(vl[n] == 4.0){
       // cum_dens_lik = normal_lcdf(vl[n] | vlfunc(a, b, t[n], tmax, log_vlmax), sigma);
        //print(cum_dens_lik);
        //if (is_nan(cum_dens_lik)) log_lik[n] = log(l_fn + (1-l_fn) * exp(cum_dens_lik));
        //else log_lik[n] = log(l_fn);
      //} 
      //print(log_lik[n]);
      //print(normal_cdf(4 | 6, 0.1));
    }
    ll_total = sum(log_lik);
    //print(ll_total);
  }
}

model{
  //for(i in 1:N){
    target += ll_total;
    //vl[i] ~ example(a, b, t[i], tmax, log_vlmax, sigma, l_fn);
  //}
  
  // priors
  a ~ normal(1,1);
  b ~ normal(1,1);
  tmax ~ normal(30,4);
  log_vlmax ~ normal(7, 2);
  sigma ~ uniform(0,1);
}