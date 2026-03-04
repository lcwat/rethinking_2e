//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;
  vector[N] d;
  vector[N] m;
  vector[N] a;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b0; // intercept
  real b1; // effect of marriage rate
  real b2; // effect of age at first marriage
  real<lower=0> sigma; // variance bound at zero
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] mu; // vector of means
  mu = b0 + b1*m + b2*a; // linear additive model
  d ~ normal(mu, sigma); // likelihood
  
  // priors
  b0 ~ normal(0, 0.2);
  b1 ~ normal(0, 0.5);
  b2 ~ normal(0, 0.5);
  sigma ~ exponential(1);
}

// generate quantities block: compute log pointwise predictive densities for 
// model comparison 
generated quantities {
  // prediction quantities, without, can't generate fitted values
  vector[N] linpred = b0 + b1*m + b2*a;
  vector[N] epred = linpred; // apply link function here if needed
  array[N] real prediction = normal_rng(epred, sigma);
  
  // likelihoods for model comparison
  vector[N] log_lik;
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(d[i] | b0 + b1*m[i] + b2*a[i], sigma);
  }
}

