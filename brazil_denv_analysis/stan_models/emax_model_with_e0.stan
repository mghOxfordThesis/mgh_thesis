// data variables 
data {
  int<lower=0> N;               // num measurements
  real response[N];     // response vector 
  real<lower=0.0> exposure[N];     // exposure vector
}

// the parameters
parameters {
  real<lower=0.0> Emax;
  real<lower=0.0> E50;
  real<lower=-1.0> E0;
  real<lower=0.0> gamma;
  real<lower=0.0> sigma;
}

// transformed parameters 
transformed parameters {
  real responseHat[N];
  for (ii in 1:N) {
    responseHat[ii] = E0 + (Emax*exposure[ii]^gamma)/(E50^gamma+exposure[ii]^gamma);
  }
}

// model
model {
  sigma ~ normal(0, 3);
  gamma ~ normal(10, 3);
  E50 ~ normal(0.25, 0.1);
  E0 ~ normal(0.3, 0.2);
  Emax ~ normal(2.5, 0.5);
  
  response ~ normal(responseHat, sigma);
}
