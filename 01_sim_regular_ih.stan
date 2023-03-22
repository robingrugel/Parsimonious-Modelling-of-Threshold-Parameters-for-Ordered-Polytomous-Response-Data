data{
  int <lower=1> I; // number of items
  int <lower=1> P; // number of examinees
  int <lower=2> K;// number of categories for all items K_1 =...= K_I = K
  array[P, I] int <lower=1, upper=K> R; // complete response matrix 
}

parameters{
  vector[P] theta; // vector of ability parameters
  array[I] real<lower=0,upper=pi()/2> alpha_raw; // raw discrimination parameter before reparametrization
  
  array[I] ordered[K-1] gamma; // ordered threshold parameters
  array[I] real mu_gamma_raw; // raw b-spline coefficient location before reparametrization
  array[I] real<lower=0,upper=pi()/2>  sigma_gamma_raw; // raw b-spline coefficient scale before reparametrization
}

transformed parameters{
  array[I] real alpha; // reparametrized discrimination parameter
  array[I] real sigma_gamma; // reparametrized standard deviations for threshold parameters
  array[I] real mu_gamma; // reparametrized means for threshold parameters
  
  for(i in 1:I) {
    // reparametrize to Cauchy(0,3)
    alpha[i] = 3 * tan(alpha_raw[i]);
    // reparametrize to Cauchy(0,5)
    sigma_gamma[i] = 5 * tan(sigma_gamma_raw[i]);
    // reparametrize to N(0,5)
    mu_gamma[i] = 5 * mu_gamma_raw[i];
  }
}

model{
  // for central normal parametrizatized
  theta ~ std_normal(); 
  mu_gamma_raw ~ std_normal();
  
  // sampling threshold parameter from normal distribution
  for(i in 1:I) {
    for(k in 1:(K-1)) {
      gamma[i,k] ~ normal(mu_gamma[i], sigma_gamma[i]);
    }
  }
  
  // sampling responses from ordered logistic distribution (multinomial with 
  // logistic differences as probability vector)
  for(p in 1:P) {
    for(i in 1:I) {
      R[p, i] ~ ordered_logistic(alpha[i] * theta[p], gamma[i]);
    }
  }
}
