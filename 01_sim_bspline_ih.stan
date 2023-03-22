data{
  int <lower=1> I; // number of items
  int <lower=1> P; // number of examinees
  int <lower=2> K; // number of categories for all items K_1 =...= K_I = K
  array[P, I] int <lower=1, upper=K> R; // complete response matrix 
  int <lower=1> M; // number of b-spline basis functions and knots
  array[I] matrix[K-1, M] B; // array of item-specific b-spline matrices
}

parameters{
  vector[P] theta; // vector of ability parameters
  array[I] real<lower=0,upper=pi()/2> alpha_raw; // raw discrimination parameter before reparametrization
  
  array[I] ordered[M] lambda; // category parameter
  array[I] real mu_lambda_raw; // raw b-spline coefficient location before reparametrization
  array[I] real<lower=0,upper=pi()/2>  sigma_lambda_raw; // raw b-spline coefficient scale before reparametrization

}

transformed parameters{
  array[I] real alpha; // reparametrized discrimination parameter
  array[I] real sigma_lambda; // reparametrized standard deviations for B-spline coefficients
  array[I] real mu_lambda; // reparametrized means for B-spline coefficients
  matrix[I, K-1] gamma; // threshold parameter matrix
  
  for(i in 1:I) {
    // reparametrize to Cauchy(0,3)
    alpha[i] = 3 * tan(alpha_raw[i]);
    // reparametrize to Cauchy(0,5)
    sigma_lambda[i] = 5 * tan(sigma_lambda_raw[i]);
    // reparametrize to N(0,5)
    mu_lambda[i] = 5 * mu_lambda_raw[i];
    // calculate spline from the provided b-spline basis B and coefficients 
    // lambda per item to model thresholds parameters
    gamma[i] = to_row_vector(B[i] * to_vector(lambda[i]));
  }
}

model{
  // for central normal parametrizatized
  theta ~ std_normal(); 
  mu_lambda_raw ~ std_normal();

  // sampling B-spline coefficients from normal distribution
  for(i in 1:I) {
    for(m in 1:M) {
      lambda[i, m] ~ normal(mu_lambda[i], sigma_lambda[i]);
    }
  }
  
  // sampling responses from ordered logistic distribution (multinomial with 
  // logistic differences as probability vector)
  for(p in 1:P) {
    for(i in 1:I) {
      R[p, i] ~ ordered_logistic(alpha[i] * theta[p], to_vector(gamma[i]));
    }
  }
}
