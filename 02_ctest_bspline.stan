functions {
  // define function for a crude estimator for the skewness of a sample by Joanes (1998)
  real skew(vector x) {
    int n = num_elements(x);
    vector[n] xdev = x - mean(x);
    real r = (sum(pow(xdev, 3.0))/n) / (pow(sum(pow(xdev, 2.0))/n , 1.5));
    return(pow((1 - 1.0/n), 1.5) * r);
  }

  // define function for a crude estimator for the kurtosis of a sample by Joanes (1998)
  real kurt(vector x) {
    int n = num_elements(x);
    vector[n] xdev = x - mean(x);
    real r = (sum(pow(xdev, 4.0))/n) / (pow(sum(pow(xdev, 2.0))/n , 2));
    return(pow((1 - 1.0/n), 2) * r - 3);
  }
}

data{
  int <lower=1> I; // number of items
  int <lower=1> P; // number of examinees
  int <lower=2> K; // number of categories for all items K_1 =...= K_I = K
  array[P, I] int <lower=1, upper=K> R; // complete response matrix
  int <lower=1> no_group; // number of groups
  array[P] int <lower=1, upper=no_group> group; // group indicator for different aggregations
  int <lower=1> M;  // number of b-spline basis functions and knots
  array[I] matrix[K-1, M] B; // array of item-specific b-spline matrices
}

parameters{
  vector[P] theta; // vector of ability parameters
  array[I] real<lower=0,upper=pi()/2> alpha_raw; // raw discrimination parameter before reparametrization
  
  array[I] ordered[M] lambda; // ordered b-spline coefficients
  real mu_lambda_raw; // raw b-spline coefficient location before reparametrization
  real <lower=0,upper=pi()/2> sigma_lambda_raw; // raw b-spline coefficient scale before reparametrization
  
  array[(no_group > 1) ? (no_group - 1) : 0] real mu_theta_raw;  // raw vector of means for groups before reparametrization
  array[(no_group > 1) ? (no_group - 1) : 0] real <lower=0,upper=pi()/2> sigma_theta_raw; // raw vector standard deviation for groups before reparametrization
}

transformed parameters{
  matrix[I, K-1] gamma; // threshold parameter matrix
  
  array[I] real alpha; // reparametrized discrimination parameter
  array[(no_group > 1) ? (no_group - 1) : 0] real sigma_theta; // reparametrized ability standard deviation for all groups
  array[(no_group > 1) ? (no_group - 1) : 0] real mu_theta; // reparametrized ability mean for all groups
  real mu_lambda; // reparametrized mean for B-spline coefficients
  real sigma_lambda; // reparametrized standard deviation for B-spline coefficients
  
  // reparamtrize to N(0,5)
  mu_lambda = 5 * mu_lambda_raw;
  // reparametrize to Cauchy(0,5)
  sigma_lambda = 5 * tan(sigma_lambda_raw);
  
  for(i in 1:I) {
      // calculate spline from the provided b-spline basis B and coefficients 
  // lambda per item to model thresholds parameters
    gamma[i] = to_row_vector(B[i] * to_vector(lambda[i]));
    alpha[i] = 3 * tan(alpha_raw[i]);
  }
  
  for(g in 1:(no_group-1)) {
    // reparametrize N(0,3)
    mu_theta[g] = 3 * mu_theta_raw[g];
    // reparametrize to Cauchy(0,3)
    sigma_theta[g] = 3 * tan(sigma_theta_raw[g]);
  }
}

model{
  // for central normal parametrizatized
  mu_theta_raw ~ std_normal();
  mu_lambda_raw ~ std_normal();
  
  // sample ability parameter regarding the group indicator
  for(p in 1:P) {
      if(group[p] == 1) {
        // use first group as N(0,1) distributed reference for group comparisons
        theta[p] ~ std_normal();
      }
      if(group[p] > 1) {
        theta[p] ~ normal(mu_theta[group[p]-1], sigma_theta[group[p]-1]);
    }
  }

  // sampling B-spline coefficients from normal distribution
  for(i in 1:I) {
    for(m in 1:M) {
      lambda[i, m] ~ normal(mu_lambda, sigma_lambda);
    }
  }
  
  // sampling responses from ordered logistic distribution (multinomial with logistic differences as probability vector)
  for(p in 1:P) {
    for(i in 1:I) {
      R[p, i] ~ ordered_logistic(alpha[i] * theta[p], to_vector(gamma[i]));
    }
  }
}

generated quantities{
  array[I, P] real log_lik; // log-likelihood matrix for model comparison
  
  // define quantities for posterior predictive checking:
  array[P, I] real Rpred; // replicated response matrix
  array[I] int <lower=0, upper=1> meanRI; // mean of item responses
  array[I] int <lower=0, upper=1> sdRI; // standard deviation of item responses
  array[I] int <lower=0, upper=1> q25RI; // 5% quantile of item responses
  array[I] int <lower=0, upper=1> q75RI; // 95% quantile of item responses
  array[I] int <lower=0, upper=1> skewRI; // skewness of item responses
  array[I] int <lower=0, upper=1> kurtRI; // kurtosis quantile of item responses
  
  // calculate log-likelihood functions for model comparison via elpd psis using loo
  for(p in 1:P) {
    for(i in 1:I) {
      log_lik[i, p] = ordered_logistic_lpmf(R[p, i] | alpha[i] * theta[p], to_vector(gamma[i]));
  }
}

  // sample replicated data for visual posterior predictive checking
  for(p in 1:P) {
    for(i in 1:I) {
      Rpred[p, i] = ordered_logistic_rng(alpha[i] * theta[p], to_vector(gamma[i]));
  }
}

  // calculating posterior predictive p-values for a test statistic T
  for(i in 1:I) {
    meanRI[i] = mean(Rpred[,i]) >= mean(R[,i]);
    sdRI[i] = sd(Rpred[,i]) >= sd(R[,i]);
    q25RI[i] = quantile(Rpred[,i], 0.25) >= quantile(R[,i], 0.25);
    q75RI[i] = quantile(Rpred[,i], 0.75) >= quantile(R[,i], 0.75);
    skewRI[i] = skew(to_vector(Rpred[,i])) >= skew(to_vector(R[,i]));
    kurtRI[i] = kurt(to_vector(Rpred[,i])) >= kurt(to_vector(R[,i]));
  }
}
