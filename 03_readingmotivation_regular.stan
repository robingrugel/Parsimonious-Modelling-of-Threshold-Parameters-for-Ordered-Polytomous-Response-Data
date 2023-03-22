functions {
  // define function for a crude estimator for the skewness of a sample by Joanes (1998)
  real skew(vector x) {
    int n = num_elements(x);
    vector[n] xdev = x - mean(x);
    real r = (sum(pow(xdev, 3.0))/n) / (pow(sum(pow(xdev, 2.0))/n , 1.5)); // with simple factor 1/n
    return(pow((1 - 1.0/n), 1.5) * r);
  }

  // define function for a crude estimator for the kurtosis of a sample by Joanes (1998)
  real kurt(vector x) {
    int n = num_elements(x);
    vector[n] xdev = x - mean(x);
    real r = (sum(pow(xdev, 4.0))/n) / (pow(sum(pow(xdev, 2.0))/n , 2)); // with simple factor 1/n
    return(pow((1 - 1.0/n), 2) * r - 3);
  }
}

data{
  int <lower=1> I; // number of items
  int <lower=1> P; // number examinees
  int <lower=2> K; // number of categories for all items K_1 =...= K_I = K
  array[2*P, I] int <lower=1, upper=K> R; // complete response matrix for both points in time
  array[P] int <lower=1, upper=2> group; // group indicator (control or treatment)
}

parameters{
  matrix[P, 2] theta; // matrix of multidimensional ability parameters
  array[I] real<lower=0,upper=pi()/2> alpha_raw; // raw discrimination parameter before reparametrization
  
  array[I] ordered[K-1] gamma; // ordered threshold parameters
  array[I] real mu_gamma_raw; // raw b-spline coefficient location before reparametrization
  array[I] real<lower=0,upper=pi()/2>  sigma_gamma_raw; // raw b-spline coefficient scale before reparametrization
  
  real mu_theta_treat_raw; // raw treatment effect before reparametrization
  real mu_all_raw; // raw overall ability intercept before reparametrization
  
  real<lower=0,upper=pi()/2> sigma_theta_control_raw; // raw standard deviation for control group before reparametrization
  real<lower=0,upper=pi()/2> sigma_theta_treat_raw; // raw standard deviation for treatment group before reparametrization
  
  cholesky_factor_corr[2] corr_control; // cholesky factor of correlation for control group
  cholesky_factor_corr[2] corr_treat; // cholesky factor of correlation for treatment group
}

transformed parameters{
  vector[2] sigma_control; // vector of standard deviations for control group
  matrix[2, 2] L_control; // L matrix for multi_normal_cholesky() for control group
  vector[2] sigma_treat; // vector of standard deviations for treatment group
  matrix[2, 2] L_treat; // L matrix for multi_normal_cholesky() for treatment group
  
  array[I] real alpha; // reparametrized discrimination parameter
  
  real sigma_theta_treat; // reparametrized standard deviation of treatment group
  real sigma_theta_control; // reparametrized standard deviation of control group
  
  array[I] real sigma_gamma; // reparametrized standard deviations for threshold parameters
  array[I] real mu_gamma; // reparametrized means for threshold parameters
  real mu_theta_treat; // reparametrized treatment effect
  real mu_all; // reparametrized overall ability intercept

  for(i in 1:I) {
    // reparametrize to Cauchy(0,3)
    alpha[i] = 3 * tan(alpha_raw[i]);
    // reparametrize to Cauchy(0,5)
    sigma_gamma[i] = 5 * tan(sigma_gamma_raw[i]);
    // reparametrize to N(0,5)
    mu_gamma[i] = 5 * mu_gamma_raw[i];
  }
  
  // reparametrize to Cauchy(0,3)
  sigma_theta_treat = 3 * tan(sigma_theta_treat_raw);
  sigma_theta_control = 3 * tan(sigma_theta_control_raw);
  
  // reparamtrize to N(0,3)
  mu_all = 3 * mu_all_raw;
  // reparamtrize to N(0,3)
  mu_theta_treat = 3 * mu_theta_treat_raw;

  // fix variance for treatment diagonal matrix
  sigma_treat[1] = 1.0;
  sigma_treat[2] = sigma_theta_treat;
  // calculate L for multi_normal_cholesky for treatment group
  L_treat = diag_pre_multiply(sigma_treat, corr_treat);

  // fix variance for control diagonal matrix
  sigma_control[1] = 1.0;
  sigma_control[2] = sigma_theta_control;
  // calculate L for multi_normal_cholesky for control group
  L_control = diag_pre_multiply(sigma_control, corr_control);
}

model{
  // for central normal parametrizatized
  mu_theta_treat_raw ~ std_normal();
  mu_all_raw ~ std_normal();
  mu_gamma_raw ~ std_normal();
  
  // lkj cholesky prior for correlation
  corr_control ~ lkj_corr_cholesky(2);
  corr_treat ~ lkj_corr_cholesky(2);

  // sampling ability parameter from multivariate normal distribution
  for(p in 1:P) {
    if(group[p] == 1) {
      // sampling for control group without treatment effect
      theta[p,] ~ multi_normal_cholesky([0, mu_all]', L_control);
    } else {
      // sampling for treatment group with treatment effect mu_theta_treat
      theta[p,] ~ multi_normal_cholesky([0, mu_all + mu_theta_treat]', L_treat);
    }
  }

  // sampling threshold parameter from normal distribution
  for(i in 1:I) {
    for(k in 1:(K-1)) {
      gamma[i, k] ~ normal(mu_gamma[i], sigma_gamma[i]);
    }
  }
  
  // sampling responses from ordered logistic distribution (multinomial with 
  // logistic differences as probability vector)
  for(p in 1:P) {
    for(i in 1:I) {
        // responses for first time point (before treatment)
        R[p, i] ~ ordered_logistic(alpha[i] * theta[p,1], to_vector(gamma[i]));
        // responses for second time point (after treatment)
        R[P + p, i] ~ ordered_logistic(alpha[i] * theta[p,2], to_vector(gamma[i]));
      }
    }
}

generated quantities{
  array[I, 2*P] real log_lik; // log-likelihood matrix for model comparison
  
  // define quantities for posterior predictive checking:
  array[2*P, I] real Rpred; // replicated response matrix
  array[I] int <lower=0, upper=1> meanRI; // mean of item responses
  array[I] int <lower=0, upper=1> sdRI; // standard deviation of item responses
  array[I] int <lower=0, upper=1> q25RI; // 25% quantile of item responses
  array[I] int <lower=0, upper=1> q75RI; // 75% quantile of item responses
  array[I] int <lower=0, upper=1> skewRI; // skewness of item responses
  array[I] int <lower=0, upper=1> kurtRI; // kurtosis quantile of item responses
  
  matrix[2,2] Corr_treat; // correlation matrix for treatment group
  matrix[2,2] Corr_control; // calculate the correlation matrix for control group
  
  matrix[2,2] Cov_treat;  // calculate the covariance matrix for treatment group
  matrix[2,2] Cov_control; // calculate the covariance matrix for control group
  
  // calculate correlation matices
  Corr_treat = multiply_lower_tri_self_transpose(corr_treat);
	Corr_control = multiply_lower_tri_self_transpose(corr_control);
	
	 // calculate covariance matices
	Cov_treat = quad_form_diag(Corr_treat, sigma_treat);
	Cov_control = quad_form_diag(Corr_control, sigma_control);
  
  // calculate log-likelihood functions for model comparison via elpd psis using loo
  for(p in 1:P) {
    for(i in 1:I) {
      // log-likelihood for first time point
      log_lik[i, p] = ordered_logistic_lpmf(R[p, i] | alpha[i] * theta[p,1], to_vector(gamma[i]));
      // log-likelihood for second time point
      log_lik[i, P + p] = ordered_logistic_lpmf(R[P + p, i] | alpha[i] * theta[p,2], to_vector(gamma[i]));
    }
  }

  // sample replicated data for visual posterior predictive checking
  for(p in 1:P) {
    for(i in 1:I) {
      // replicated data for first time point
      Rpred[p, i] = ordered_logistic_rng(alpha[i] * theta[p,1], to_vector(gamma[i]));
      // replicated data for second time point
      Rpred[P + p, i] = ordered_logistic_rng(alpha[i] *  theta[p,2], to_vector(gamma[i]));
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
