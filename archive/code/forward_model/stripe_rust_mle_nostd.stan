// Stan model for Maximum Likelihood Estimation of Hop Model Parameters
// Scaled via invariant global standard deviations calculated in transformed data

functions {
  /**
   * Calculates neighborhood dispersal pressure for a specific block.
   * Excludes self-infection (j != i).
   */
  vector calc_dispersal(matrix dist,
                        matrix wind,
                        vector source_strength,
                        real kappa) {
    int N = size(source_strength);
    matrix[N, N] kernel;
    vector[N] dispersal_pressure;
    
    // Build the power-law dispersal kernel
    for (i in 1 : N) {
      for (j in 1 : N) {
        if (i == j) {
          kernel[i, j] = 0;
        } else {
          kernel[i, j] = pow(1 + dist[i, j], -kappa);
        }
      }
    }
    
    // Weight by wind and calculate arriving pressure
    dispersal_pressure = (wind .* kernel) * source_strength;
    
    return dispersal_pressure;
  }
}

data {
  int<lower=1> N; // Number of observations
  
  // Outcomes
  array[N] real<lower=0> y; // Infected plants (Outcome)
  
  // Lagged predictors
  vector[N] y_lag; // Infected last month
  
  // Wind and Distance Matrices
  matrix[N, N] dist_mat;
  matrix[N, N] wind_mat;

  // Prior Hyperparameters
  real beta_mu;    real beta_sigma;
  real delta_mu;   real delta_sigma;
  real gamma_mu;   real gamma_sigma;
  real kappa_mu;   real kappa_sigma;
  real phi_alpha;  real phi_beta;
}

transformed data {
  // 1. Compute static scale for autoinfection
  vector[N] auto_raw = y_lag .* (1 - y_lag);
}

parameters {
  real beta;                              // Global Intercept
  real<lower=0> delta;                    // Auto-infection magnitude
  real<lower=0> gamma;                    // Dispersal magnitude
  real<lower=0> phi;                      // Variance parameter
  real<lower=0> kappa;                    // Dispersal kernel decay
  real<lower=0, upper=1> alpha;          // Zero-inflation term
}

transformed parameters {
  vector[N] logit_mu;
  vector[N] mu;
  
  // Calculate dynamic dispersal based on current iteration's kappa
  vector[N] dispersal_raw = calc_dispersal(dist_mat, wind_mat, y_lag, kappa);
  
  // Isotropic linear predictor
  logit_mu = beta + delta * auto_raw + gamma * dispersal_raw;
  mu = inv_logit(logit_mu);
}

model {
  // Priors
  beta ~ normal(beta_mu, beta_sigma);
  delta ~ normal(delta_mu, delta_sigma);
  gamma ~ normal(gamma_mu, gamma_sigma);
  kappa ~ normal(kappa_mu, kappa_sigma);
  phi ~ inv_gamma(phi_alpha, phi_beta);

  // Vectorized Likelihood Evaluation
  for (n in 1:N) {
    if (y[n] == 0) {
      target += log(alpha);
    } else {
      target += log1m(alpha) + beta_proportion_lpdf(y[n] | mu[n], phi);
    }
  }
}

generated quantities {
  vector[N] deviance_residual;
  real deviance = 0;

  for (n in 1:N) {
    real log_lik_saturated;
    real log_lik_model;
    real d_sq; 

    if (y[n] == 0) {
      log_lik_saturated = 0.0; 
      log_lik_model = log(alpha);
    } else {
      log_lik_saturated = log1m(alpha) + beta_proportion_lpdf(y[n] | y[n], phi);
      log_lik_model = log1m(alpha) + beta_proportion_lpdf(y[n] | mu[n], phi);
    }

    d_sq = 2.0 * (log_lik_saturated - log_lik_model);
    if (d_sq < 0) d_sq = 0.0;

    if (y[n] >= mu[n]) {
      deviance_residual[n] = sqrt(d_sq);
    } else {
      deviance_residual[n] = -sqrt(d_sq);
    }

    deviance += d_sq;
  }
}
