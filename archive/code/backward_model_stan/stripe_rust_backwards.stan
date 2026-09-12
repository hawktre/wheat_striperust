// Stan model for Mixed Membershp approach
functions {
  /**
   * Calculates neighborhood dispersal pressure restricted to a specific hidden component.
   * source_in_comp is an indicator vector (1 if node belongs to component k, 0 otherwise).
   */
  vector calc_dispersal_backward(matrix dist,
                                 matrix wind,
                                 vector source_strength,
                                 real kappa,
                                 vector source_in_comp) {
    int N = size(source_strength);
    matrix[N, N] kernel;
    vector[N] dispersal_pressure;
    
    // Element-wise mask ensures infection only emanates from source group k
    vector[N] masked_source = source_strength .* source_in_comp;
    
    for (i in 1 : N) {
      for (j in 1 : N) {
        if (i == j) {
          kernel[i, j] = 0;
        } else {
          kernel[i, j] = pow(1 + dist[i, j], -kappa);
        }
      }
    }
    
    dispersal_pressure = (wind .* kernel) * masked_source;
    return dispersal_pressure;
  }
}
data {
  int<lower=1> N; // Number of observed nodes
  int<lower=1> K; // Number of latent source combinations (ncol of combos)
  
  array[N] real<lower=0> y; // Infected plant proportion this week
  vector[N] y_lag;          // Infected last week
  
  matrix[N, N] dist_mat;
  matrix[N, N] wind_mat;
  
  // Logical indicator matrix mapping group memberships: N rows x K combinations
  matrix[N, K] source_in_component;
  
}
parameters {
  real beta; 
  real<lower=0> delta; 
  real<lower=0> gamma; 
  real<lower=0> phi; 
  real<lower=0> kappa; 
  real<lower=0, upper=1> alpha;
  simplex[K] pi_vec; // Mixture weights for the K source configurations
}
model {
  // Priors
  beta  ~ normal(0, 2.5);
  delta ~ normal(0, 2.5);
  gamma ~ normal(0, 2.5);
  kappa ~ lognormal(1, 1);
  phi   ~ inv_gamma(0.1, 0.1);
  pi_vec ~ dirichlet(rep_vector(1.0, K));

  // Marginal log-likelihood calculation
  for (n in 1:N) {
    vector[K] lp_components;
    real auto_lag = (y_lag[n] * (1.0 - y_lag[n]));
    
    for (k in 1:K) {
      // 1. Calculate dispersal localized strictly to component k
      vector[N] disp_k = calc_dispersal_backward(dist_mat, wind_mat, y_lag, kappa, source_in_component[, k]);
      
      // 2. Link expected infection proportion
      real mu_nk = inv_logit(beta + delta * auto_lag + gamma * disp_k);
      mu_nk = mu_nk * (1.0 - 2e-5) + 1e-5; // Protection boundary buffer
      
      // 3. Complete density calculation for component k
      if (y[n] == 0) {
        lp_components[k] = log(pi_vec[k]) + log(alpha);
      } else {
        lp_components[k] = log(pi_vec[k]) + log1m(alpha) + beta_proportion_lpdf(y[n] | mu_nk, phi);
      }
    }
    
    // Sum out latent profiles using the log_sum_exp identity
    target += log_sum_exp(lp_components);
  }
}
generated quantities {
  // Calculate posterior responsibilities (p_mat) for each observation n and configuration k
  matrix[N, K] p_mat;
  
  for (n in 1:N) {
    vector[K] log_prob_k;
    real auto_scaled = (y_lag[n] * (1.0 - y_lag[n]));
    
    for (k in 1:K) {
      vector[N] disp_raw_k = calc_dispersal_backward(dist_mat, wind_mat, y_lag, kappa, source_in_component[, k]);
      real disp_scaled_k = disp_raw_k[n];
      real mu_nk = inv_logit(beta + delta * auto_scaled + gamma * disp_scaled_k);
      
      if (y[n] == 0) {
        log_prob_k[k] = log(pi_vec[k]) + log(alpha);
      } else {
        log_prob_k[k] = log(pi_vec[k]) + log1m(alpha) + beta_proportion_lpdf(y[n] | mu_nk, phi);
      }
    }
    // Normalize log probability profile across components to yield valid probabilities
    p_mat[n, ] = to_row_vector(softmax(log_prob_k));
  }
}
