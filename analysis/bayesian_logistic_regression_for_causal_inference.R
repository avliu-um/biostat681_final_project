library(mvtnorm)
library(invgamma)
library(pgdraw)

bayesian_causal_inference_logistic_regression = function(X, y, treatment_indicator, 
                                                         a_0 = 1, b_0 = 1, a_1 = 1, b_1 = 1, 
                                                         total_steps = 10000) {
  control_indicator = (!treatment_indicator)
  X = cbind(1, X)
  Xt = t(X)
  n = dim(X)[1]
  p = dim(X)[2]
  n_1 = sum(treatment_indicator)
  n_0 = sum(control_indicator)
  
  # initialization
  # sigma_1_sq = rinvgamma(n = 1, shape = a_1, scale = b_1)
  # sigma_0_sq = rinvgamma(n = 1, shape = a_0, scale = b_0)
  # beta_1 = rnorm(p, mean = 0, sd = sqrt(sigma_1_sq))
  # beta_0 = rnorm(p, mean = 0, sd = sqrt(sigma_0_sq))
  # p_1 = 1 / (1 + exp(- X %*% beta_1))
  # p_0 = 1 / (1 + exp(- X %*% beta_0))
  
  sigma_1_sq = 1
  sigma_0_sq = 1
  beta_1 = rep(0, p)
  beta_0 = rep(0, p)
  p_1 = 1 / (1 + exp(- X %*% beta_1))
  p_0 = 1 / (1 + exp(- X %*% beta_0))
  
  y_1_complete = rep(0, n)
  y_0_complete = rep(0, n)
  y_1_complete[treatment_indicator] = y[treatment_indicator]
  y_0_complete[control_indicator] = y[control_indicator]
  y_1_complete[control_indicator] = rbinom(n = n_0, size = 1, prob = p_1[control_indicator])
  y_0_complete[treatment_indicator] = rbinom(n = n_1, size = 1, prob = p_0[treatment_indicator])
  
  g_1 = draw_g(n, X %*% beta_1)
  g_0 = draw_g(n, X %*% beta_0)
  
  # setup
  tau_p_samples = rep(0, total_steps)
  sigma_samples = matrix(0, nrow = total_steps, ncol = 2)
  beta_1_samples = matrix(0, nrow = total_steps, ncol = p)
  beta_0_samples = matrix(0, nrow = total_steps, ncol = p)
  
  # gibbs sampler
  pb = txtProgressBar(style = 3)
  for (s in 1:total_steps) {
    kappa_1 = y_1_complete - 0.5
    kappa_0 = y_0_complete - 0.5
    
    Sigma_1 = Xt %*% diag(g_1) %*% X
    diag(Sigma_1) = diag(Sigma_1) + 1 / sigma_1_sq
    Sigma_1 = solve(Sigma_1)
    mu_1 = Sigma_1 %*% Xt %*% kappa_1
    
    Sigma_0 = Xt %*% diag(g_0) %*% X
    diag(Sigma_0) = diag(Sigma_0) + 1 / sigma_0_sq
    Sigma_0 = solve(Sigma_0)
    mu_0 = Sigma_0 %*% Xt %*% kappa_0
    
    # draw g
    g_1 = draw_g(n, X %*% beta_1)
    g_0 = draw_g(n, X %*% beta_0)
    
    # draw beta
    beta_1 = t(rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1))
    beta_0 = t(rmvnorm(n = 1, mean = mu_0, sigma = Sigma_0))
    
    # draw sigma
    sigma_1_sq = rinvgamma(n = 1, shape = a_1 + p / 2, scale = b_1 + sum(beta_1^2) / 2)
    sigma_0_sq = rinvgamma(n = 1, shape = a_0 + p / 2, scale = b_0 + sum(beta_0^2) / 2)
    
    # impute missing data
    p_1 = 1 / (1 + exp(- X %*% beta_1))
    p_0 = 1 / (1 + exp(- X %*% beta_0))
    y_1_complete[control_indicator] = rbinom(n = n_0, size = 1, prob = p_1[control_indicator])
    y_0_complete[treatment_indicator] = rbinom(n = n_1, size = 1, prob = p_0[treatment_indicator])
    
    # store parameters
    tau_p_samples[s] = mean(y_1_complete) - mean(y_0_complete)
    sigma_samples[s, ] = c(sigma_0_sq, sigma_1_sq)
    beta_1_samples[s, ] = beta_1
    beta_0_samples[s, ] = beta_0
    
    setTxtProgressBar(pb, s / total_steps)
  }
  close(pb)
  return(list(tau_p_samples = tau_p_samples, 
              sigma_samples = sigma_samples,
              beta_1_samples = beta_1_samples, 
              beta_0_samples = beta_0_samples))
}

draw_g = function(n, eta) {
  g_sample = rep(0, n)
  for (i in 1:n) {
    g_sample[i] = pgdraw::pgdraw(1, eta[i])
  }
  return(g_sample)
}
