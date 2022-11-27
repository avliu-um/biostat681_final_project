games = read.csv("games_200_22_clean.csv")
games = games[, -1]
games$medium = factor(games$medium)
games$time_control = relevel(factor(games$time_control), ref = "standard")
games$home_country_diff = relevel(factor(games$home_country_diff), ref = "0")
games$players_sex = relevel(factor(games$players_sex), ref = "MM")
cleaned_games = model.matrix(~ . - 1, data = games)
cleaned_games = cleaned_games[, -1]
cleaned_games = as.data.frame(cleaned_games)

A = cleaned_games$mediumonline
X = as.matrix(cleaned_games[, 2:10])
X[, c("age_diff", "rating_diff")] = scale(X[, c("age_diff", "rating_diff")])
y = cleaned_games$upset

###############################################################################
# Bayesian logistic regression
source("bayesian_logistic_regression_for_causal_inference.R")
set.seed(2022)
blr = bayesian_causal_inference_logistic_regression(X, y, A, 
                                                    a_0 = 1, b_0 = 1, a_1 = 1, b_1 = 1, 
                                                    total_steps = 10000)
hist(blr$tau_p_samples[5001:1000])
quantile(blr$tau_p_samples[5001:10000], c(0.025, 0.975))

beta_0 = colMeans(blr$beta_0_samples)
beta_1 = colMeans(blr$beta_1_samples)
hist(1 / (1 + exp(-cbind(1, X) %*% beta_0)))
hist(1 / (1 + exp(-cbind(1, X) %*% beta_1)))

################################################################################
# Bayesian linear regression from hw1
X = cbind(1, X)
n = length(A)
p = ncol(X)

treatment = which(A == 1)
control = which(A == 0)
n_1 = length(treatment)
n_0 = length(control)

set.seed(2022)
a = 1
b = 1
y_1_complete = rep(0, n)
y_1_complete[treatment] = y[treatment]
y_0_complete = rep(0, n)
y_0_complete[control] = y[control]
beta_0 = rep(0, p)
beta_1 = rep(0, p)
sigma_0_sq = 1
sigma_1_sq = 1
y_1_complete[control] = rnorm(n_0, X[treatment, ] %*% beta_1, sqrt(sigma_1_sq))
y_0_complete[treatment] = rnorm(n_1, X[control, ] %*% beta_0, sqrt(sigma_0_sq))
inv_tXX = solve(t(X) %*% X)
tau_p_hat = rep(0, 10000)
# mcmc
for (s in 1:10000) {
  mu_1 = inv_tXX %*% t(X) %*% y_1_complete
  mu_0 = inv_tXX %*% t(X) %*% y_0_complete
  s_1 = inv_tXX * sigma_1_sq
  s_0 = inv_tXX * sigma_0_sq
  rss_1 = sum((y_1_complete - X %*% beta_1)^2)
  rss_0 = sum((y_0_complete - X %*% beta_0)^2)
  
  beta_1 = t(rmvnorm(1, mu_1, s_1))
  beta_0 = t(rmvnorm(1, mu_0, s_0))
  sigma_1_sq = 1 / rgamma(1, a + n / 2, b + rss_1 / 2)
  sigma_0_sq = 1 / rgamma(1, a + n / 2, b + rss_0 / 2)
  
  y_1_complete[control] = rnorm(n_0, X[treatment, ] %*% beta_1, sqrt(sigma_1_sq))
  y_0_complete[treatment] = rnorm(n_1, X[control, ] %*% beta_0, sqrt(sigma_0_sq))
  
  tau_p_hat[s] = mean(y_1_complete) - mean(y_0_complete)
}
mean(tau_p_hat)
quantile(tau_p_hat, c(0.025, 0.975))