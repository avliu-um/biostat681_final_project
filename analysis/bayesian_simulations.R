setwd("~/Library/CloudStorage/OneDrive-Umich/Assignments/BIOSTAT 681/project/biostat681_final_project/analysis")
p = 5
n = 500
beta_0 = runif(p, -0.5, 1.5)
beta_1 = runif(p, -1.5, 0.5)
A = sample(c(0, 1), n, replace = T)
n_0 = sum(A == 0)
n_1 = sum(A == 1)
X = matrix(0, nrow = n, ncol = p - 1)
X[A == 0, ] = rnorm(n_0 * (p - 1), mean = -0.2)
X[A == 1, ] = rnorm(n_1 * (p - 1), mean = 0.2)
X_augmented = cbind(1, X)

p_0 = 1 / (1 + exp(- X_augmented %*% beta_0))
p_1 = 1 / (1 + exp(- X_augmented %*% beta_1))

y_0_complete = rbinom(n = n, size = 1, prob = p_0)
y_1_complete = rbinom(n = n, size = 1, prob = p_1)
y = rep(0, n)
y[A == 0] = y_0_complete[A == 0]
y[A == 1] = y_1_complete[A == 1]

source("bayesian_logistic_regression_for_causal_inference.R")
set.seed(2022)
blr = bayesian_causal_inference_logistic_regression(X, y, A, 
                                                    a_0 = 1, b_0 = 1, a_1 = 1, b_1 = 1, 
                                                    total_steps = 1000)

mean(y_1_complete - y_0_complete)
hist(blr$tau_p_samples)
mean(blr$tau_p_samples)

beta_0
colMeans(blr$beta_0_samples)
p_0_est = 1 / (1 + exp(- X_augmented %*% colMeans(blr$beta_0_samples)))

beta_1
colMeans(blr$beta_1_samples)
p_1_est = 1 / (1 + exp(- X_augmented %*% colMeans(blr$beta_1_samples)))

plot(p_0, p_0_est)
plot(p_1, p_1_est)

