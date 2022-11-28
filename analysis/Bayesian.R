################################################################################
# a subset of covariates
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
# all covariates
games = read.csv("games_200_22_complete_clean.csv")
games = games[, -1]
games$medium = factor(games$medium)
games$time_control = relevel(factor(games$time_control), ref = "standard")
games$home_country_diff = relevel(factor(games$home_country_diff), ref = "0")
games$players_sex = relevel(factor(games$players_sex), ref = "MM")
games$first_game_diff = relevel(factor(games$first_game_diff), ref = "0")
cleaned_games = model.matrix(~ . - 1, data = games)
cleaned_games = cleaned_games[, -1]
cleaned_games = as.data.frame(cleaned_games)

A = cleaned_games$mediumonline
X = as.matrix(cleaned_games[, 2:19])
X[, c("age_diff", "rating_diff", "n_tour_diff", "n_game_diff", "score_diff", "n_games_higher", "n_upset_higher", "n_peo", "n_upset_peo")] = 
  scale(X[, c("age_diff", "rating_diff", "n_tour_diff", "n_game_diff", "score_diff", "n_games_higher", "n_upset_higher", "n_peo", "n_upset_peo")])
y = cleaned_games$upset

source("bayesian_logistic_regression_for_causal_inference.R")
set.seed(2022)
blr = bayesian_causal_inference_logistic_regression(X, y, A, 
                                                    a_0 = 1, b_0 = 1, a_1 = 1, b_1 = 1, 
                                                    total_steps = 10000)
hist(blr$tau_p_samples[5001:1000])
quantile(blr$tau_p_samples[5001:10000], c(0.025, 0.975))
save(blr, file = "blr.rdata")

beta_0 = colMeans(blr$beta_0_samples)
beta_1 = colMeans(blr$beta_1_samples)
hist(1 / (1 + exp(-cbind(1, X) %*% beta_0)))
hist(1 / (1 + exp(-cbind(1, X) %*% beta_1)))