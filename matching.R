library(tidyverse)
library(MatchIt)
library(fastDummies)
library(broom)
library(marginaleffects)

data <- read.csv('~/Desktop/biostat681_final_project/games_200_22_complete_clean.csv', stringsAsFactors = TRUE)
data <- data[,2:ncol(data)]
data$medium <- ifelse(data$medium == "online", 1, 0) ### 1 = online; 0 = offline

## Check initial imbalance
m.out0 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = NULL, distance = "mahalanobis")

summary(m.out0)

#### ESTIMATE ACET
m.out1 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "nearest", distance = "glm", replace = FALSE,
                  estimand = "ATT")
summary(m.out1, un = FALSE)
plot(summary(m.out1))

m.out2 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "nearest", distance = "glm", replace = TRUE,
                  estimand = "ATT")
summary(m.out2, un = FALSE)
plot(summary(m.out2))

m.out3 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "nearest", distance = "mahalanobis", replace = FALSE,
                  estimand = "ATT")
summary(m.out3, un = FALSE)
plot(summary(m.out3))

m.out4 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "nearest", distance = "mahalanobis", replace = TRUE,
                  estimand = "ATT")
summary(m.out4, un = FALSE)
plot(summary(m.out4))

m.out5 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "full", distance = "glm",
                  estimand = "ATT")
summary(m.out5, un = FALSE)
plot(summary(m.out5))

m.out6 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "full", distance = "mahalanobis",
                  estimand = "ATT")
summary(m.out6, un = FALSE)
plot(summary(m.out6))

smd0 <- summary(m.out0)$sum.all[,"Std. Mean Diff."]
smd0 <- unname(smd0)
smd1 <- summary(m.out1)$sum.matched[,"Std. Mean Diff."]
smd1 <- unname(smd1[2:length(smd1)])
smd2 <- summary(m.out2)$sum.matched[,"Std. Mean Diff."]
smd2 <- unname(smd2[2:length(smd2)])
smd3 <- summary(m.out3)$sum.matched[,"Std. Mean Diff."]
smd3 <- unname(smd3)
smd4 <- summary(m.out4)$sum.matched[,"Std. Mean Diff."]
smd4 <- unname(smd4)
smd5 <- summary(m.out5)$sum.matched[,"Std. Mean Diff."]
smd5 <- unname(smd5[2:length(smd5)])
smd6 <- summary(m.out6)$sum.matched[,"Std. Mean Diff."]
smd6 <- unname(smd6)

smd_summary <- data.frame(variables = names(summary(m.out0)$sum.all[,"Std. Mean Diff."]),
                          '0' = round(smd0,2),
                          '1' = round(smd1,2),
                          '2' = round(smd2,2),
                          '3' = round(smd3,2),
                          '4' = round(smd4,2),
                          '5' = round(smd5,2),
                          '6' = round(smd6,2))

### METHOD 2
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #Do 1:1 PS matching with replacement
  m <- matchit(medium ~ time_control + age_diff + home_country_diff +
                 rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                 score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
               data = boot_data, method = "nearest", distance = "glm", replace = TRUE,
               estimand = "ATT")
  
  #Extract matched dataset
  md <- match.data(m, data = boot_data)
  
  #Fit outcome model
  fit <- lm(upset ~ medium * (time_control + age_diff + home_country_diff +
                                 rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                                 score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo),
             data = md, weights = weights)
  
  ## G-computation ##
  #Subset to treated units for ATT; skip for ATE
  md1 <- subset(md, medium == 1)
  
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, medium = 1))
  Ep1 <- weighted.mean(p1, md1$weights)
  
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, medium = 0))
  Ep0 <- weighted.mean(p0, md1$weights)
  
  #Risk ratio
  return(Ep1 - Ep0)
}

library("boot")
set.seed(2022)
boot_out <- boot(data, boot_fun, R = 1000)
boot_out
boot.ci(boot_out, type = "perc")


### METHOD 4
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #Do 1:1 PS matching with replacement
  m <- matchit(medium ~ time_control + age_diff + home_country_diff +
                 rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                 score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
               data = boot_data, method = "nearest", distance = "mahalanobis", replace = TRUE,
               estimand = "ATT")
  
  #Extract matched dataset
  md <- match.data(m, data = boot_data)
  
  #Fit outcome model
  fit <- lm(upset ~ medium * (time_control + age_diff + home_country_diff +
                                rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                                score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo),
            data = md, weights = weights)
  
  ## G-computation ##
  #Subset to treated units for ATT; skip for ATE
  md1 <- subset(md, medium == 1)
  
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, medium = 1))
  Ep1 <- weighted.mean(p1, md1$weights)
  
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, medium = 0))
  Ep0 <- weighted.mean(p0, md1$weights)
  
  #Risk ratio
  return(Ep1 - Ep0)
}

library("boot")
set.seed(2022)
boot_out <- boot(data, boot_fun, R = 1000)
boot_out
boot.ci(boot_out, type = "perc")


### METHOD 5
m.out5 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "full", distance = "glm",
                  estimand = "ATT")
m.data <- match.data(m.out5)
fit <- lm(upset ~ medium * (time_control + age_diff + home_country_diff +
                              rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                              score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo), 
          data = m.data, weights = weights)

comp <- comparisons(fit,
                    variables = "medium",
                    vcov = ~subclass,
                    newdata = subset(m.data, medium == 1),
                    wts = "weights")
summary(comp)


### METHOD 5 for ACE
m.out5 <- matchit(medium ~ time_control + age_diff + home_country_diff +
                    rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                    score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo,
                  data = data, method = "full", distance = "glm",
                  estimand = "ATE")
summary(m.out5)
m.data <- match.data(m.out5)
fit <- lm(upset ~ medium * (time_control + age_diff + home_country_diff +
                              rating_diff + players_sex + first_game_diff + n_tour_diff + n_game_diff + 
                              score_diff + n_games_higher + n_upset_higher + n_peo + n_upset_peo), 
          data = m.data, weights = weights)

comp <- comparisons(fit,
                    variables = "medium",
                    vcov = ~subclass,
                    wts = "weights")
summary(comp)