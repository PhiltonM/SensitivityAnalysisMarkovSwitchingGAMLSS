# Set the working directory to the project root or appropriate directory
setwd("c:/Users/philt/OneDrive - stud.uni-goettingen.de/Studium/Master/WiSe 2425/CTA/SensitivityAnalysisMarkovSwitchingGAMLSS")

# Ensure the test directory path is correct
testthat::test_dir("code")

if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
renv::activate()
required_packages <- c("fda", "boot", "numDeriv", "gamlss", "gamboostLSS", "ggplot2", "dplyr", "reshape2")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}
library(fda)
library(boot)
library(numDeriv)
library(gamlss)
library(gamboostLSS)
library(ggplot2)
library(dplyr)
library(reshape2)

source("code/functions.R")
source("code/sim_data.R")

init_state_probs_list <- list(
  matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2),
  matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2))

true_gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2)
k <- 2
x_lin <- linear_data[,2:101]
y_lin <- linear_data$y
states_lin <- linear_data$states
data <- data.frame(y = y_lin, x = x_lin, states = states_lin)
formula <-as.formula(paste("y ~", paste( names(x_lin), collapse = " + ")))
formula

cv_linear <- cv_msgamlss(x = x_lin, y=y_lin, init_probs_list = init_state_probs_list, formula = formula, states = states_lin, type ="MSGLMLSS")

x_nonlin <- nonlinear_data[,2:101]
y_nonlin <- nonlinear_data$y
states_nonlin <- nonlinear_data$states
data_nonlin <- data.frame(y = y_nonlin, x = x_nonlin, states = states_nonlin)
formula <-as.formula(paste("y ~", paste("bbs(", names(x_nonlin), ")", collapse = " + ")))
formula

cv_nonlinear <- cv_msgamlss(x = x_nonlin, y=y_nonlin, init_probs_list = init_state_probs_list, formula = formula, states = states_nonlin, type ="MSGAMLSS")


cv_res_linear <- cv_results(cv_linear, init_state_probs_list, true_gamma, k)
gamma_lin_df <- cv_res_linear$gamma_df
results_lin_df <- cv_res_linear$results_df
agg_results_lin_df <- cv_res_linear$agg_results_df
print(agg_results_lin_df)
print(gamma_lin_df)
print(results_lin_df)

cv_res_nonlinear <- cv_results(cv_nonlinear, init_state_probs_list, true_gamma, k)
gamma_nonlin_df <- cv_res_nonlinear$gamma_df
results_nonlin_df <- cv_res_nonlinear$results_df
agg_results_nonlin_df <- cv_res_nonlinear$agg_results_df
print(agg_results_nonlin_df)
print(gamma_nonlin_df)
print(results_nonlin_df)

ggplot(linear_data, aes(x = X1, y = y, color = as.factor(states))) +
  geom_point() +
  labs(title = "Linear Data: Y vs. X", x = "X", y = "Y", color = "States") +
  theme_minimal()

linear_data$time <- 1:nrow(data)
ggplot(linear_data, aes(x = time, y = y)) +
  geom_line(aes(color = as.factor(states), group = 1)) +
  geom_point(aes(color = as.factor(states))) +
  labs(title = "Linear Data", x = "Time", y = "Y", color = "States") +
  theme_minimal()

gamma_lin_df$init_state_probs <- rep(sapply(init_state_probs_list, function(x_lin) paste(x_lin, collapse = ",")), each = 10)

ggplot(gamma_lin_df, aes(x = variable, y = value, color = init_state_probs)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs( y = "State Transition probabilities", x ="", color = "Initial state transition probs.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, color = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")),
        legend.title = element_text(size = 14)) +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22])))

ggplot(gamma_lin_df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(x = "", y = "State transition probability") +
  theme_minimal() +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22]))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))

ggplot(agg_results_lin_df, aes(x = init_state_probs, y = mean_likelihood)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Likelihood for Each Initial State Probability", x = "Initial State Probability", y = "Mean Likelihood") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


coefficients_list <- lapply(cv_linear$model_results, function(res) {
  lapply(res, function(mod) {
    coef(mod$mod[[1]]$mu)
  })
})

coefficients_df <- do.call(rbind, lapply(seq_along(coefficients_list), function(i) {
  do.call(rbind, lapply(seq_along(coefficients_list[[i]]), function(j) {
    data.frame(variable = names(coefficients_list[[i]][[j]]),
                value = coefficients_list[[i]][[j]],
                init_state_prob = paste(init_state_probs_list[[i]], collapse = ","))
  }))
}))

coefficients_df <- coefficients_df %>%
  group_by(init_state_prob, variable) %>%
  summarise(mean_value = mean(value), .groups = 'drop')

coefficients_df <- coefficients_df %>%   
  mutate(variable = ifelse(variable == "X1", "X1", "X2,...,X100")) %>%
  group_by(init_state_prob, variable) %>%
  summarise(mean_value = sum(mean_value), .groups = 'drop')

ggplot(coefficients_df, aes(x = variable, y = mean_value, fill = init_state_prob)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Coefficients from Cross-Validation", x = "Variable", y = "Mean Coefficient Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
coefficients_df


nonlin_coefficients_list <- lapply(cv_nonlinear$model_results, function(res) {
  lapply(res, function(mod) {
    coef(mod$mod[[1]]$mu)
  })
})
nonlin_coefficients_list

nonlin_coefficients_df <- do.call(rbind, lapply(seq_along(nonlin_coefficients_list), function(i) {
  do.call(rbind, lapply(seq_along(nonlin_coefficients_list[[i]]), function(j) {
    data.frame(variable = names(nonlin_coefficients_list[[i]][[j]]),
                value = nonlin_coefficients_list[[i]][[j]],
                init_state_prob = paste(init_state_probs_list[[i]], collapse = ","))
  }))
}))

nonlin_coefficients_df <- nonlin_coefficients_df %>%
  group_by(init_state_prob, variable) %>%
  summarise(mean_value = mean(value), .groups = 'drop')

nonlin_coefficients_df <- nonlin_coefficients_df %>%   
  mutate(variable = ifelse(variable == "X1", "X1", "X2,...,X100")) %>%
  group_by(init_state_prob, variable) %>%
  summarise(mean_value = sum(mean_value), .groups = 'drop')

ggplot(coefficients_df, aes(x = variable, y = mean_value, fill = init_state_prob)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Coefficients from Cross-Validation", x = "Variable", y = "Mean Coefficient Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
nonlin_coefficients_df
