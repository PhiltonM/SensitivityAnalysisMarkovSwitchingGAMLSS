# Load necessary libraries
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

# Source additional R scripts
source("code/functions.R")
source("code/sim_data.R")

# Define initial state probabilities and true gamma
init_probs_list <- list(
  matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2),
  matrix(c(0.85, 0.15, 0.15, 0.85), nrow = 2),
  matrix(c(0.75, 0.25, 0.25, 0.75), nrow = 2),
  matrix(c(0.65, 0.35, 0.35, 0.65), nrow = 2),
  matrix(c(0.55, 0.45, 0.45, 0.55), nrow = 2),
  matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2),
  matrix(c(0.45, 0.55, 0.55, 0.45), nrow = 2),
  matrix(c(0.35, 0.65, 0.65, 0.35), nrow = 2),
  matrix(c(0.25, 0.75, 0.75, 0.25), nrow = 2),
  matrix(c(0.15, 0.85, 0.85, 0.15), nrow = 2),
  matrix(c(0.05, 0.95, 0.95, 0.05), nrow = 2)
)

true_gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2)
k <- 2

# Prepare linear data
x_lin <- linear_data[, 2:101]
y_lin <- linear_data$y
states_lin <- linear_data$states
data <- data.frame(y = y_lin, x = x_lin, states = states_lin)
formula <- as.formula(paste("y ~", paste(names(x_lin), collapse = " + ")))

# Perform cross-validation for linear data
cv_linear <- cv_msgamlss(x = x_lin, y = y_lin, init_probs_list = init_probs_list, formula = formula, states = states_lin, type = "MSGLMLSS", max_iter = 100, mstop = c(500, 500))

# Prepare nonlinear data
x_nonlin <- nonlinear_data[, 2:101]
y_nonlin <- nonlinear_data$y
states_nonlin <- nonlinear_data$states
data_nonlin <- data.frame(y = y_nonlin, x = x_nonlin, states = states_nonlin)
formula <- as.formula(paste("y ~", paste("bbs(", names(x_nonlin), ")", collapse = " + ")))

# Perform cross-validation for nonlinear data
cv_nonlinear <- cv_msgamlss(x = x_nonlin, y = y_nonlin, init_probs_list = init_probs_list, formula = formula, states = states_nonlin, type = "MSGAMLSS", max_iter = 400, mstop = c(160, 160))

# Process results for linear data
cv_res_linear <- cv_results(cv_linear, init_probs_list, true_gamma, k)
gamma_lin_df <- cv_res_linear$gamma_df
MAE_lin <- gamma_lin_df[105:130, ]
gamma_lin_df <- gamma_lin_df[1:104, ]
results_lin_df <- cv_res_linear$results_df
agg_results_lin_df <- cv_res_linear$agg_results_df

# Calculate absolute differences
abs_diff <- numeric(nrow(agg_results_lin_df))
for (i in 1:nrow(agg_results_lin_df)) {
  abs_diff[i] <- sum(abs(init_probs_list[[i]][1, 1] - true_gamma[1, 1]), abs(init_probs_list[[i]][1, 2] - true_gamma[1, 2]), abs(init_probs_list[[i]][2, 1] - true_gamma[2, 1]), abs(init_probs_list[[i]][2, 2] - true_gamma[2, 2])) / 4
}
abs_diff <- rev(abs_diff)
agg_results_lin_df$abs_diff <- abs_diff

# Print results for linear data
print(agg_results_lin_df)
print(gamma_lin_df)
print(results_lin_df)

# Process results for nonlinear data
cv_res_nonlinear <- cv_results(cv_nonlinear, init_probs_list, true_gamma, k)
gamma_nonlin_df <- cv_res_nonlinear$gamma_df
MAE_nonlin <- gamma_nonlin_df[105:130, ]
gamma_nonlin_df <- gamma_nonlin_df[1:104, ]
results_nonlin_df <- cv_res_nonlinear$results_df
agg_results_nonlin_df <- cv_res_nonlinear$agg_results_df
agg_results_nonlin_df$abs_diff <- abs_diff

# Print results for nonlinear data
print(agg_results_nonlin_df)
print(gamma_nonlin_df)
print(results_nonlin_df)

# Plot results for linear data
ggplot(linear_data, aes(x = X1, y = y, color = as.factor(states))) +
  geom_point() +
  labs(title = "", x = "X", y = "Y", color = "States") +
  theme_minimal()

linear_data$time <- 1:nrow(data)
ggplot(linear_data, aes(x = time, y = y)) +
  geom_line(aes(color = as.factor(states), group = 1)) +
  geom_point(aes(color = as.factor(states))) +
  labs(title = "", x = "Time", y = "Y", color = "States") +
  theme_minimal()

gamma_lin_df$init_state_probs <- rep(sapply(init_probs_list, function(x_lin) paste(x_lin, collapse = ",")), each = k)

ggplot(gamma_lin_df, aes(x = variable, y = value, color = init_state_probs)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(y = "State Transition probabilities", x = "", color = "Initial state transition probs.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, color = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")),
        legend.title = element_text(size = 12)) +
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

ggplot(agg_results_lin_df, aes(x = init_state_probs, y = mean_accuracy)) +
  geom_bar(stat = "identity", aes(fill = "Accuracy")) +
  geom_point(aes(y = abs_diff, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 3) +
  geom_line(aes(y = abs_diff, group = 1, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 1) +
  labs(title = "", x = "Initial State Probability", y = "Accuracy / Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(name = "", values = c("Accuracy" = "#6A8DA5")) +
  scale_color_manual(name = "", values = c("Absolute Difference \nbetween True Gamma \nand Initial guess" = "#FF6961"))

ggplot(agg_results_lin_df, aes(x = init_state_probs, y = mean_abs_diff)) +
  geom_bar(stat = "identity", aes(fill = "MAE per Parameter")) +
  geom_point(aes(y = abs_diff, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 3) +
  geom_line(aes(y = abs_diff, group = 1, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 1) +
  labs(title = "", x = "Initial State Probability", y = "MAE / Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(name = "", values = c("MAE per Parameter" = "#6A8DA5")) +
  scale_color_manual(name = "", values = c("Absolute Difference \nbetween True Gamma \nand Initial guess" = "#FF6961"))

ggplot(agg_results_lin_df, aes(x = init_state_probs, y = mean_final_i)) +
  geom_bar(stat = "identity", fill = "#6A8DA5") +
  labs(title = "", x = "Initial State Probability", y = "Mean Convergence iteration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot results for nonlinear data
ggplot(nonlinear_data, aes(x = X1, y = y, color = as.factor(states))) +
  geom_point() +
  labs(title = "", x = "X", y = "Y", color = "States") +
  theme_minimal()

nonlinear_data$time <- 1:nrow(data)
ggplot(nonlinear_data, aes(x = time, y = y)) +
  geom_line(aes(color = as.factor(states), group = 1)) +
  geom_point(aes(color = as.factor(states))) +
  labs(title = "", x = "Time", y = "Y", color = "States") +
  theme_minimal()

gamma_nonlin_df$init_state_probs <- rep(sapply(init_probs_list, function(x_lin) paste(x_lin, collapse = ",")), each = k)

ggplot(gamma_nonlin_df, aes(x = variable, y = value, color = init_state_probs)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(y = "State Transition probabilities", x = "", color = "Initial state transition probs.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, color = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")),
        legend.title = element_text(size = 12)) +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22])))

ggplot(gamma_nonlin_df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(x = "", y = "State transition probability") +
  theme_minimal() +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22]))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_accuracy)) +
  geom_bar(stat = "identity", aes(fill = "Accuracy")) +
  geom_point(aes(y = abs_diff, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 3) +
  geom_line(aes(y = abs_diff, group = 1, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 1) +
  labs(title = "", x = "Initial State Probability", y = "Accuracy / Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(name = "", values = c("Accuracy" = "#6A8DA5")) +
  scale_color_manual(name = "", values = c("Absolute Difference \nbetween True Gamma \nand Initial guess" = "#FF6961"))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_abs_diff)) +
  geom_bar(stat = "identity", aes(fill = "MAE per Parameter")) +
  geom_point(aes(y = abs_diff, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 3) +
  geom_line(aes(y = abs_diff, group = 1, color = "Absolute Difference \nbetween True Gamma \nand Initial guess"), size = 1) +
  labs(title = "", x = "Initial State Probability", y = "MAE / Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(name = "", values = c("MAE per Parameter" = "#6A8DA5")) +
  scale_color_manual(name = "", values = c("Absolute Difference \nbetween True Gamma \nand Initial guess" = "#FF6961"))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_final_i)) +
  geom_bar(stat = "identity", fill = "#6A8DA5") +
  labs(title = "", x = "Initial State Probability", y = "Mean Convergence iteration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


