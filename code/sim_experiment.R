source("code/functions.R")

init_state_probs_list <- list(
  matrix(c(0.99, 0.01, 0.01, 0.99), nrow = 2),
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
  matrix(c(0.05, 0.95, 0.95, 0.05), nrow = 2),
  matrix(c(0.01, 0.99, 0.99, 0.01), nrow = 2))

true_gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2)
k <- 2
x <- linear_data[,2:101]
y <- linear_data$y
states <- nonlinear_data$states
data <- data.frame(y = y, x = x, states = states)
formula <-as.formula(paste("y ~", paste("bols(", names(x), ")", collapse = " + ")))
formula

results <- cross_validate_msgamlss(x = x, y=y, init_probs_list = init_state_probs_list, formula = formula)

k <- 20
x <- nonlinear_data[,2]
y <- nonlinear_data$y
states <- nonlinear_data$states
data <- data.frame(y = y, x = x, states = states)

results_nonlin <- CrossValidateMarkovSwitchingGAMLSS(data, k = k, init_state_probs_list = init_state_probs_list, type = "MSGAMLSS")

save(results_nonlin, file = "results_nonlin.RData")
print(results)

#Extracting result data
gamma_list <- matrix(NA, length(init_state_probs_list)*k, 5)
i <- 1
while(i <= length(init_state_probs_list)*k) {
  for(model in 1:length(init_state_probs_list)){
    for(j in 1:k){
      gamma_list[i,1] <- results$model_results[[model]][[j]]$gamma[1,1]
      gamma_list[i,2] <- results$model_results[[model]][[j]]$gamma[1,2]
      gamma_list[i,3] <- results$model_results[[model]][[j]]$gamma[2,1]
      gamma_list[i,4] <- results$model_results[[model]][[j]]$gamma[2,2]
      gamma_list[i,5] <- (abs(gamma_list[i,1] - true_gamma[1,1]) + abs(gamma_list[i,2] - true_gamma[1,2]) + abs(gamma_list[i,3] - true_gamma[2,1]) + abs(gamma_list[i,4] - true_gamma[2,2]))/4
      i <- i + 1
    }
  }
}
colnames(gamma_list) <- c("gamma_11", "gamma_12", "gamma_21", "gamma_22", "abs_diff")
gamma_df <- as.data.frame(gamma_list)
gamma_df <- reshape2::melt(gamma_df)
print(gamma_df)

gammanonlin_list <- matrix(NA, length(init_state_probs_list)*k, 5)
i <- 1
while(i <= length(init_state_probs_list)*k) {
  for(model in 1:length(init_state_probs_list)){
    for(j in 1:k){
      gammanonlin_list[i,1] <- results_nonlin$model_results[[model]][[j]]$gamma[1,1]
      gammanonlin_list[i,2] <- results_nonlin$model_results[[model]][[j]]$gamma[1,2]
      gammanonlin_list[i,3] <- results_nonlin$model_results[[model]][[j]]$gamma[2,1]
      gammanonlin_list[i,4] <- results_nonlin$model_results[[model]][[j]]$gamma[2,2]
      gammanonlin_list[i,5] <- (abs(gammanonlin_list[i,1] - true_gamma[1,1]) + abs(gammanonlin_list[i,2] - true_gamma[1,2]) + abs(gammanonlin_list[i,3] - true_gamma[2,1]) + abs(gammanonlin_list[i,4] - true_gamma[2,2]))/4
      i <- i + 1
    }
  }
}
colnames(gammanonlin_list) <- c("gamma_11", "gamma_12", "gamma_21", "gamma_22", "abs_dif")
gammanonlin_df <- as.data.frame(gammanonlin_list)
gammanonlin_df <- reshape2::melt(gammanonlin_df)
print(gammanonlin_df)

# Extracting influence estimates for x1


results_df <- data.frame(
  init_state_probs = rep(sapply(init_state_probs_list, function(x) paste(x, collapse = ",")), each = 10),
  gamma_11 = gamma_list[, 1],
  gamma_12 = gamma_list[, 2],
  gamma_21 = gamma_list[, 3],
  gamma_22 = gamma_list[, 4],
  abs_dif = gamma_list[, 5],
  mse = unlist(results$mse),
  accuracy = unlist(results$accuracy),
  final_i = unlist(results$final_iterations)
)

print(results_df)

results_nonlin_df <- data.frame(
  init_state_probs = rep(sapply(init_state_probs_list, function(x) paste(x, collapse = ",")), each = 10),
  gamma_11 = gammanonlin_list[, 1],
  gamma_12 = gammanonlin_list[, 2],
  gamma_21 = gammanonlin_list[, 3],
  gamma_22 = gammanonlin_list[, 4],
  abs_diff = gammanonlin_list[, 5],
  mse = unlist(results_nonlin$mse),
  accuracy = unlist(results_nonlin$accuracy),
  final_i = unlist(results_nonlin$final_iterations)
)

print(results_nonlin_df)


# Aggregating MSE, Accuracy, and final_i for each init_state_prob
agg_results_df <- results_df %>%
  group_by(init_state_probs) %>%
  summarise(
    mean_mse = mean(mse),
    mean_accuracy = mean(accuracy),
    mean_final_i = round(mean(final_i)),
    mean_abs_diff = mean(abs_dif)
  )
print(agg_results_df)

agg_results_nonlin_df <- results_nonlin_df %>%
  group_by(init_state_probs) %>%
  summarise(
    mean_mse = mean(mse),
    mean_accuracy = mean(accuracy),
    mean_final_i = round(mean(final_i)),
    mean_abs_diff = mean(abs_diff)
  )
print(agg_results_nonlin_df)
#Plots 
# Plotting linear data
ggplot(linear_data, aes(x = X1, y = y, color = as.factor(states))) +
  geom_point() +
  labs(title = "Linear Data: Y vs. X", x = "X", y = "Y", color = "States") +
  theme_minimal()

# Plotting nonlinear data
data_nonlin <- data.frame(y = nonlinear_data$y, x = nonlinear_data[,2], states = nonlinear_data$states)
ggplot(data_nonlin, aes(x = x, y = y, color = as.factor(states))) +
  geom_point() +
  labs(title = "Nonlinear Data: Y vs. X", x = "X", y = "Y", color = "States") +
  theme_minimal()

# Plotting y against time for linear data
linear_data$time <- 1:nrow(data)
ggplot(linear_data, aes(x = time, y = y)) +
  geom_line(aes(color = as.factor(states), group = 1)) +
  geom_point(aes(color = as.factor(states))) +
  labs(title = "Linear Data", x = "Time", y = "Y", color = "States") +
  theme_minimal()

# Plotting y against time for nonlinear data
data_nonlin$time <- 1:nrow(data_nonlin)
ggplot(data_nonlin, aes(x = time, y = y)) +
  geom_line(aes(color = as.factor(states), group = 1)) +
  geom_point(aes(color = as.factor(states))) +
  labs(title = "Nonlinear Data", x = "Time", y = "Y", color = "States") +
  theme_minimal()


gamma_df$init_state_probs <- rep(sapply(init_state_probs_list, function(x) paste(x, collapse = ",")), each = 10)

ggplot(gamma_df, aes(x = variable, y = value, color = init_state_probs)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs( y = "State Transition probabilities", x ="", color = "Initial state transition probs.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, color = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")),
        legend.title = element_text(size = 14)) +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22])))

gammanonlin_df$init_state_probs <- rep(sapply(init_state_probs_list, function(x) paste(x, collapse = ",")), each = 10)

ggplot(gammanonlin_df, aes(x = variable, y = value, color = init_state_probs)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(x = "", y = "State transition probability", color = "Initial state transition probs.") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, color = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")),
        legend.title = element_text(size = 14)) +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22])))

# Plot 1: Boxplots for each column in gamma_df
ggplot(gamma_df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(x = "", y = "State transition probability") +
  theme_minimal() +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22]))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))

ggplot(gammanonlin_df, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable)) +
  scale_fill_manual(values = c("#FF6961", "#6A8DA5", "#6A8DA5", "#FF6961")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#FF6961") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#6A8DA5") +
  labs(x = "", y = "State transition probability") +
  theme_minimal() +
  scale_x_discrete(labels = c(expression(gamma[11]), expression(gamma[12]), expression(gamma[21]), expression(gamma[22]))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))

# Plot 2: MSE for each init_state_prob
ggplot(agg_results_df, aes(x = init_state_probs, y = mean_mse)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean MSE for Each Initial State Probability", x = "Initial State Probability", y = "Mean MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_mse)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean MSE for Each Initial State Probability", x = "Initial State Probability", y = "Mean MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot 3: Accuracy for each init_state_prob
ggplot(agg_results_df, aes(x = init_state_probs, y = mean_accuracy)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Accuracy for Each Initial State Probability", x = "Initial State Probability", y = "Mean Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_accuracy)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Accuracy for Each Initial State Probability", x = "Initial State Probability", y = "Mean Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot 4: Final iterations for each init_state_prob
ggplot(agg_results_df, aes(x = init_state_probs, y = mean_final_i)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Final Iterations for Each Initial State Probability", x = "Initial State Probability", y = "Mean Final Iterations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_final_i)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Final Iterations for Each Initial State Probability", x = "Initial State Probability", y = "Mean Final Iterations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plot5
ggplot(agg_results_df, aes(x = init_state_probs, y = mean_abs_diff)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean MSE for Each Initial State Probability", x = "Initial State Probability", y = "Mean MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(agg_results_nonlin_df, aes(x = init_state_probs, y = mean_abs_diff)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean MSE for Each Initial State Probability", x = "Initial State Probability", y = "Mean MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

x <- linear_data[,2:101]
y <- linear_data$y
states <- linear_data$states
data <- data.frame(y = y, x = x, states = states)

test <- FitMarkovSwitchingGAMLSS(
  x = as.matrix(x),
  y = linear_data$y,
  N = 2,
  type = "MSGAMLSS",
  stat = FALSE,
  m.stop = rep(100, 2),
  max.iter = 10,
  conv.tol = 1e-02,
  conv.print = TRUE,
  formula = formula)
formula <-as.formula(paste("y ~", paste("bols(", names(x), ")", collapse = " + ")))
formula
