library(testthat)
library(dplyr)

# Test FitMarkovSwitchingGAMLSS function
test_that("FitMarkovSwitchingGAMLSS works correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  y <- rnorm(50)
  formula <- y ~ x1 + x2
  result <- FitMarkovSwitchingGAMLSS(x, y, formula = formula)
  
  expect_type(result, "list")
  expect_true("delta" %in% names(result))
  expect_true("gamma" %in% names(result))
  expect_true("mod" %in% names(result))
  expect_true("llh" %in% names(result))
})

# Test PredictStateSequence function
test_that("PredictStateSequence works correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  y <- rnorm(50)
  formula <- y ~ x1 + x2
  model <- FitMarkovSwitchingGAMLSS(x, y, formula = formula)
  predicted_states <- PredictStateSequence(model, y)
  
  expect_type(predicted_states, "integer")
  expect_length(predicted_states, length(y))
})

# Test compute_likelihood function
test_that("compute_likelihood works correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  y <- rnorm(50)
  formula <- y ~ x1 + x2
  model <- FitMarkovSwitchingGAMLSS(x, y, formula = formula)
  test_data <- data.frame(y = y, x = x)
  likelihood <- compute_likelihood(test_data, model$delta, model$gamma, model)
  
  expect_type(likelihood, "double")
  expect_true(likelihood > 0)
})

# Test cv_msgamlss function
test_that("cv_msgamlss works correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  y <- rnorm(50)
  states <- sample(1:2, 50, replace = TRUE)
  init_probs_list <- list(matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2))
  formula <- y ~ x1 + x2
  result <- cv_msgamlss(x, y, states, init_probs_list, formula = formula)
  
  expect_type(result, "list")
  expect_true("model_results" %in% names(result))
  expect_true("final_probs" %in% names(result))
  expect_true("accuracy" %in% names(result))
  expect_true("likelihood" %in% names(result))
})

# Test cv_results function
test_that("cv_results works correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), ncol = 2)
  y <- rnorm(50)
  states <- sample(1:2, 50, replace = TRUE)
  init_probs_list <- list(matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2))
  formula <- y ~ x1 + x2
  cv_result <- cv_msgamlss(x, y, states, init_probs_list, formula = formula)
  true_gamma <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2)
  result <- cv_results(cv_result, init_probs_list, true_gamma, 2)
  
  expect_type(result, "list")
  expect_true("gamma_df" %in% names(result))
  expect_true("results_df" %in% names(result))
  expect_true("agg_results_df" %in% names(result))
})
