## Function that fits Markov-switching GAMLSS using the msgamboostLSS algorithm
FitMarkovSwitchingGAMLSS <- function(x,
                                     y,
                                     N = 2,
                                     type = "MSGAMLSS",
                                     stat = FALSE,
                                     m.stop = rep(480, 2),
                                     max.iter = 50,
                                     conv.tol = 1e-02,
                                     conv.print = TRUE,
                                     init_state_trans_prob = NULL,
                                     formula = NULL) { 

  delta = NULL
  if(is.null(init_state_trans_prob)){
    gamma = matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2)
  }
  else{
    gamma = as.matrix(init_state_trans_prob)
  }
  mod = fv.mu = fv.sigma = list()
  term = FALSE
  old = 0
  circ.crit = 0
  for(j in 1:N) {
    fv.mu[[j]] = rep(c(3, 7)[j], length(y))
    fv.sigma[[j]] = rep(3, length(y))
  }
  while(term == FALSE) {
    for(i in 1:max.iter) {
      delta.next = delta
      gamma.next = gamma
      fv.mu.next = fv.mu
      fv.sigma.next = fv.sigma
      if(is.null(delta)) {
        delta = solve(t(diag(N) - gamma + 1), rep(1, N))
      }
      allprobs = matrix(NA, length(y), N)
      lalpha = lbeta = matrix(NA, N, length(y))
      for(j in 1:N) {
        allprobs[,j] = apply(as.matrix(y), 2, dNO, mu = fv.mu.next[[j]], sigma=exp(fv.sigma.next[[j]]))
      }
      allprobs = ifelse(!is.na(allprobs), allprobs, 1)
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      lscale = log(sumfoo)
      foo = foo / sumfoo
      lalpha[, 1] = log(foo) + lscale
      for(j in 2:length(y)) {
        foo = foo %*% gamma * allprobs[j,]
        sumfoo = sum(foo)
        lscale = lscale + log(sumfoo)
        foo = foo / sumfoo
        lalpha[, j] = log(foo) + lscale
      }
      foo = rep(1 / N, N)
      lbeta[, length(y)] = rep(0, N)
      lscale = log(N)
      foo = foo / sum(foo)
      for(j in (length(y) - 1):1) {
        foo = gamma %*% (allprobs[j + 1,] * foo)
        lbeta[, j] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)
      }                   
      lallprobs = log(allprobs)
      llh = max(lalpha[, length(y)]) + log(sum(exp(lalpha[, length(y)] - max(lalpha[, length(y)]))))
      weights = matrix(NA, N, length(y))
      for(j in 1:length(y)) {
        weights[,j] = exp(lalpha[,j] + lbeta[, j] - llh)
      }
      for(j in 1:N) {
        for(k in 1:N) {
          gamma.next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(length(y) - 1)] + lallprobs[2:length(y), k] + lbeta[k, 2:length(y)] - llh))
        }
      }
      gamma.next = gamma.next / apply(gamma.next, 1, sum)
      if(stat == TRUE) {
        delta.next = solve(t(diag(N) - gamma.next + 1), rep(1, N))
      }else{
        delta.next = exp(lalpha[, 1] + lbeta[, 1] - llh)
        delta.next = delta.next / sum(delta.next)
      }
      conv.crit = sum(abs(gamma - gamma.next)) + sum(abs(delta - delta.next))
      ind = weights
      for(j in 1:N){
        if(type == "MSGLMLSS") {
            mod[[j]] = glmboostLSS(formula , weights = ind[j,], data = data.frame(x, y), families = as.families("NO"), method = "noncyclic", control = boost_control(mstop = m.stop[j], nu = 0.1))
        }
        if(type == "MSGAMLSS") {
          mod[[j]] = gamboostLSS(formula, weights = ind[j,], data = data.frame(x, y), families = as.families("NO"), method = "noncyclic", control = boost_control(mstop = m.stop[j], nu = 0.1))
        }
        fv.mu.next[[j]] = as.vector(fitted(mod[[j]]$mu))
        fv.sigma.next[[j]] = as.vector(fitted(mod[[j]]$sigma))
      }
      conv.crit = abs(llh - old)
      if(conv.print == TRUE) {
        cat("Iteration = ", i, ", criterion = ", round(conv.crit, 3), "\r", sep = "")
      }
      if(conv.crit < conv.tol | (conv.crit < 1 & abs(conv.crit - circ.crit) < 1e-09) | i == max.iter) {
        if(i == max.iter) {
          print(paste("No convergence within", max.iter, "iterations"))
        }else{
          print(paste("Convergence after", i, "iterations"))
        }
        term = TRUE
        break
      }  
      delta = delta.next
      gamma = gamma.next
      fv.mu = fv.mu.next
      fv.sigma = fv.sigma.next
      old = llh
      final_i <- i
    }
  }
  return(list(delta = delta.next, gamma = gamma.next, mod = mod, m.stop = m.stop, llh = llh, state.probs = weights, final_i = final_i))
}

## Function to predict the optimal state sequence using the Viterbi algorithm
PredictStateSequence <- function(model, y) {
  N <- length(model$delta)
  T <- length(y)
  gamma <- model$gamma
  fv.mu <- lapply(model$mod, function(m) as.vector(fitted(m$mu)))
  fv.sigma <- lapply(model$mod, function(m) as.vector(fitted(m$sigma)))
  
  delta <- log(model$delta)
  allprobs <- matrix(NA, T, N)
  
  for (j in 1:N) {
    allprobs[, j] <- apply(as.matrix(y), 2, dNO, mu = fv.mu[[j]], sigma = exp(fv.sigma[[j]]))
  }
  
  allprobs <- log(ifelse(!is.na(allprobs), allprobs, 1))
  
  viterbi <- matrix(NA, N, T)
  backpointer <- matrix(NA, N, T)
  
  viterbi[, 1] <- delta + allprobs[1, ]
  
  for (t in 2:T) {
    for (j in 1:N) {
      max_prob <- -Inf
      max_state <- NA
      for (i in 1:N) {
        prob <- viterbi[i, t - 1] + log(gamma[i, j]) + allprobs[t, j]
        if (prob > max_prob) {
          max_prob <- prob
          max_state <- i
        }
      }
      viterbi[j, t] <- max_prob
      backpointer[j, t] <- max_state
    }
  }
  
  best_path <- numeric(T)
  best_path[T] <- which.max(viterbi[, T])
  
  for (t in (T - 1):1) {
    best_path[t] <- backpointer[best_path[t + 1], t + 1]
  }
  
  return(best_path)
}

# Load necessary libraries
library(dplyr)

# Function to compute the likelihood of the test data
compute_likelihood <- function(test_data, delta, gamma, model) {
  n <- nrow(test_data)
  k <- length(delta)
  mu <- list(exp(fitted(model$mod[[1]]$mu)), fitted(model$mod[[2]]$mu))
  sigma <- list(exp(fitted(model$mod[[1]]$sigma)), exp(fitted(model$mod[[2]]$sigma)))
  # Initialize the likelihood matrix
  likelihood <- matrix(0, n, k)
  
  # Transform coefficients into parameters

  
  # Compute the initial likelihood
  for (i in 1:k) {
    likelihood[1, i] <- delta[i] * dnorm(test_data$y[1], mean = mu[[i]][[1]], sd = exp(sigma[[i]][[1]]))
  }
  
  # Compute the likelihood for the rest of the data
  for (t in 2:n) {
    for (j in 1:k) {
      likelihood[t, j] <- sum(likelihood[t-1, ] * gamma[, j]) * dnorm(test_data$y[t], mean = mu[[j]][[t]], sd = sigma[[j]][[t]])
    }
  }
  
  # Compute the total likelihood
  total_likelihood <- sum(likelihood[n, ])
  
  return(total_likelihood)
}

# K-Fold Cross Validation
cv_msgamlss <- function(x, y, states, init_probs_list, max_iter = 10, mstop = c(100, 100), type = "MSGLMLSS", formula = NULL) {
  n <- length(y)
  k <- 2
  folds <- rep(1:2, length.out = n)
  folds[seq(1, n, by = 2)] <- 1
  folds[seq(2, n, by = 2)] <- 2
  mse_list <- list()
  accuracy_list <- list()
  likelihood_list <- list()  # List to store likelihoods
  model_results_list <- list()
  final_probs_list <- list()
  final_iterations_list <- list()
  
  total_iterations <- length(init_probs_list) * k
  current_iteration <- 0
  
  for (init_probs in init_probs_list) {
    mse <- numeric(k)
    accuracy <- numeric(k)
    likelihood <- numeric(k)  # Vector to store likelihoods for each fold
    model_results <- list()
    final_probs <- list()
    final_iterations <- numeric(k)
    
    for (i in 1:k) {
      train_indices <- which(folds != i)
      test_indices <- which(folds == i)
      
      train_x <- x[train_indices, ]
      train_y <- y[train_indices]
      test_x <- x[test_indices, ]
      test_y <- y[test_indices]
      test_states <- states[test_indices]
      
      model <- FitMarkovSwitchingGAMLSS(
        x = as.matrix(train_x),
        y = train_y,
        init_state_trans_prob = init_probs,
        max.iter = max_iter,
        m.stop = mstop,
        type = type,
        formula = formula
      )
      
      predicted_states <- PredictStateSequence(model, test_y)
      accuracy[i] <- mean(predicted_states == test_states)
      
      # Compute likelihood for the test data
      test_data <- data.frame(y = test_y, x = test_x)
      likelihood[i] <- compute_likelihood(test_data, model$delta, model$gamma, model)
      
      model_results[[i]] <- model
      final_probs[[i]] <- model$gamma
      final_iterations[i] <- model$final_i
      
      current_iteration <- current_iteration + 1
      cat("Progress: ", round((current_iteration / total_iterations) * 100, 2), "%\n")
    }
    
    model_results_list[[toString(init_probs)]] <- model_results
    final_probs_list[[toString(init_probs)]] <- final_probs
    final_iterations_list[[toString(init_probs)]] <- final_iterations
    accuracy_list[[toString(init_probs)]] <- accuracy
    likelihood_list[[toString(init_probs)]] <- likelihood  # Store likelihoods
  }
  
  return(list(
    model_results = model_results_list,
    final_probs = final_probs_list,
    final_iterations = final_iterations_list,
    accuracy = accuracy_list,
    likelihood = likelihood_list  # Return likelihoods
  ))
}

cv_results <- function(results, init_state_probs_list, true_gamma, k) {
    gamma_list <- matrix(NA, length(init_state_probs_list) * k, 5)
    i <- 1
    while (i <= length(init_state_probs_list) * k) {
        for (model in 1:length(init_state_probs_list)) {
            for (j in 1:k) {
                gamma_list[i, 1] <- results$model_results[[model]][[j]]$gamma[1, 1]
                gamma_list[i, 2] <- results$model_results[[model]][[j]]$gamma[1, 2]
                gamma_list[i, 3] <- results$model_results[[model]][[j]]$gamma[2, 1]
                gamma_list[i, 4] <- results$model_results[[model]][[j]]$gamma[2, 2]
                gamma_list[i, 5] <- (abs(gamma_list[i, 1] - true_gamma[1, 1]) + abs(gamma_list[i, 2] - true_gamma[1, 2]) + abs(gamma_list[i, 3] - true_gamma[2, 1]) + abs(gamma_list[i, 4] - true_gamma[2, 2])) / 4
                i <- i + 1
            }
        }
    }
    colnames(gamma_list) <- c("gamma_11", "gamma_12", "gamma_21", "gamma_22", "abs_diff")
    gamma_df <- as.data.frame(gamma_list)
    gamma_df <- reshape2::melt(gamma_df)
    
    results_df <- data.frame(
        init_state_probs = rep(sapply(init_state_probs_list, function(x) paste(x, collapse = ",")), each = k),
        gamma_11 = gamma_list[, 1],
        gamma_12 = gamma_list[, 2],
        gamma_21 = gamma_list[, 3],
        gamma_22 = gamma_list[, 4],
        abs_dif = gamma_list[, 5],
        accuracy = unlist(results$accuracy),
        final_i = unlist(results$final_iterations),
        likelihood = unlist(results$likelihood)  # Add likelihood to results_df
    )
    
    agg_results_df <- results_df %>%
        group_by(init_state_probs) %>%
        summarise(
            mean_accuracy = mean(accuracy),
            mean_final_i = round(mean(final_i)),
            mean_abs_diff = mean(abs_dif),
            mean_likelihood = mean(likelihood)  # Add mean likelihood to aggregated results
        )
    
    list(gamma_df = gamma_df, results_df = results_df, agg_results_df = agg_results_df)
}