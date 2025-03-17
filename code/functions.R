
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
renv::activate()
library(fda)
library(boot)
library(numDeriv)
library(gamlss)
library(gamboostLSS)
library(ggplot2)
library(dplyr)
required_packages <- c("fda", "boot", "numDeriv", "gamlss", "gamboostLSS", "ggplot2", "dplyr")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}
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
            mod[[j]] = glmboostLSS(y ~. , weights = ind[j,], data = data.frame(x, y), families = as.families("NO"), method = "noncyclic", control = boost_control(mstop = m.stop[j], nu = 0.1))
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

# K-Fold Cross Validation

cross_validate_msgamlss <- function(x, y, init_probs_list, type = "MSGLMLSS", formula = NULL) {
  n <- length(y)
  k <- 2
  folds <- rep(1:2, length.out = n)
  folds[seq(1, n, by = 2)] <- 1
  folds[seq(2, n, by = 2)] <- 2
  mse_list <- list()
  accuracy_list <- list()
  model_results_list <- list()
  final_probs_list <- list()
  final_iterations_list <- list()
  
  total_iterations <- length(init_probs_list) * k
  current_iteration <- 0
  
  for (init_probs in init_probs_list) {
    mse <- numeric(k)
    accuracy <- numeric(k)
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
      
      model <- FitMarkovSwitchingGAMLSS(
        x = as.matrix(train_x),
        y = train_y,
        init_state_trans_prob = init_probs,
        max.iter = 500,
        formula = formula)
      


      model_results[[i]] <- model
      final_probs[[i]] <- model$gamma
      final_iterations[i] <- model$final_i
      
      current_iteration <- current_iteration + 1
      cat("Progress: ", round((current_iteration / total_iterations) * 100, 2), "%\n")
    }
    
    model_results_list[[toString(init_probs)]] <- model_results
    final_probs_list[[toString(init_probs)]] <- final_probs
    final_iterations_list[[toString(init_probs)]] <- final_iterations
  }
  
  return(list(

    model_results = model_results_list,
    final_probs = final_probs_list,
    final_iterations = final_iterations_list
  ))
}
