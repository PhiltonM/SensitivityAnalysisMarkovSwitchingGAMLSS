# Linear setting
set.seed(1896)
n <- 500
states <- 2
p <- 100

# Transition probability matrix
gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = states, byrow = TRUE)

# Initial state probabilities
delta <- c(0.5, 0.5)

# Simulate state sequence
state_seq <- numeric(n)
state_seq[1] <- sample(1:states, 1, prob = delta)
for (t in 2:n) {
  state_seq[t] <- sample(1:states, 1, prob = gamma[state_seq[t - 1], ])
}

# Simulate explanatory variables
x <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)
mu_1 <- exp(2 + 2 * x[t, 1])
mu_2 <- exp(2 - 2 * x[t, 1])
size_1 <- exp(2 * x[t, 1])
size_2 <- exp(-2 * x[t, 1])

# Simulate response variable
y <- numeric(n)
for (t in 1:n) {
  if (state_seq[t] == 1) {
    y[t] <- rnbinom(1, mu = mu_1, size = size_1)
  } else {
    y[t] <- rnbinom(1, mu = mu_2, size = size_2)
  }
}
# Save x and y from the linear setting in a dataframe
linear_data <- data.frame(y = y, x)
linear_data$states <- state_seq


# Non-linear setting
set.seed(1896)
n <- 500
states <- 2
p <- 100

# Transition probability matrix
gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = states, byrow = TRUE)

# Initial state probabilities
delta <- c(0.5, 0.5)

# Simulate state sequence
state_seq <- numeric(n)
state_seq[1] <- sample(1:states, 1, prob = delta)
for (t in 2:n) {
  state_seq[t] <- sample(1:states, 1, prob = gamma[state_seq[t - 1], ])
}

# Simulate explanatory variables
x <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)
mu_1 <- 2 + 2 * sin(pi * (x[t, 1] - 0.5))
mu_2 <- -2 - sin(pi * (x[t, 1] - 0.5))
size_1 <- exp(sin(pi * (x[t, 1] - 0.5)))
size_2 <- exp(-2 * sin(pi * (x[t, 1] - 0.5)))

# Simulate response variable
y <- numeric(n)
for (t in 1:n) {
  if (state_seq[t] == 1) {
    y[t] <- rnorm(1,
                  mean = mu_1,
                  sd = size_1)
  } else {
    y[t] <- rnorm(1,
                  mean = mu_2,
                  sd = size_2)
  }
}
nonlinear_data <- data.frame(y = y, x)
nonlinear_data$states <- state_seq


