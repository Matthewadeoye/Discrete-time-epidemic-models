# Approximate Bayesian Computation Rejection Sampling for Discrete-time Deterministic SIR model

start_time <- Sys.time()
set.seed(8)

# Initial conditions and parameter values
N <- 1000    # Total population size
S0 <- 900    # Initial number of susceptibles
I0 <- 100    # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.0001    # Infection rate
gamma <- 0.05     # Recovery rate
step_size <- 1    # Step size for time discretization

# Define the Discrete-Time Deterministic SIR Model
DT_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  Steps <- seq(minTime, maxTime, by = step_size)   # Time discretization
  S <- numeric(length(Steps))   # Vector to store susceptibles
  I <- numeric(length(Steps))   # Vector to store infectives
  R <- numeric(length(Steps))   # Vector to store removed
  
  # Assign initial conditions
  S[1] <- S0
  I[1] <- I0
  R[1] <- N - S0 - I0
  
  # Loop and update compartments
  for (t in 2:length(Steps)) {
    S[t] <- S[t-1] - (beta * S[t-1] * I[t-1]) * step_size
    I[t] <- I[t-1] + (beta * S[t-1] * I[t-1] - gamma * I[t-1]) * step_size
    R[t] <- R[t-1] + gamma * I[t-1] * step_size
  }
  
  # Actual data on removals (O_t)
  data <- data.frame(Steps = Steps, Removals = c(0, rpois(length(Steps) - 1, diff(R))))
  return(data$Removals)
}

# Define the observed data (O_t)
Observed_data <- DT_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Define model
model <- function(params) {
  # Simulate data based on the parameters
  Simulated_data <- DT_model(N, S0, I0, minTime, maxTime, params[1], params[2], step_size)
  return(Simulated_data)
}

# Define the acceptance threshold
epsilon <- 20

# Define the number of iterations
num_iterations <- 200000

# Saving parameter samples
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define a function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Define prior distributions
prior_dist <- function(n) {
  beta_samples <- runif(n, 0, 0.001)  # Uniform prior between 0 and 0.001
  gamma_samples <- runif(n, 0, 0.1)     # Uniform prior between 0 and 0.1
  return(cbind(beta_samples, gamma_samples))
}


# ABC-rejection algorithm
for (i in 1:num_iterations) {
  # Generate parameter proposals from the prior distributions
  proposed_params <- prior_dist(1)
  
  # Generate synthetic data based on the proposed parameters
  synthetic_data <- model(proposed_params)
  
  # Calculate the distance metric between synthetic and observed
  metric <- distance(synthetic_data, Observed_data)
  
  # Accept or reject the proposed parameters based on the distance
  if (!is.na(metric) && metric <= epsilon) {
    chain[i, ] <- proposed_params
  }
}

# Remove rows with all zeros and NA values from the parameter samples
chain <- chain[apply(chain, 1, function(row) !all(row == 0)), ]
chain <- chain[complete.cases(chain), ]

# Plot the approximate posterior samples
plot(chain[, 1], type = "l", main = "Chain of Parameter 1")
plot(chain[, 2], type = "l", main = "Chain of Parameter 2")

infection.rate <- chain[, 1]
recovery.rate <- chain[, 2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

# Display the samples
head(chain)
tail(chain)

end_time <- Sys.time()
end_time - start_time


#########################################################################################

# Approximate Bayesian Computation Rejection Sampling for Discrete-time Stochastic SIR model

start_time <- Sys.time()
set.seed(8)

# Initial conditions and parameter values
N <- 1000    # Total population size
S0 <- 900    # Initial number of susceptibles
I0 <- 100    # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.0001    # Infection rate
gamma <- 0.05     # Recovery rate
step_size <- 1    # Step size for time discretization

# Define the Discrete-Time Stochastic SIR Model
ST_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  Steps <- seq(minTime, maxTime, by = step_size)   # Time discretization
  S <- numeric(length(Steps))   # Vector to store susceptibles
  I <- numeric(length(Steps))   # Vector to store infectives
  R <- numeric(length(Steps))   # Vector to store removed
  
  # Assign initial conditions
  S[1] <- S0
  I[1] <- I0
  R[1] <- N - S0 - I0
  
  # Loop through discretized time and update the compartments at each step
  for (t in 2:length(Steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-(beta) * I[t - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, S[t - 1], p_I)
    delta_I <- rbinom(1, I[t - 1], p_R)
    
    # Update the compartments
    S[t] <- S[t - 1] - delta_S
    I[t] <- I[t - 1] + delta_S - delta_I
    R[t] <- R[t - 1] + delta_I
  }
  
  # Actual data on removals (O_t)
  data <- data.frame(Steps = Steps, Removals = c(0, rpois(length(Steps) - 1, diff(R))))
  return(data$Removals)
}

# Define the observed data (O_t)
Observed_data <- ST_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Define model
model <- function(params) {
  # Simulate data based on the parameters
  Simulated_data <- ST_model(N, S0, I0, minTime, maxTime, params[1], params[2], step_size)
  return(Simulated_data)
}

# Define the acceptance threshold
epsilon <- 20

# Define the number of iterations
num_iterations <- 200000

# Saving parameter samples
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define a function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Define prior distributions
prior_dist <- function(n) {
  beta_samples <- runif(n, 0, 0.001)  # Uniform prior between 0 and 0.001
  gamma_samples <- runif(n, 0, 0.1)     # Uniform prior between 0 and 0.1
  return(cbind(beta_samples, gamma_samples))
}


# ABC-rejection algorithm
for (i in 1:num_iterations) {
  # Generate parameter proposals from the prior distributions
  proposed_params <- prior_dist(1)
  
  # Generate synthetic data based on the proposed parameters
  synthetic_data <- model(proposed_params)
  
  # Calculate the distance metric between synthetic and observed
  metric <- distance(synthetic_data, Observed_data)
  
  # Accept or reject the proposed parameters based on the distance
  if (!is.na(metric) && metric <= epsilon) {
    chain[i, ] <- proposed_params
  }
}

# Remove rows with zeros and NA values
chain <- chain[apply(chain, 1, function(row) !all(row == 0)), ]
chain <- chain[complete.cases(chain), ]

# Plot the chains
plot(chain[, 1], type = "l", main = "Chain of Parameter 1")
plot(chain[, 2], type = "l", main = "Chain of Parameter 2")

infection.rate <- chain[, 1]
recovery.rate <- chain[, 2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

# Display the samples
head(chain)
tail(chain)

end_time <- Sys.time()
end_time - start_time