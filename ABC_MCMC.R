# Approximate Bayesian Computation Markov chain Monte Carlo for Discrete-time Deterministic SIR model

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

# Compute observed data (O_t)
Observed_data <- DT_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)


# Define model
model <- function(params) {
  # Simulate data based on the parameters
  Simulated_data <- DT_model(N, S0, I0, minTime, maxTime, params[1], params[2], step_size)
  return(Simulated_data)
}


# Define acceptance threshold
epsilon <- 60

# Define number of iterations
num_iterations <- 150000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Generate arbitrary initial parameter values
chain[1, ] <- runif(2)

# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the prior distributions
  proposed_params <- rnorm(2, chain[i - 1, ], c(0.6, 0.5))
  
 
  # Generate synthetic data based on the proposed parameters
  synthetic_data <- model(proposed_params)
  
  # Calculate the distance metric between synthetic and observed data
  metric <- distance(synthetic_data, Observed_data)
  
  # Compute based on the distance
  if (!is.na(metric) && metric <= epsilon) {
    likelihood <- 1
  } else {
    likelihood <- 0
  }
  
  # M-H probability
  proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = c(0.6, 0.5), log = TRUE))
  proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = c(0.6, 0.5), log = TRUE))
  
  prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.3, log = TRUE) +
    dnorm(proposed_params[2], mean = 0.04, sd = 0.3, log = TRUE)
  
  prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.3, log = TRUE) +
    dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.3, log = TRUE)
  
  mh.prob <- exp(proposal_proposed - proposal_current + prior_proposed - prior_current) * likelihood
  
  # Accept or reject the proposal
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    chain[i, ] <- proposed_params
  } else {
    chain[i, ] <- chain[i - 1, ]
  }
}

# Plot the ABC-MCMC chain
plot(chain[, 1], type = "l", main = "Chain of Parameter 1")
plot(chain[, 2], type = "l", main = "Chain of Parameter 2")

head(chain)
tail(chain)

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)


###########################################################################################################

# Approximate Bayesian Computation Markov chain Monte Carlo for Discrete-time Stochastic SIR model

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

# Define acceptance threshold
epsilon <- 60

# Define number of iterations
num_iterations <- 150000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Generate arbitrary initial parameter values
chain[1, ] <- runif(2)

# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the prior distributions
  proposed_params <- rnorm(2, chain[i - 1, ], c(0.6, 0.5))
  
  
  # Generate synthetic data based on the proposed parameters
  synthetic_data <- model(proposed_params)
  
  # Calculate the distance metric between synthetic and observed data
  metric <- distance(synthetic_data, Observed_data)
  
  # Compute based on the distance
  if (!is.na(metric) && metric <= epsilon) {
    likelihood <- 1
  } else {
    likelihood <- 0
  }
  
  # M-H probability
  proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = c(0.6, 0.5), log = TRUE))
  proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = c(0.6, 0.5), log = TRUE))
  
  prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.3, log = TRUE) +
    dnorm(proposed_params[2], mean = 0.04, sd = 0.3, log = TRUE)
  
  prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.3, log = TRUE) +
    dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.3, log = TRUE)
  
  mh.prob <- exp(proposal_proposed - proposal_current + prior_proposed - prior_current) * likelihood
  
  # Accept or reject the proposal
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    chain[i, ] <- proposed_params
  } else {
    chain[i, ] <- chain[i - 1, ]
  }
}

# Plot chain
plot(chain[, 1], type = "l", main = "Chain of Parameter 1")
plot(chain[, 2], type = "l", main = "Chain of Parameter 2")

head(chain)
tail(chain)

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)
