# Approximate Bayesian Computation Markov chain Monte Carlo for Discrete-time Deterministic SIR model

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

# Define acceptance threshold
epsilon <- 20

# Define number of iterations
num_iterations <- 200000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Initialize parameter values
chain[1,]<- c(runif(1, 0, 0.001), runif(1, 0, 0.1))


# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], chain[1,]))
  
  # M-H probability
  #proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = c(0.6, 0.5), log = TRUE))
  #proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = c(0.6, 0.5), log = TRUE))
  
  
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.001, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.001, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  #If using the alternative prior distributions
  #prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.03, log = TRUE) +
  # dnorm(proposed_params[2], mean = 0.04, sd = 0.03, log = TRUE)
  # prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.03, log = TRUE) +
  #  dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.03, log = TRUE)
  
  
  
  #mh.prob <- exp(proposal_proposed - proposal_current + prior_proposed - prior_current)
  
   mh.prob <- exp(prior_proposed - prior_current)
  
  # Accept (and simulate) or reject the proposals
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    
    # Generate synthetic data based on the proposed parameters
    synthetic_data <- model(proposed_params)
    
    # Calculate the distance metric between synthetic and observed data
    metric <- distance(synthetic_data, Observed_data)
    
    # Evaluate distance metric
    if (!is.na(metric) && metric <= epsilon) {
      chain[i, ] <- proposed_params 
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  } else {
    chain[i, ] <- chain[i - 1, ]  # Keep the current parameter values in the chain
  }
}


# Plot the ABC-MCMC chain
plot(chain[, 1], type = "l", ylab="Beta", main = "Chain of Beta for Deterministic SIR model")
plot(chain[, 2], type = "l", ylab="Gamma", main = "Chain of Gamma for Deterministic SIR model")

head(chain)
tail(chain)

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time


# Plot prior and posterior distributions
#prior_samples <- prior_dist(num_iterations)

#par(mfrow = c(1, 2)) 
#hist(prior_samples[, 1], xlab = "Infection rate", main = "Prior distribution of Beta", col = "lightblue")
hist(infection.rate, xlab = "Infection rate", main = "Posterior distribution of Beta", col = "lightgreen")

#par(mfrow = c(1, 2))
#hist(prior_samples[, 2], xlab = "Recovery rate", main = "Prior distribution of Gamma", col = "lightblue")
hist(recovery.rate, xlab = "Recovery rate", main = "Posterior distribution of Gamma", col = "lightgreen")


###########################################################################################################


# Approximate Bayesian Computation Markov chain Monte Carlo for Discrete-time Stochastic SIR model

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

# Define acceptance threshold
epsilon <- 20

# Define number of iterations
num_iterations <- 200000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Define function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}

# Initialize parameter values
chain[1,]<- c(runif(1, 0, 0.001), runif(1, 0, 0.1))


# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], chain[1,]))
  
  # M-H probability
  #proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = c(0.6, 0.5), log = TRUE))
  #proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = c(0.6, 0.5), log = TRUE))
  
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.001, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.001, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  
  #If using the alternative prior distributions
  #prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.03, log = TRUE) +
  # dnorm(proposed_params[2], mean = 0.04, sd = 0.03, log = TRUE)
  # prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.03, log = TRUE) +
  #  dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.03, log = TRUE)
  
  
  
  #mh.prob <- exp(proposal_proposed - proposal_current + prior_proposed - prior_current)
  
   mh.prob <- exp(prior_proposed - prior_current)
  
  # Accept (and simulate) or reject the proposals
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    
    # Generate synthetic data based on the proposed parameters
    synthetic_data <- model(proposed_params)
    
    # Calculate the distance metric between synthetic and observed data
    metric <- distance(synthetic_data, Observed_data)
    
    # Evaluate distance metric
    if (!is.na(metric) && metric <= epsilon) {
      chain[i, ] <- proposed_params 
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  } else {
    chain[i, ] <- chain[i - 1, ]  # Keep the current parameter values in the chain
  }
}


# Plot the ABC-MCMC chain
plot(chain[, 1], type = "l", ylab="Beta", main = "Chain of Beta for Stochastic SIR model")
plot(chain[, 2], type = "l", ylab="Gamma", main = "Chain of Gamma for Stochastic SIR model")

head(chain)
tail(chain)

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time


# Plot prior and posterior distributions
#prior_samples <- prior_dist(num_iterations)

#par(mfrow = c(1, 2)) 
#hist(prior_samples[, 1], xlab = "Infection rate", main = "Prior distribution of Beta", col = "lightblue")
hist(infection.rate, xlab = "Infection rate", main = "Posterior distribution of Beta", col = "lightgreen")

#par(mfrow = c(1, 2))
#hist(prior_samples[, 2], xlab = "Recovery rate", main = "Prior distribution of Gamma", col = "lightblue")
hist(recovery.rate, xlab = "Recovery rate", main = "Posterior distribution of Gamma", col = "lightgreen")

