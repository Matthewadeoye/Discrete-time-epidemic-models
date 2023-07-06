# INCORPORATING ABC-REJECTION INTO ABC-MCMC TO GUESS INITIAL PARAMETER VALUES AND TOLERANCE

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
posterior.samples <- matrix(0, nrow = num_iterations, ncol = 2)

# Define a function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}


# Define prior distributions
prior_dist <- function(n) {
  beta_samples <- runif(n, 0, 0.05)  # Uniform prior between 0 and 0.05
  gamma_samples <- runif(n, 0, 0.5)     # Uniform prior between 0 and 0.5
  return(cbind(beta_samples, gamma_samples))
}


grid<- matrix(0, nrow = num_iterations, ncol = 3)

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
    posterior.samples[i, ] <- proposed_params
    grid[i,]<- c(posterior.samples[i,1], posterior.samples[i,2], metric)
  }
}

# Remove rows with all zeros and NA values from the parameter samples
posterior.samples <- posterior.samples[apply(posterior.samples, 1, function(row) !all(row == 0)), ]
posterior.samples <- posterior.samples[complete.cases(posterior.samples), ]

# Remove rows with all zeros and NA values from grid
grid <- grid[apply(grid, 1, function(row) !all(row == 0)), ]
grid <- grid[complete.cases(grid), ]
colnames(grid) <- c("Beta","Gamma", "Distance metric")


infection.rate <- posterior.samples[, 1]
recovery.rate <- posterior.samples[, 2]

# Plot the approximate posterior samples
hist(infection.rate, xlab= "Infection rate", main = "Posterior samples of Beta for deterministic SIR model")
hist(recovery.rate, xlab= "Recovery rate", main = "Posterior samples of Gamma for deterministic SIR model")


#Posterior means
mean(infection.rate)
mean(recovery.rate)

# Display the samples
head(posterior.samples)
tail(posterior.samples)

# Find the parameters with minimum distance metric
column_index<- 3
smallest.distance <- min(grid[,column_index])
row.number <- which(grid[, column_index] == smallest.distance, arr.ind = TRUE)[1]
# Print the row index
print(smallest.distance)
print(row.number)
grid[row.number, 1]
grid[row.number, 2]
grid[row.number, 3]

#ABC-MCMC begins here

epsilon <- smallest.distance
#epsilon <- 20

# Define number of iterations
num_iterations <- 150000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)

# Initialize parameter values
chain[1,]<- c(grid[row.number, 1], grid[row.number, 2])

# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], chain[1,]))
  
  # M-H probability
  #proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = chain[1,], log = TRUE))
  #proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = chain[1,], log = TRUE))
  
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.05, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.5, log = TRUE)
  
  
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.05, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.5, log = TRUE) 
  
  
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
prior_samples <- prior_dist(num_iterations)

par(mfrow = c(1, 2)) 
hist(prior_samples[, 1], xlab = "Infection rate", main = "Prior distribution of Beta", col = "lightblue")
hist(infection.rate, xlab = "Infection rate", main = "Posterior distribution of Beta", col = "lightgreen")

par(mfrow = c(1, 2))
hist(prior_samples[, 2], xlab = "Recovery rate", main = "Prior distribution of Gamma", col = "lightblue")
hist(recovery.rate, xlab = "Recovery rate", main = "Posterior distribution of Gamma", col = "lightgreen")

#Obtaining the range of posterior samples
min(infection.rate)   
max(infection.rate)   
min(recovery.rate)    
max(recovery.rate)    

#create the contour plot
x<- seq(0.00003,0.0003, 0.00005)

y<- seq(0.02, 0.1, 0.005)

u <- as.matrix(expand.grid(x, y))

z<- matrix(apply(u, 1, function(v) distance(model(c(v[1], v[2])),Observed_data)),
           nrow = length(x)) 

filled.contour(x,y,z, col=rainbow(39, 1, rev = FALSE), 
               xlab="Infection rate (Beta)",
               ylab="Recovery rate (Gamma)",
               main="Distance metric between observed and simulated data")
