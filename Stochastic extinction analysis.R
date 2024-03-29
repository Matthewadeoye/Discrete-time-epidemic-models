#S-S
#Stochastic SIR #6, #8, #4, #3, #10, #222(works better)
start_time <- Sys.time()
set.seed(222)

# Initial conditions and parameter values
N <- 500    # Total population size
S0 <- 499    # Initial number of susceptibles
I0 <- 1    # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.1    # Infection rate
gamma <- 0.05     # Recovery rate
step_size <- 1    # Step size for time discretization

# Define the Discrete-Time Stochastic SIR Model
ST_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  minTime <- minTime           # start time
  maxTime <- maxTime           # end time
  step_size <- step_size       # step size
  Steps <- seq(minTime, maxTime, by = step_size)   # time discretization
  beta <- beta                 # infection rate
  gamma <- gamma               # recovery rate
  S <- numeric(length(Steps))  # create empty vector to store simulations for susceptibles
  I <- numeric(length(Steps))  # create empty vector to store simulations for infectives
  R <- numeric(length(Steps))  # create empty vector to store simulations for removed
  S[1] <- S0                # initial proportion of susceptibles
  I[1] <- I0                # initial proportion of infectives
  R[1] <- N - S0 - I0       # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(Steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-(beta/N) * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1]), p_I)
    delta_I <- rbinom(1, round(I[tym - 1]), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S 
    I[tym] <- I[tym - 1] + delta_S  - delta_I 
    R[tym] <- R[tym - 1] + delta_I 
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
epsilon <- 25

# Define the number of iterations
num_iterations <- 70000

# Saving parameter samples
posterior.samples <- matrix(0, nrow = num_iterations, ncol = 2)

# Define a function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}


# Define prior distributions
prior_dist <- function(n) {
  beta_samples <- runif(n, 0, 0.2)  # Uniform prior between 0 and 0.2
  gamma_samples <- runif(n, 0, 0.1)     # Uniform prior between 0 and 0.1
  return(cbind(beta_samples, gamma_samples))
}

# Alternative prior distributions
#prior_dist <- function(n, params) {
#  beta_samples <- abs(rnorm(n, mean = 0.0003, sd = 0.03))  # Gaussian prior for beta
#  gamma_samples <- abs(rnorm(n, mean = 0.04, sd = 0.03))   # Gaussian prior for gamma
#  return(cbind(beta_samples, gamma_samples))
#}


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

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time

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


##ABC-MCMC begins here

#epsilon <- smallest.distance
epsilon <- 1.60

# Define number of iterations
num_iterations <- 250000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)


# Initialize parameter values
chain[1,]<- c(grid[row.number, 1], grid[row.number, 2])

# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], c(0.03,0.01)))
  
  # M-H probability
  #proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = chain[1,], log = TRUE))
  #proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = chain[1,], log = TRUE))
  
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.2, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.2, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  
  #If using the alternative prior distributions
  # prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.03, log = TRUE) +
  #  dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.03, log = TRUE)
  #prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.03, log = TRUE) +
  # dnorm(proposed_params[2], mean = 0.04, sd = 0.03, log = TRUE)
  
  
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

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time


#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area
#par(mar=c(4,2,2,2))
hist(infection.rate, freq=F, xlab = expression(beta), main =bquote("Marginal posterior for " * bold(beta)), col = "white", ylab="Density")
abline(v=0.1, col="red", lwd=2,lty=1)
abline(h=5, col="blue", lwd=2,lty=1)
plot(infection.rate, type = "l", ylab=expression(beta), main =bquote("Trace plot for " * bold(beta)), xlab="MCMC iteration", ylim=c(0,0.3))
#mtext("(a)", side = 3, line = -26, outer = TRUE)
plot(infection.rate,recovery.rate,  pch = 16, col = "#00000005", main="Joint posterior", xlab=expression(beta),ylab=expression(gamma))
hist(recovery.rate, freq=F, xlab = expression(gamma), main =bquote("Marginal posterior for " * bold(gamma)), col = "white",ylab="Density")
abline(v=0.05, col="red", lwd=2,lty=1)
abline(h=10, col="blue", lwd=2,lty=1)
plot(recovery.rate, type = "l", ylab=expression(gamma), main =bquote("Trace plot for " * bold(gamma)), xlab="MCMC iteration",ylim=c(0,0.3))
#mtext("(b)", side = 1, line = -1, outer = TRUE)
plot(recovery.rate,infection.rate,  pch = 16, col = "#00000005", main="Joint posterior", xlab=expression(gamma),ylab=expression(beta))
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

add_legend("topright", legend=c("True value", "Prior density"), lty=1, 
           col=c("red", "blue"),
           horiz=TRUE, bty='n', cex=1.1)

#df<- data.frame(Beta=infection.rate,Gamma=recovery.rate)
#ggplot(df, aes(Beta, Gamma)) + stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + scale_fill_viridis_c(name = "density") + geom_point(shape = '.') + guides(fill = FALSE)  # This line removes the fill legend


library(bayestestR)
mean(infection.rate)
ci(infection.rate, method = "HDI")
ci(infection.rate, method = "ETI")

mean(recovery.rate)
ci(recovery.rate, method = "HDI")
ci(recovery.rate, method = "ETI")


# Create a kernel density estimate of the posterior distribution
posterior_density1 <- density(infection.rate)

# Find the mode using numerical optimization
mode_estimate1 <- posterior_density1$x[which.max(posterior_density1$y)]

print(mode_estimate1)

# Create a kernel density estimate of the posterior distribution
posterior_density2 <- density(recovery.rate)

# Find the mode using numerical optimization
mode_estimate2 <- posterior_density2$x[which.max(posterior_density2$y)]

print(mode_estimate2)

SS.R0<- infection.rate/recovery.rate
hist(SS.R0, freq=F, xlab =expression("R"[0]), breaks=37, xlim=c(0,10), main = "(a)", col = "white",ylab="Density")
abline(v=2, col="red", lwd=2,lty=1)


################################################################################################################################

#S-S
#Stochastic SIR #6, #8, #4, #3, #10, #222(works better)
start_time <- Sys.time()
set.seed(222)

# Initial conditions and parameter values
N <- 500    # Total population size
S0 <- 499    # Initial number of susceptibles
I0 <- 1    # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.1    # Infection rate
gamma <- 0.05     # Recovery rate
step_size <- 1    # Step size for time discretization

# Define the Discrete-Time Stochastic SIR Model
ST_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  minTime <- minTime           # start time
  maxTime <- maxTime           # end time
  step_size <- step_size       # step size
  Steps <- seq(minTime, maxTime, by = step_size)   # time discretization
  beta <- beta                 # infection rate
  gamma <- gamma               # recovery rate
  S <- numeric(length(Steps))  # create empty vector to store simulations for susceptibles
  I <- numeric(length(Steps))  # create empty vector to store simulations for infectives
  R <- numeric(length(Steps))  # create empty vector to store simulations for removed
  S[1] <- S0                # initial proportion of susceptibles
  I[1] <- I0                # initial proportion of infectives
  R[1] <- N - S0 - I0       # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(Steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-(beta/N) * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1]), p_I)
    delta_I <- rbinom(1, round(I[tym - 1]), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S 
    I[tym] <- I[tym - 1] + delta_S  - delta_I 
    R[tym] <- R[tym - 1] + delta_I 
  }
  
  # Actual data on removals (O_t)
  data <- data.frame(Steps = Steps, Removals = c(0, rpois(length(Steps) - 1, diff(R))))
  return(data$Removals)
}

# Define the observed data (O_t)
Observed_data <- ST_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

#S-D
#Stochastic SIR #6, #8, #4, #3, #10, #222(works better)
start_time <- Sys.time()
set.seed(222)

# Initial conditions and parameter values
N <- 500    # Total population size
S0 <- 499    # Initial number of susceptibles
I0 <- 1    # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.1    # Infection rate
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
    S[t] <- S[t-1] - ((beta/N) * S[t-1] * I[t-1]) * step_size
    I[t] <- I[t-1] + ((beta/N) * S[t-1] * I[t-1] - gamma * I[t-1]) * step_size
    R[t] <- R[t-1] + gamma * I[t-1] * step_size
  }
  
  # Actual data on removals (O_t)
  data <- data.frame(Steps = Steps, Removals = c(0, rpois(length(Steps) - 1, diff(R))))
  return(data$Removals)
}

# Define the observed data (O_t)
#Observed_data <- DT_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Define model
model <- function(params) {
  # Simulate data based on the parameters
  Simulated_data <- DT_model(N, S0, I0, minTime, maxTime, params[1], params[2], step_size)
  return(Simulated_data)
}


# Define the acceptance threshold
epsilon <- 25

# Define the number of iterations
num_iterations <- 70000

# Saving parameter samples
posterior.samples <- matrix(0, nrow = num_iterations, ncol = 2)

# Define a function to calculate the distance metric
distance <- function(x, y) {
  return(sqrt(sum(((x - y)^2) / (y + 1))))
}


# Define prior distributions
prior_dist <- function(n) {
  beta_samples <- runif(n, 0, 0.2)  # Uniform prior between 0 and 0.2
  gamma_samples <- runif(n, 0, 0.1)     # Uniform prior between 0 and 0.1
  return(cbind(beta_samples, gamma_samples))
}

# Alternative prior distributions
#prior_dist <- function(n, params) {
#  beta_samples <- abs(rnorm(n, mean = 0.0003, sd = 0.03))  # Gaussian prior for beta
#  gamma_samples <- abs(rnorm(n, mean = 0.04, sd = 0.03))   # Gaussian prior for gamma
#  return(cbind(beta_samples, gamma_samples))
#}


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

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time

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


##ABC-MCMC begins here

#epsilon <- smallest.distance
epsilon <- 1.60

# Define number of iterations
num_iterations <- 250000

# Initialize the MCMC chain
chain <- matrix(0, nrow = num_iterations, ncol = 2)


# Initialize parameter values
chain[1,]<- c(grid[row.number, 1], grid[row.number, 2])

# ABC-MCMC algorithm
for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], c(0.009,0.01)))
  
  # M-H probability
  #proposal_proposed <- sum(dnorm(proposed_params, mean = chain[i - 1, ], sd = chain[1,], log = TRUE))
  #proposal_current <- sum(dnorm(chain[i - 1, ], mean = proposed_params, sd = chain[1,], log = TRUE))
  
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.2, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.2, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  
  #If using the alternative prior distributions
  # prior_current <- dnorm(chain[i - 1, 1], mean = 0.0003, sd = 0.03, log = TRUE) +
  #  dnorm(chain[i - 1, 2], mean = 0.04, sd = 0.03, log = TRUE)
  #prior_proposed <- dnorm(proposed_params[1], mean = 0.0003, sd = 0.03, log = TRUE) +
  # dnorm(proposed_params[2], mean = 0.04, sd = 0.03, log = TRUE)
  
  
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

infection.rate<- chain[,1]
recovery.rate<- chain[,2]

#Posterior means
mean(infection.rate)
mean(recovery.rate)

end_time <- Sys.time()
end_time - start_time


#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area
#par(mar=c(4,2,2,2))
hist(infection.rate, freq=F, xlab = expression(beta), main =bquote("Marginal posterior for " * bold(beta)), col = "white", ylab="Density")
abline(v=0.1, col="red", lwd=2,lty=1)
abline(h=5, col="blue", lwd=2,lty=1)
plot(infection.rate, type = "l", ylab=expression(beta), main =bquote("Trace plot for " * bold(beta)), xlab="MCMC iteration", ylim=c(0,0.3))
#mtext("(a)", side = 3, line = -26, outer = TRUE)
plot(infection.rate,recovery.rate,  pch = 16, col = "#00000005", main="Joint posterior", xlab=expression(beta),ylab=expression(gamma))
hist(recovery.rate, freq=F, xlab = expression(gamma), main =bquote("Marginal posterior for " * bold(gamma)), col = "white",ylab="Density")
abline(v=0.05, col="red", lwd=2,lty=1)
abline(h=10, col="blue", lwd=2,lty=1)
plot(recovery.rate, type = "l", ylab=expression(gamma), main =bquote("Trace plot for " * bold(gamma)), xlab="MCMC iteration",ylim=c(0,0.3))
#mtext("(b)", side = 1, line = -1, outer = TRUE)
plot(recovery.rate,infection.rate,  pch = 16, col = "#00000005", main="Joint posterior", xlab=expression(gamma),ylab=expression(beta))
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

add_legend("topright", legend=c("True value", "Prior density"), lty=1, 
           col=c("red", "blue"),
           horiz=TRUE, bty='n', cex=1.1)

#df<- data.frame(Beta=infection.rate,Gamma=recovery.rate)
#ggplot(df, aes(Beta, Gamma)) + stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + scale_fill_viridis_c(name = "density") + geom_point(shape = '.') + guides(fill = FALSE)  # This line removes the fill legend


# Compute the lower and upper quantiles for the credible interval
lower_quantile <- quantile(infection.rate, probs = 0.025)  # 2.5th percentile
upper_quantile <- quantile(infection.rate, probs = 0.975)  # 97.5th percentile

# Compute the 95% credible interval
Infectioncredible_interval <- c(lower_quantile, upper_quantile)

print(Infectioncredible_interval)

# Compute the lower and upper quantiles for the credible interval
lower_quantile <- quantile(recovery.rate, probs = 0.025)  # 2.5th percentile
upper_quantile <- quantile(recovery.rate, probs = 0.975)  # 97.5th percentile

# Compute the 95% credible interval
Recoverycredible_interval <- c(lower_quantile, upper_quantile)

print(Recoverycredible_interval)

library(bayestestR)
mean(infection.rate)
ci(infection.rate, method = "HDI")
ci(infection.rate, method = "ETI")

mean(recovery.rate)
ci(recovery.rate, method = "HDI")
ci(recovery.rate, method = "ETI")


# Create a kernel density estimate of the posterior distribution
posterior_density1 <- density(infection.rate)

# Find the mode using numerical optimization
mode_estimate1 <- posterior_density1$x[which.max(posterior_density1$y)]

print(mode_estimate1)

# Create a kernel density estimate of the posterior distribution
posterior_density2 <- density(recovery.rate)

# Find the mode using numerical optimization
mode_estimate2 <- posterior_density2$x[which.max(posterior_density2$y)]

print(mode_estimate2)

SD.R0<- infection.rate/recovery.rate
hist(SD.R0, freq=F, xlab =expression("R"[0]), breaks=3, xlim=c(0,10), main = "(a)", col = "white",ylab="Density")
abline(v=2, col="red", lwd=2,lty=1)
