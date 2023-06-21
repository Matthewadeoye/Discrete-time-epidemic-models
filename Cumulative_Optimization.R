set.seed(8)

# Define the Discrete-Time Deterministic SIR Model
Deterministic_DT_SIR_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  steps <- seq(minTime, maxTime, by = step_size) # Time discretization
  S <- numeric(length(steps)) # Vector to store susceptibles
  I <- numeric(length(steps)) # Vector to store infectives
  R <- numeric(length(steps)) # Vector to store removed
  
  # Set initial conditions
  S[1] <- S0
  I[1] <- I0
  R[1] <- N - S0 - I0
  
  # Loop and update compartments
  for (t in 2:length(steps)) {
    S[t] <- S[t-1] - (beta * S[t-1] * I[t-1]) * step_size
    I[t] <- I[t-1] + (beta * S[t-1] * I[t-1] - gamma * I[t-1]) * step_size
    R[t] <- R[t-1] + gamma * I[t-1] * step_size
  }
  
  # Return simulated trajectories
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Generate simulated data using the SIR model
N <- 1000   # Total population size
S0 <- 900   # Initial number of susceptibles
I0 <- 100   # Initial number of infectives
minTime <- 1
maxTime <- 365
beta <- 0.0001   # Infection rate
gamma <- 0.05    # Recovery rate
step_size <- 0.01   # Step size for time discretization

simulated_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Observed data on removals (cumulative with Binomial or Poisson noise)
observed_data <- data.frame(
  Steps = simulated_data$Steps,
  Removals= c(0,cumsum(rpois(length(simulated_data$Steps)-1, diff(simulated_data$R))))
  #Removals = cumsum(rbinom(length(simulated_data$R), size = 1, prob = simulated_data$R / N)) # Cumulative observed removals using Binomial 
  #Removals = cumsum(rpois(length(simulated_data$R), lambda = simulated_data$R)) # Cumulative observed removals with Poisson noise
)


# The objective function (sum of squared differences)
objective_function <- function(params) {
  beta <- params[1]
  gamma <- params[2]
  
  beta<- exp(beta)
  gamma<- exp(gamma)
  # Simulate SIR model using current parameter values
  model_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  
  # Compute sum of squared differences between simulated removals and observed removals
  OLS = sum((model_data$R - observed_data$Removals)^2)
  
  #Alternatively using MLE
  #NLL= -sum(log(dnorm(observed_data$Removals, mean=model_data$R, sd=0.1*mean(observed_data$Removals)))) #Assuming a Gaussian error term
  #NLL2 = -sum(log(dpois(round(observed_data$Removals/10000),round(model_data$R/10000)))) #Assuming observed Removals follows a Poisson distribution
  return(OLS)
}

# Perform optimization to estimate model parameters
initial_params <- c(log(0.0000001), log(0.5))  # Initial parameter values
#Optim
#opt_result <- optim(par = initial_params, fn = objective_function, method='Nelder-Mead')
#print(opt_result)
#estimated_params <- opt_result$par
#estimated_params<- exp(estimated_params)

#nlm
opt_result <- nlm(p = initial_params, f = objective_function)
print(opt_result)
estimated_params <- opt_result$estimate
estimated_params<- exp(estimated_params)


# Simulate SIR model using estimated parameters
predicted_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, estimated_params[1], estimated_params[2], step_size)
Initpredicted_data<- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, exp(initial_params[1]), exp(initial_params[2]), step_size)

# Plot observed and simulated removals  
library(ggplot2)
ggplot() +
  geom_point(data = observed_data, aes(x = Steps, y = Removals, color = "Observed"), size = 1) +
  geom_line(data = predicted_data, aes(x = Steps, y = R, color = "Predicted at optimized values"), size = 1) +
  geom_line(data = Initpredicted_data, aes(x = Steps, y =R, color = "Predicted at initial values"), size = 1) +
  labs(x = "Time", y = "Cumulative removals", title = "Observed vs predicted cumulative removals") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted at optimized values" = "red", "Predicted at initial values" = "black")) +
  theme_minimal()


estimated_params

data.frame(observed_data,predicted_data$R, Initpredicted_data$R)



#Manual Search
# Define the objective function
objective_function <- function(params) {
  beta <- params[1]
  gamma <- params[2]
  
  # Simulate SIR model using current parameter values
  model_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  
  # Compute sum of squared differences between simulated removals and observed removals
  sum((model_data$R - observed_data$Removals)^2)
}

# Create a grid of parameter values
beta_values <- seq(0.000001, 0.001, length.out = 50)
gamma_values <- seq(0.01, 0.095, length.out = 50)
grid <- expand.grid(beta = beta_values, gamma = gamma_values)

# Evaluate the objective function for each parameter combination
z <- apply(grid, 1, FUN = objective_function)

# Add parameter values and objective function values to the grid
grid$z <- z

# Plot the objective function
library(ggplot2)

ggplot(data = grid) +
  geom_contour(aes(x = beta, y = gamma, z = z), bins = 20) +
  labs(x = "Beta", y = "Gamma", z = "Objective Function") +
  theme_minimal()

# Find the row with minimum sum of square difference
row_index <- which(grid$z == min(grid$z))
# Print the row index
print(row_index)

