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
step_size <- 1   # Step size for time discretization

simulated_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Observed data on removals (I’m assuming a Gaussian noise on the simulated removals)
observed_data <- data.frame(
  Steps = simulated_data$Steps,
  Removals = simulated_data$R + rnorm(length(simulated_data$R), mean = 0, sd = 35) # Simulated removals with added noise
)

# The objective function (sum of squared differences)
objective_function <- function(params) {
  beta <- params[1]
  gamma <- params[2]
  
  # Simulate SIR model using current parameter values
  model_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  
  # Compute sum of squared differences between simulated removals and observed removals
  sum((model_data$R - observed_data$Removals)^2)
}

# Perform optimization to estimate model parameters
initial_params <- c(0.000001, 0.9)  # Initial parameter values
opt_result <- optim(par = initial_params, fn = objective_function, method='Nelder-Mead')
estimated_params <- opt_result$par


# Simulate SIR model using estimated parameters
predicted_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, estimated_params[1], estimated_params[2], step_size)

# Plot observed and simulated removals  
library(ggplot2)
ggplot() +
  geom_line(data = observed_data, aes(x = Steps, y = Removals, color = "Observed"), size = 1) +
  geom_line(data = predicted_data, aes(x = Steps, y = R, color = "Simulated"), size = 1) +
  labs(x = "Time", y = "Removed Individuals", title = "Observed vs Simulated Removals") +
  scale_color_manual(values = c("Observed" = "blue", "Simulated" = "red")) +
  theme_minimal()

estimated_params
data.frame(observed_data,predicted_data$R)
