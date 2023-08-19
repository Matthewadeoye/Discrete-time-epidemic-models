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
    S[t] <- S[t-1] - ((beta/N) * S[t-1] * I[t-1]) * step_size
    I[t] <- I[t-1] + ((beta/N) * S[t-1] * I[t-1] - gamma * I[t-1]) * step_size
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
beta <- 0.1   # Infection rate
gamma <- 0.05    # Recovery rate
step_size <- 1   # Step size for time discretization

simulated_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)

# Observed data on removals (cumulative with Poisson noise)
observed_data <- data.frame(
  Steps = simulated_data$Steps,
  Removals= c(0, rpois(length(simulated_data$Steps)-1, diff(simulated_data$R))))



# The objective function (sum of squared differences)
objective_function <- function(params) {
  beta <- params[1]
  gamma <- params[2]
  
  beta<- exp(beta)
  gamma<- exp(gamma)
  
  # Simulate SIR model using current parameter values
  model_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  
  # Compute sum of squared differences between simulated removals and observed removals
  OLS = sum((c(0,diff(model_data$R)) - observed_data$Removals)^2)
  
  #Alternatively using RMSE
  #RMSE <- sqrt(mean((c(0, diff(model_data$R)) - observed_data$Removals)^2))
  #return(RMSE)
  
  return(OLS)
}

# Perform optimization to estimate model parameters
initial_params <- c(log(0.0000001), log(0.5))  # Initial parameter values
#Optim
opt_result <- optim(par = initial_params, fn = objective_function, method='Nelder-Mead')
print(opt_result)
estimated_params <- opt_result$par
estimated_params<- exp(estimated_params)

#nlm
#opt_result <- nlm(p = initial_params, f = objective_function)
#print(opt_result)
#estimated_params <- opt_result$estimate
#estimated_params<- exp(estimated_params)


# Simulate SIR model using estimated parameters
predicted_data <- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, estimated_params[1], estimated_params[2], step_size)
predicted_data$ActualRemovals<- c(0,diff(predicted_data$R))

Initpredicted_data<- Deterministic_DT_SIR_model(N, S0, I0, minTime, maxTime, exp(initial_params[1]), exp(initial_params[2]), step_size)
Initpredicted_data$ActualRemovals<- c(0,diff(Initpredicted_data$R))

# Plot observed and simulated removals  
library(ggplot2)
ggplot() +
  geom_point(data = observed_data, aes(x = Steps, y = Removals, color = "Observed removals"), size = 1) +
  geom_line(data = predicted_data, aes(x = Steps, y = ActualRemovals, color = "Predicted using optimized values"), linewidth = 1) +
  geom_line(data = Initpredicted_data, aes(x = Steps, y = ActualRemovals, color = "Predicted using initial values"), linewidth = 1) +
  ylim(0, 50) +
  labs(x = "Day", y = "Removals", title = "") +
  scale_color_manual(values = c("Observed removals" = "blue", "Predicted using optimized values" = "red", "Predicted using initial values" = "black")) +
  theme_minimal()

estimated_params  #[1] \beta=0.09975901  \gamma=0.04880698


#data.frame(observed_data$Removals, predicted_data$ActualRemovals, Initpredicted_data$ActualRemovals)
