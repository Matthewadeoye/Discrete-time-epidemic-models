# Stochastic simulation
Stochastic_DT_SIR_model <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  minTime <- minTime           # start time
  maxTime <- maxTime           # end time
  step_size <- step_size       # step size
  steps <- seq(minTime, maxTime, by = step_size)   # time discretization
  beta <- beta                 # infection rate
  gamma <- gamma               # recovery rate
  S <- numeric(length(steps))  # create empty vector to store simulations for susceptibles
  I <- numeric(length(steps))  # create empty vector to store simulations for infectives
  R <- numeric(length(steps))  # create empty vector to store simulations for removed
  S[1] <- S0                   # initial number of susceptibles
  I[1] <- I0                   # initial number of infectives
  R[1] <- N - S0 - I0          # initial number of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, S[tym - 1], p_I)
    delta_I <- rbinom(1, I[tym - 1], p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S
    I[tym] <- I[tym - 1] + delta_S - delta_I
    R[tym] <- R[tym - 1] + delta_I
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


Deterministic_DT_SIR_model2 <- function(N, S0, I0, minTime, maxTime, beta, gamma, step_size) {
  minTime <- minTime           # start time
  maxTime <- maxTime           # end time
  step_size <- step_size       # step size
  steps <- seq(minTime, maxTime, by = step_size)   # time discretization
  beta <- beta                 # infection rate
  gamma <- gamma               # recovery rate
  S <- numeric(length(steps))  # create empty vector to store simulations for susceptibles
  I <- numeric(length(steps))  # create empty vector to store simulations for infectives
  R <- numeric(length(steps))  # create empty vector to store simulations for removed
  S[1] <- S0                   # initial number of susceptibles
  I[1] <- I0                   # initial number of infectives
  R[1] <- N - S0 - I0          # initial number of Removed
  
  # Loop through discretized time and Update the compartments at each step
  for (tym in 2:length(steps)){
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 1000
S0 <- 990
I0 <- 10
minTime <- 0
maxTime <- 50
beta <- 0.0005
gamma <- 0.1
step_size <- 0.0001

# Plot the results
library(ggplot2)
# Empty plot frame 
Stoch <- ggplot() +
  xlim(minTime, maxTime) +
  xlab("Time") +
  ylab("Population") +
  ggtitle(paste("Deterministic and stochastic discrete-time SIR simulations (beta =",beta, ", gamma =",gamma,")"))

# Overlaying plots
plots <- 10
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  Stoch <- Stoch +
    geom_line(data = stochastic, aes(x = Steps, y = S), color = "#add8e6") +
    geom_line(data = stochastic, aes(x = Steps, y = I), color = "#ffcccb") +
    geom_line(data = stochastic, aes(x = Steps, y = R), color = "#90ee90")
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
Bind_plots <- Stoch +
  geom_line(data = deterministic, aes(x = Steps, y = S, color = "Susceptibles")) +
  geom_line(data = deterministic, aes(x = Steps, y = I, color = "Infectives")) +
  geom_line(data = deterministic, aes(x = Steps, y = R, color = "Removed")) +
  scale_color_manual(values = c("Susceptibles" = "blue", "Infectives" = "red", "Removed" = "#013220")) +
  labs(color="Compartments")
Bind_plots
