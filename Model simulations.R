#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area

#A

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 500
S0 <- 499
I0 <- 1
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1),cex.main=1.0, xlab = "Day", ylab = "Proportion", main = "N = 500, I(0) = 1")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)


#B

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 1000
S0 <- 999
I0 <- 1
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1), xlab = "Day",cex.main=1.0, ylab = "Proportion", main = "N = 1000, I(0) = 1")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)



#C

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 50000
S0 <- 49999
I0 <- 1
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1), xlab = "Day", cex.main=1.0,ylab = "Proportion", main = "N = 50000, I(0) = 1")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)


#D

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 500
S0 <- 450
I0 <- 50
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1), xlab = "Day",cex.main=1.0, ylab = "Proportion", main = "N = 500, I(0) = 50")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)

#E

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 1000
S0 <- 900
I0 <- 100
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1), xlab = "Day",cex.main=1.0, ylab = "Proportion", main = "N = 1000, I(0) = 100")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)

#F

#Frequency-dependent formulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta * I[tym - 1] * step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma * step_size)              # Probability of recovery
    
    # Sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, round(S[tym - 1] * N), p_I)
    delta_I <- rbinom(1, round(I[tym - 1] * N), p_R)
    
    # Update the compartments
    S[tym] <- S[tym - 1] - delta_S / N
    I[tym] <- I[tym - 1] + delta_S / N - delta_I / N
    R[tym] <- R[tym - 1] + delta_I / N
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}


#Frequency-dependent formulation
# Deterministic simulation
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
  S[1] <- S0 / N               # initial proportion of susceptibles
  I[1] <- I0 / N               # initial proportion of infectives
  R[1] <- (N - S0 - I0) / N    # initial proportion of Removed
  
  # Loop through discretized time and update the compartments at each step
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    S[tym] <- S[tym - 1] - S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size))
    I[tym] <- I[tym - 1] + S[tym - 1] * (1 - exp(-beta * I[tym - 1] * step_size)) - I[tym - 1] * (1 - exp(-gamma * step_size))
    R[tym] <- R[tym - 1] + I[tym - 1] * (1 - exp(-gamma * step_size))
  }
  
  return(data.frame(Steps = steps, S = S, I = I, R = R))
}

# Initial values and parameters
N <- 50000
S0 <- 45000
I0 <- 5000
minTime <- 1
maxTime <- 365
beta <- 0.1
gamma <- 0.05
step_size <- 0.03

#par(mfrow=c(1,1)) # Create a 2x3 plotting area
# Empty plot frame
plot(0, type = "n", xlim = c(minTime, maxTime), ylim = c(0, 1), xlab = "Day",cex.main=1.0, ylab = "Proportion", main = "N = 50000, I(0) = 5000")
colors <- c("#add8e6", "#ffcccb", "#90ee90")

# Overlaying plots
plots <- 50
for (i in 1:plots) {
  stochastic <- Stochastic_DT_SIR_model(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
  lines(stochastic$Steps, stochastic$S, col = colors[1])
  lines(stochastic$Steps, stochastic$I, col = colors[2])
  lines(stochastic$Steps, stochastic$R, col = colors[3])
}

deterministic <- Deterministic_DT_SIR_model2(N, S0, I0, minTime, maxTime, beta, gamma, step_size)
lines(deterministic$Steps, deterministic$S, col = "blue", lwd = 2.0)
lines(deterministic$Steps, deterministic$I, col = "red", lwd = 2.0)
lines(deterministic$Steps, deterministic$R, col = "#013220", lwd = 2.0)

#legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
add_legend("bottomright", legend=c("Susceptible", "Infected", "Recovered"), lty=1, 
           col=c("blue", "red", "green"),
           horiz=F, bty='n', cex=0.8)


