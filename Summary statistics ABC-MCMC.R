#ABC-MCMC using Summary Statistics constructed from linear regression 

#Bayesian Inference for Deterministic SIR

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
num_iterations <- 50000

# Saving parameter samples
posterior.samples <- matrix(0, nrow = num_iterations, ncol = 2)

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


#Constructing summary statistics with linear regression

N<- 10000
bdata<- numeric(N)
gdata<- numeric(N)
designMatrix<- matrix(0, nrow=N, ncol=365)
colnames(designMatrix) <- paste(rep("X=", 365),seq(1,365), sep="")

for(i in 1:N){
  params<- prior_dist(1)
  bdata[i]<- params[1]
  gdata[i]<- params[2]
  
  designMatrix[i,]<- model(params)
}

designMatrix<- as.data.frame(designMatrix)
Observed_data<- as.data.frame(t(Observed_data))
colnames(Observed_data) <- paste(rep("X=", 365),seq(1,365), sep="")

#Summary statistic for Beta
data1<- cbind(bdata,designMatrix)
regression1<- lm(bdata~ ., data=data1)
Observed_statisticBeta <- predict(regression1, newdata = Observed_data)


#Summary statistic for Gamma
data2<- cbind(gdata,designMatrix)
regression2<- lm(gdata~ ., data=data2)
Observed_statisticGamma <- predict(regression2, newdata = Observed_data)

Observed_statisticBeta
Observed_statisticGamma


# Define a function to calculate the distance metric (L2 norm)
distance <- function(x1, x2, y1, y2) {
  e<- sqrt(((x1 - y1)^2) + ((x2 - y2)^2))
  return(e)
}


# ABC-MCMC algorithm

epsilon <- 0.0007

# Define number of iterations
num_iterations <- 10000

chain <- matrix(0, nrow = num_iterations, ncol = 2)


# Initialize parameter values
chain[1,]<- c(grid[row.number, 1], grid[row.number, 2])

for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], chain[1,]))
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.001, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.001, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  
  mh.prob <- exp(prior_proposed - prior_current)
  
  # Accept (and simulate) or reject the proposals
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    
    # Generate synthetic data based on the proposed parameters
    synthetic_data <- model(proposed_params)
    
    synthetic_data<- as.data.frame(t(synthetic_data))
    colnames(synthetic_data) <- paste(rep("X=", 365),seq(1,365), sep="")
    
    synthetic_statisticBeta<- predict(regression1, newdata = synthetic_data)
    synthetic_statisticGamma<- predict(regression2, newdata = synthetic_data)
    
    
    # Calculate the distance metric between synthetic and observed summary statistics
    metric<- distance(Observed_statisticBeta, Observed_statisticGamma, synthetic_statisticBeta, synthetic_statisticGamma)
    
    
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

#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area
#par(mar=c(4,2,2,2))
hist(infection.rate, freq=F, xlab = expression(beta), main = "", col = "white", ylab="Density")
abline(v=0.0001, col="red", lwd=2,lty=1)
abline(h=1, col="blue", lwd=2,lty=1)
plot(chain[, 1], type = "l", ylab=expression(beta), main = "", xlab="Index", ylim=c(0,0.00035))
mtext("(a)", side = 3, line = -26, outer = TRUE)
plot(infection.rate,recovery.rate, xlab=expression(beta),ylab=expression(gamma))
hist(recovery.rate, freq=F, xlab = expression(gamma), main = "", col = "white",ylab="Density")
abline(v=0.05, col="red", lwd=2,lty=1)
abline(h=1, col="blue", lwd=2,lty=1)
plot(chain[, 2], type = "l", ylab=expression(gamma), main = "", xlab="Index",ylim=c(0.02,0.15))
mtext("(b)", side = 1, line = -1, outer = TRUE)
plot(recovery.rate,infection.rate, xlab=expression(gamma),ylab=expression(beta))
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


# Create the contour plot
x <- seq(0.00003, 0.0003, 0.00005)
y <- seq(0.02, 0.1, 0.005)
u <- as.matrix(expand.grid(x, y))

summary <- function(x, y) {
  synthetic_data <- model(c(x, y))
  synthetic_data <- as.data.frame(t(synthetic_data))
  colnames(synthetic_data) <- paste(rep("X=", 365), seq(1, 365), sep = "")
  
  synthetic_statisticBeta <- predict(regression1, newdata = synthetic_data)
  synthetic_statisticGamma <- predict(regression2, newdata = synthetic_data)
  
  return(list(synthetic_statisticBeta, synthetic_statisticGamma)) # Return the results as a list
}

z <- matrix(apply(u, 1, function(v) distance(summary(v[1], v[2])[[1]], Observed_statisticBeta, summary(v[1], v[2])[[2]], Observed_statisticGamma)), 
            nrow = length(x))

filled.contour(x,y,z, col=rainbow(39, 1, rev = FALSE), 
               xlab="Infection rate (Beta)",
               ylab="Recovery rate (Gamma)",
               main="Distance metric between observed and simulated data")



##################################################################################################################

#Bayesian Inference for Stochastic SIR

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
num_iterations <- 50000

# Saving parameter samples
posterior.samples <- matrix(0, nrow = num_iterations, ncol = 2)

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


#Constructing summary statistics with linear regression

N<- 10000
bdata<- numeric(N)
gdata<- numeric(N)
designMatrix<- matrix(0, nrow=N, ncol=365)
colnames(designMatrix) <- paste(rep("X=", 365),seq(1,365), sep="")

for(i in 1:N){
  params<- prior_dist(1)
  bdata[i]<- params[1]
  gdata[i]<- params[2]
  
  designMatrix[i,]<- model(params)
}

designMatrix<- as.data.frame(designMatrix)
Observed_data<- as.data.frame(t(Observed_data))
colnames(Observed_data) <- paste(rep("X=", 365),seq(1,365), sep="")

#Summary statistic for Beta
data1<- cbind(bdata,designMatrix)
regression1<- lm(bdata~ ., data=data1)
Observed_statisticBeta <- predict(regression1, newdata = Observed_data)


#Summary statistic for Gamma
data2<- cbind(gdata,designMatrix)
regression2<- lm(gdata~ ., data=data2)
Observed_statisticGamma <- predict(regression2, newdata = Observed_data)

Observed_statisticBeta
Observed_statisticGamma


# Define a function to calculate the distance metric (L2 norm)
distance <- function(x1, x2, y1, y2) {
  e<- sqrt(((x1 - y1)^2) + ((x2 - y2)^2))
  return(e)
}


# ABC-MCMC algorithm

epsilon <- 0.0007

# Define number of iterations
num_iterations <- 10000

chain <- matrix(0, nrow = num_iterations, ncol = 2)


# Initialize parameter values
chain[1,]<- c(grid[row.number, 1], grid[row.number, 2])

for (i in 2:num_iterations) {
  # Generate parameter proposals from the proposal kernels
  proposed_params <- abs(rnorm(2, chain[i - 1, ], chain[1,]))
  
  prior_proposed<- dunif(proposed_params[1], min = 0, max = 0.001, log = TRUE) + 
    dunif(proposed_params[2], min = 0, max = 0.1, log = TRUE)
  prior_current<- dunif(chain[i - 1, 1], min = 0, max = 0.001, log = TRUE) + 
    dunif(chain[i - 1, 2], min = 0, max = 0.1, log = TRUE) 
  
  
  mh.prob <- exp(prior_proposed - prior_current)
  
  # Accept (and simulate) or reject the proposals
  if (!is.na(mh.prob) && runif(1) < mh.prob) {
    
    # Generate synthetic data based on the proposed parameters
    synthetic_data <- model(proposed_params)
    
    synthetic_data<- as.data.frame(t(synthetic_data))
    colnames(synthetic_data) <- paste(rep("X=", 365),seq(1,365), sep="")
    
    synthetic_statisticBeta<- predict(regression1, newdata = synthetic_data)
    synthetic_statisticGamma<- predict(regression2, newdata = synthetic_data)
    
    
    # Calculate the distance metric between synthetic and observed summary statistics
    metric<- distance(Observed_statisticBeta, Observed_statisticGamma, synthetic_statisticBeta, synthetic_statisticGamma)
    
    
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

#Plots
par(mfrow=c(2,3)) # Create a 2x3 plotting area
#par(mar=c(4,2,2,2))
hist(infection.rate, freq=F, xlab = expression(beta), main = "", col = "white", ylab="Density")
abline(v=0.0001, col="red", lwd=2,lty=1)
abline(h=1, col="blue", lwd=2,lty=1)
plot(chain[, 1], type = "l", ylab=expression(beta), main = "", xlab="Index", ylim=c(0,0.00035))
mtext("(a)", side = 3, line = -26, outer = TRUE)
plot(infection.rate,recovery.rate, xlab=expression(beta),ylab=expression(gamma))
hist(recovery.rate, freq=F, xlab = expression(gamma), main = "", col = "white",ylab="Density")
abline(v=0.05, col="red", lwd=2,lty=1)
abline(h=1, col="blue", lwd=2,lty=1)
plot(chain[, 2], type = "l", ylab=expression(gamma), main = "", xlab="Index",ylim=c(0.02,0.15))
mtext("(b)", side = 1, line = -1, outer = TRUE)
plot(recovery.rate,infection.rate, xlab=expression(gamma),ylab=expression(beta))
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



# Create the contour plot
x <- seq(0.00003, 0.0003, 0.00005)
y <- seq(0.02, 0.1, 0.005)
u <- as.matrix(expand.grid(x, y))

summary <- function(x, y) {
  synthetic_data <- model(c(x, y))
  synthetic_data <- as.data.frame(t(synthetic_data))
  colnames(synthetic_data) <- paste(rep("X=", 365), seq(1, 365), sep = "")
  
  synthetic_statisticBeta <- predict(regression1, newdata = synthetic_data)
  synthetic_statisticGamma <- predict(regression2, newdata = synthetic_data)
  
  return(list(synthetic_statisticBeta, synthetic_statisticGamma)) # Return the results as a list
}

z <- matrix(apply(u, 1, function(v) distance(summary(v[1], v[2])[[1]], Observed_statisticBeta, summary(v[1], v[2])[[2]], Observed_statisticGamma)), 
            nrow = length(x))

filled.contour(x,y,z, col=rainbow(39, 1, rev = FALSE), 
               xlab="Infection rate (Beta)",
               ylab="Recovery rate (Gamma)",
               main="Distance metric between observed and simulated data")
