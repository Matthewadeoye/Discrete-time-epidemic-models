
#First trial & correct
Deterministic_DT_SIR_model1<- function(N,S0,I0,time,beta,gamma){
  S<- numeric() #create empty vector to store simulations for susceptibles
  I<- numeric() #create empty vector to store simulations for infectives
  R<- numeric() #create empty vector to store simulations for removed
  S[1]<- S0  #initial number of susceptibles
  I[1]<- I0   #initial number of infectives
  R[1]<- N-S0-I0      #initial number of Removed
  t<- time           #vector of time
  b<- beta           #infection rate
  g<- gamma          #recovery rate
  #loop and update compartments
  for(tym in 2:length(t)){
    S[tym]<- S[tym-1]-b*S[tym-1]*I[tym-1]
    I[tym]<- I[tym-1]+b*S[tym-1]*I[tym-1]-g*I[tym-1]
    R[tym]<- R[tym-1]+g*I[tym-1]
  }
  #plot results
  plot(t, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Discrete-time deterministic SIR Model")
  lines(t, I, col = "red")
  lines(t, R, col = "green")
  legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)
}
Deterministic_DT_SIR_model1(N=11000,S0=10000,I0=1000,time=(1:400),beta=0.0001,gamma=0.05,step_size=1)


#Second trial assuming infections occur at the point of a Poisson process
Deterministic_DT_SIR_model2<- function(N,S0,I0,time,beta,gamma,step_size){
S<- numeric() #create empty vector to store simulations for susceptibles
I<- numeric() #create empty vector to store simulations for infectives
R<- numeric() #create empty vector to store simulations for removed
b<- beta
g<- gamma
t <- time  # A vector of total time period
h <- step_size   # Time discretization step size 
total_steps<- length(time)  #Total number of simulations 

# The initial susceptibe, infected and recovered at the start
S[1] <- S0
I[1] <- I0
R[1] <- N-S0-I0

#Loop through discretized time and Update the compartments at each step
for (tym in 2:total_steps){
#Assuming infections occur at the point of a Poisson process
S[tym]<- S[tym-1] - S[tym-1]*(1-exp(-b*I[tym-1]*h))
I[tym]<- I[tym-1] + S[tym-1]*(1-exp(-b*I[tym-1]*h))-I[tym-1]*(1-exp(-g*h))
R[tym]<- R[tym-1] + I[tym-1]*(1-exp(-g*h))
}
#plot results
plot(t, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Discrete-time deterministic SIR Model")
lines(t, I, col = "red")
lines(t, R, col = "green")
legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)
}
Second_DT_SIR_model(N=11000,S0=10000,I0=1000,time=(1:400),beta=0.0001,gamma=0.05,step_size=1)




#Stochastic simulation
Stochastic_DT_SIR_model<- function(N,S0,I0,time,beta,gamma,step_size){
t <- time  # Total period in days
h <- step_size  # Time discretization step size
total_steps <- length(t)  # Total number of time steps

# Initialize compartments
S <- numeric()
I <- numeric()
R <- numeric()

# Initial values
b <- beta
g <- gamma
S[1] <- S0
I[1] <- I0
R[1] <- N-S0-I0

# Simulation loop
for (tym in 2:total_steps) {
  # Assuming infections occur at the point of a Poisson process
  p_I <- 1 - exp(-b * I[tym-1] * h)  # Probability of infection
  p_R <- 1 - exp(-g * h)  # Probability of recovery
  
  # sample infections and recoveries from a binomial trial
  delta_S <- rbinom(1, S[tym-1], p_I)
  delta_I <- rbinom(1, I[tym-1], p_R)
  
  # Update the compartments
  S[tym] <- S[tym-1] - delta_S
  I[tym] <- I[tym-1] + delta_S - delta_I
  R[tym] <- R[tym-1] + delta_I
}

plot(t, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Stochastic SIR Model")
lines(t, I, col = "red")
lines(t, R, col = "green")
legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)
}
Stochastic_DT_SIR_model(N=11000,S0=10000,I0=1000,time=(1:400),beta=0.0001,gamma=0.05,step_size=1)
  

