
#First version
Deterministic_DT_SIR_model1<- function(N,S0,I0,time,beta,gamma,step_size){
  S<- numeric() #create empty vector to store simulations for susceptibles
  I<- numeric() #create empty vector to store simulations for infectives
  R<- numeric() #create empty vector to store simulations for removed
  S[1]<- S0  #initial number of susceptibles
  I[1]<- I0   #initial number of infectives
  R[1]<- N-S0-I0      #initial number of Removed
  t<- time           #vector of time
  h<- step_size      #step size
  steps<- seq(min(t),max(t),by=h)   #time discretization
  b<- beta           #infection rate
  g<- gamma          #recovery rate
  #loop and update compartments
  for(tym in 2:length(steps)){
    S[tym]<- S[tym-1]-(b*S[tym-1]*I[tym-1])*h
    I[tym]<- I[tym-1]+(b*S[tym-1]*I[tym-1]-g*I[tym-1])*h
    R[tym]<- R[tym-1]+g*I[tym-1]*h
  }
  #results
  simulations <- data.frame(S = S, I = I, R = R)
  
  graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Deterministic discrete-time SIR Model1")
lines(steps, I, col = "red")
lines(steps, R, col = "green")
legend("topright", legend = c("Susceptible", "Infected", "Removed"), col = c("blue", "red", "green"), lty = 1)
  
return (list(simulations, graph))
}
Deterministic_DT_SIR_model1(N=11000,S0=10000,I0=1000,time=(1:1000),beta=0.0001,gamma=0.05,step_size=0.03)


#Second version assuming infections occur at the point of a Poisson process
Deterministic_DT_SIR_model2<- function(N,S0,I0,time,beta,gamma,step_size){
S<- numeric() #create empty vector to store simulations for susceptibles
I<- numeric() #create empty vector to store simulations for infectives
R<- numeric() #create empty vector to store simulations for removed
b<- beta
g<- gamma
t<- time           #vector of time
h<- step_size      #step size
steps<- seq(min(t),max(t),by=h)   #time discretization

# The initial susceptibe, infected and recovered at the start
S[1] <- S0
I[1] <- I0
R[1] <- N-S0-I0

#Loop through discretized time and Update the compartments at each step
for (tym in 2:length(steps)){
#Assuming infections occur at the point of a Poisson process
S[tym]<- S[tym-1] - S[tym-1]*(1-exp(-b*I[tym-1]*h))
I[tym]<- I[tym-1] + S[tym-1]*(1-exp(-b*I[tym-1]*h))-I[tym-1]*(1-exp(-g*h))
R[tym]<- R[tym-1] + I[tym-1]*(1-exp(-g*h))
}

#results
simulations <- data.frame(S = S, I = I, R = R)

graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Deterministic discrete-time SIR Model2")
lines(steps, I, col = "red")
lines(steps, R, col = "green")
legend("topright", legend = c("Susceptible", "Infected", "Removed"), col = c("blue", "red", "green"), lty = 1)
return (list(simulations, graph))
}
Deterministic_DT_SIR_model2(N=11000,S0=10000,I0=1000,time=(1:1000),beta=0.0001,gamma=0.05,step_size=0.03)



#Stochastic simulation
Stochastic_DT_SIR_model<- function(N,S0,I0,time,beta,gamma,step_size){
t<- time           #vector of time
h<- step_size      #step size
steps<- seq(min(t),max(t),by=h)   #time discretization

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
for (tym in 2:length(steps)) {
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
simulations <- data.frame(S = S, I = I, R = R)


graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Stochastic discrete-time SIR Model")
lines(steps, I, col = "red")
lines(steps, R, col = "green")
legend("topright", legend = c("Susceptible", "Infected", "Removed"), col = c("blue", "red", "green"), lty = 1)
  
 return (list(simulations, graph))
}
Stochastic_DT_SIR_model(N=11000,S0=10000,I0=1000,time=(1:1000),beta=0.0001,gamma=0.05,step_size=0.03)
  

