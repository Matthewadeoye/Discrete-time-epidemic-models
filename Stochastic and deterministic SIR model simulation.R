
#First version
Deterministic_DT_SIR_model1<- function(N,S0,I0,minTime,maxTime,beta,gamma,step_size,plot=TRUE){
  minTime<- minTime           #start time
  maxTime<- maxTime         #end time
  step_size<- step_size      #step size
  steps<- seq(minTime,maxTime, by=step_size)   #time discretization
  beta<- beta           #infection rate
  gamma<- gamma          #recovery rate
  S<- numeric(length(steps)) #create empty vector to store simulations for susceptibles
  I<- numeric(length(steps)) #create empty vector to store simulations for infectives
  R<- numeric(length(steps)) #create empty vector to store simulations for removed
  S[1]<- S0  #initial number of susceptibles
  I[1]<- I0   #initial number of infectives
  R[1]<- N-S0-I0      #initial number of Removed
  #loop and update compartments
  for(tym in 2:length(steps)){
    S[tym]<- S[tym-1]-(beta*S[tym-1]*I[tym-1])*step_size
    I[tym]<- I[tym-1]+(beta*S[tym-1]*I[tym-1]-gamma*I[tym-1])*step_size
    R[tym]<- R[tym-1]+gamma*I[tym-1]*step_size
  }
  #results
  simulations <- data.frame(Steps=steps, S = S, I = I, R = R)
  if(plot){
  graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Deterministic discrete-time SIR Model1",ylim=c(0,N))
  lines(steps, I, col = "red")
  lines(steps, R, col = "green")
  legend("topright", legend = c("Susceptible", "Infected", "Removed"), col = c("blue", "red", "green"), lty = 1)
  
  return (list(simulations, graph))
}else{
  return(list(simulations=simulations))
}
}
Deterministic_DT_SIR_model1(N=1000,S0=900,I0=100,minTime=1,maxTime=100,beta=0.001,gamma=0.05, step_size=0.03,plot=T)



#Second version assuming infections occur at the point of a Poisson process
Deterministic_DT_SIR_model2<- function(N,S0,I0,minTime,maxTime,beta,gamma,step_size,plot=TRUE){
  minTime<- minTime           #start time
  maxTime<- maxTime           #end time
  step_size<- step_size      #step size
  steps<- seq(minTime,maxTime, by=step_size)   #time discretization
  beta<- beta            #infection rate
  gamma<- gamma          #recovery rate
  S<- numeric(length(steps)) #create empty vector to store simulations for susceptibles
  I<- numeric(length(steps)) #create empty vector to store simulations for infectives
  R<- numeric(length(steps)) #create empty vector to store simulations for removed
  S[1]<- S0  #initial number of susceptibles
  I[1]<- I0   #initial number of infectives
  R[1]<- N-S0-I0      #initial number of Removed
  #Loop through discretized time and Update the compartments at each step
  for (tym in 2:length(steps)){
    #Assuming infections occur at the point of a Poisson process
    S[tym]<- S[tym-1] - S[tym-1]*(1-exp(-beta*I[tym-1]*step_size))
    I[tym]<- I[tym-1] + S[tym-1]*(1-exp(-beta*I[tym-1]*step_size))-I[tym-1]*(1-exp(-gamma*step_size))
    R[tym]<- R[tym-1] + I[tym-1]*(1-exp(-gamma*step_size))
  }
  #results
  simulations <- data.frame(Steps=steps,S = S, I = I, R = R)
  if(plot){
  graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Deterministic discrete-time SIR Model2",ylim=c(0,N))
  lines(steps, I, col = "red")
  lines(steps, R, col = "green")
  legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)
  return (list(simulations, graph))
}else{
  return(list(simulations=simulations))
}
}
Deterministic_DT_SIR_model2(N=1000,S0=900,I0=100,minTime=1,maxTime=1000,beta=0.0001,gamma=0.05,step_size=0.03,plot=T)



#Stochastic simulation
Stochastic_DT_SIR_model<- function(N,S0,I0,minTime,maxTime,beta,gamma,step_size,plot=TRUE){
  minTime<- minTime           #start time
  maxTime<- maxTime         #end time
  step_size<- step_size      #step size
  steps<- seq(minTime,maxTime, by=step_size)   #time discretization
  beta<- beta           #infection rate
  gamma<- gamma          #recovery rate
  S<- numeric(length(steps)) #create empty vector to store simulations for susceptibles
  I<- numeric(length(steps)) #create empty vector to store simulations for infectives
  R<- numeric(length(steps)) #create empty vector to store simulations for removed
  S[1]<- S0  #initial number of susceptibles
  I[1]<- I0   #initial number of infectives
  R[1]<- N-S0-I0      #initial number of Removed
  #Loop through discretized time and Update the compartments at each step
  # Simulation loop
  for (tym in 2:length(steps)) {
    # Assuming infections occur at the point of a Poisson process
    p_I <- 1 - exp(-beta*I[tym-1]*step_size)  # Probability of infection
    p_R <- 1 - exp(-gamma*step_size)  # Probability of recovery
    
    # sample infections and recoveries from a binomial trial
    delta_S <- rbinom(1, S[tym-1], p_I)
    delta_I <- rbinom(1, I[tym-1], p_R)
    
    # Update the compartments
    S[tym] <- S[tym-1] - delta_S
    I[tym] <- I[tym-1] + delta_S - delta_I
    R[tym] <- R[tym-1] + delta_I
  }
  simulations <- data.frame(Steps=steps,S = S, I = I, R = R)
  if(plot){
  graph<- plot(steps, S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Stochastic discrete-time SIR Model",ylim=c(0,N))
  lines(steps, I, col = "red")
  lines(steps, R, col = "green")
  legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)
  return (list(simulations, graph))
}else{
  return(list(simulations=simulations))
}
}
Stochastic_DT_SIR_model(N=1000,S0=900,I0=100,minTime=1,maxTime=365,beta=0.0001,gamma=0.05,step_size=0.03,plot=T)
