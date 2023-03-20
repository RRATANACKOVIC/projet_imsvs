

bac.model <- function(t, pop, param) {
  
  X  <- pop[1]  # biomass concentration
  P  <- pop[2]  # ethanol concentration
  S  <- pop[3]  # substrate concentration
  
  FA <- pop[4]  # higher alcohol
  E  <- pop[5]  # esters
  A  <- pop[6]  # aldehydes
  
  mu_max  = param[1]  # max growth rate
  q_max   = param[2]  # carrying capacity
  Yxs     = param[3]  # conjugation rate
  Yps     = param[4]  # death rate
  Ksx     = param[5]
  Ksp     = param[6]
  
  Yfa     = param[7]
  Ye      = param[8]
  Ya      = param[9]
  Ka      = param[10]
  

  
  dX  <- mu_max * (S/(S+Ksx)) * X
  dP  <- q_max * (S/(S+Ksp)) * X
  dS  <- -(1/Yxs) * dX - (1/Yps)*dP
  
  dFA <- Yfa * mu_max*(S/(S+Ksx)) * X 
  dE  <- Ye * mu_max*(S/(S+Ksx)) * X 
  dA  <- Ya * mu_max*(S/(S+Ksx)) * X - Ka * A * X
  
  res<-c(dX, dP, dS, dFA, dE, dA)
  
  list(res)
}

require('rODE')

# Load the deSolve package
library(deSolve)
library('rODE')
# Define parameters
mu_max = 0.424
q_max  = 2.042
Yxs    = 0.125
Yps    = 0.53
Ksx    = 150
Ksp    = 150

Yfa    = 12.55
Ye     = 12.23
Ya     = 7.63
Ka     = 0.023

dt     = 0.1
Tmax   = 9

X0     = 1.5
P0     = 0
S0     = 10

FA0    = 0
E0     = 0
A0     = 0

# Define model time-steps
Time=seq(from=0,to=Tmax,by=dt)

# Define vectors for initial conditions and parameters
Init.cond=c(X0,P0,S0, FA0, E0, A0) 
param=c(mu_max,q_max,Yxs,Yps,Ksx,Ksp, Yfa, Ye, Ya, Ka)

# Execute
result <- ode(Init.cond, Time, bac.model, param)

# Name columns for convenience
colnames(result) <- c("Time", "X", "P", "S", "FA", "E", "A")

# Quick look
head(result)

# Log-transform
#result[,2:4] = log10(result[,2:4])

# Plot each line
plot(Time,result[,"X"],type="l",col="red",xlab="Time (d)",ylab="concentration [g/L]", ylim = c(0,1.1*max(result)))
lines(Time,result[,"P"],col="blue")
lines(Time,result[,"S"],col="green")

lines(Time,result[,"FA"],col="purple")
lines(Time,result[,"E"],col="cyan")
lines(Time,result[,"A"],col="orange")

legend("topright",legend=c("X","P","S"),col=c("red","blue","green"),lty=1)

