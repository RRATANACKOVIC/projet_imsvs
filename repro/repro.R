

bac.model <- function(t, pop, param) {
  
  X  <- pop[1]  # biomass concentration
  P  <- pop[2]  # ethanol concentration
  S  <- pop[3]  # substrate concentration
  
  FA <- pop[4]  # higher alcohol
  E  <- pop[5]  # esters
  A  <- pop[6]  # aldehydes
  
  mu_max  = param[1]  # max growth rate
  q_max   = param[2]  
  Ksx     = param[3]
  Ksp     = param[4]
  Yxs     = param[5]  
  Yps     = param[6]  

  
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
library(rODE)
library(rjson)

#reading json input file

list_param <- fromJSON(file = "parameters.json")

# Define parameters


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
param = list_param$constants[[1]]$values
# Execute
result <- rk(Init.cond, Time, bac.model, param, method = "rk45dp7")

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
legend("topright",legend=c("X","P","S"),col=c("red","blue","green"),lty=1)


plot(Time,result[,"FA"],type="l",col="purple",xlab="Time (d)",ylab="concentration [g/L]")
lines(Time,result[,"E"],col="cyan")
lines(Time,result[,"A"],col="orange")
legend("bottomright",legend=c("FA","E","A"),col=c("purple","cyan","orange"),lty=1)


