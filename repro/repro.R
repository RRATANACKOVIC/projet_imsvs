


###################################################################
###################################################################
## EXERCICE 1
###################################################################

bac.model <- function(t, pop, param) {
  
  X <- pop[1]  # biomass concentration
  P <- pop[2]  # ethanol concentration
  S <- pop[3]  # substrate concentration
  
  mu_max  = param[1]  # max growth rate
  q_max   = param[2]  # carrying capacity
  Yxs     = param[3]  # conjugation rate
  Yps     = param[4]  # death rate
  Ksx     = param[5]
  Ksp     = param[6]
  
  # Three differential equations needed, one for each bacteria population
  #B=R+D+C   # total bacteria
  
  dX <- mu_max*(S/(S+Ksx))*X
  dP <- q_max*(S/(S+Ksp))*X
  dS <- -(1/Yxs)*dX - (1/Yps)*dP
  
  res<-c(dX, dP, dS)
  
  list(res)
}


# Load the deSolve package
library(deSolve)

# Define parameters
mu_max = 0.424
q_max  = 2.042
Yxs    = 0.125
Yps    = 0.53
Ksx    = 150
Ksp    = 150

dt     = 0.5
Tmax   = 9

X0     = 1.5
P0     = 0
S0     = 10

# Define model time-steps
Time=seq(from=0,to=Tmax,by=dt)

# Define vectors for initial conditions and parameters
Init.cond=c(X0,P0,S0) 
param=c(mu_max,q_max,Yxs,Yps,Ksx,Ksp)

# Execute
result <- lsoda(Init.cond, Time, bac.model, param)

# Name columns for convenience
colnames(result) <- c("Time", "X", "P", "S")

# Quick look
head(result)

# Log-transform
#result[,2:4] = log10(result[,2:4])

# Plot each line
plot(Time,result[,"X"],type="l",col="red",xlab="Time (d)",ylab="concentration [g/L]", ylim = c(0,1.1*max(result)))
lines(Time,result[,"P"],col="blue")
lines(Time,result[,"S"],col="green")
legend("topright",legend=c("X","P","S"),col=c("red","blue","green"),lty=1)

