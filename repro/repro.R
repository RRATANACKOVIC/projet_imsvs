


###################################################################
###################################################################
## EXERCICE 1
###################################################################

bac.model <- function(t, pop, param) {
  
  # You need three bacteria populations: Recipient (no plasmid), Donor (with plasmid),
  # Conjugate (acquired plasmid from Donor)
  X <- pop[1]  # Recipient
  P <- pop[2]  # Donor
  S <- pop[3]  # Conjugate
  
  # Parameters required for: max growth rate, carrying capacity, conjugation rate,
  # death rate
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

dt     = 1
Tmax   = 15

X0     = 5
P0     = 10
S0     = 0

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
result[,2:4] = log10(result[,2:4])

# Plot each line
plot(Time,result[,"X"],type="l",col="red",xlab="Time",ylab="Bacteria",ylim=c(0,log10(Bmax)))
lines(Time,result[,"P"],col="blue")
lines(Time,result[,"S"],col="green")
legend("bottomright",legend=c("X","P","S"),col=c("red","blue","green"),lty=1)


# BONUS QUESTION ##
Tmax=500  # set a long simulation time
Time=seq(from=0,to=Tmax,by=dt)
result <- lsoda(Init.cond, Time, bac.model, param)
colnames(result) <- c("Time", "X", "P", "S")
result[(min(which(result[,"X"]<1))),1] #find the first time-step where R < 1
# note: the advantage of returning a value of result[,1] here instead of just
# min(which(...))is that this solution works if you change the dt and Tmax parameters!



###################################################################
###################################################################
## EXERCICE 2
###################################################################

bac.model <- function(t, pop, param) {
  
  # You need three bacteria populations: Recipient (no plasmid), Donor (with plasmid),
  # Conjugated (acquired plasmid from Donor), Antibiotic concentration
  R <- pop[1]   # Recipient
  D <- pop[2]   # Donor
  C <- pop[3]   # Conjugate
  A <- pop[4]   # Antibiotic
  
  # Parameters required for: max growth rate, carrying capacity, conjugation rate,
  # death rate, max growth rate for plasmid-carriers, max effect of antibiotic,
  # EC50 for antibiotic, antibiotic removal rate
  mu=param[1]       # max growth rate
  Bmax=param[2]     # carrying capacity
  beta=param[3]     # conjugation rate
  gamma=param[4]    # death rate
  mu_p=param[5]     # max growth rate for plasmid-carriers
  emax=param[6]     # max effect of antibiotic
  EC50=param[7]     # EC50 for antibiotic
  H=param[8]        # Hill coefficient
  gamma_a=param[9]  # antibiotic removal rate
  
  # Four differential equations needed, one for each bacteria population, one for
  # antibiotic
  B=R+D+C   # total bacteria
  
  e_t = emax*A^H/(A^H+EC50^H) # antibiotic effect at time t
  
  dR <- mu*(1-B/Bmax)*R - beta*(D+C)/B*R - gamma*R - e_t*R
  dD <- mu_p*(1-B/Bmax)*D - gamma*D
  dC <- mu_p*(1-B/Bmax)*C + beta*(D+C)/B*R - gamma*C
  
  dA <- -gamma_a*A
  
  res<-c(dR, dD, dC, dA)
  
  list(res)
}


# Load the deSolve package
library(deSolve)

# Define parameters
mu=1.5
Bmax=1e9
beta=0.1
gamma=0.5

mu_p=1
emax=3
EC50=2
H=3
gamma_a=0.5

dt=1
Tmax=30

R0=1e4
D0=1e4
C0=0
A0=10

# Define model time-steps
Time=seq(from=0,to=Tmax,by=dt)

# Define vectors for initial conditions and parameters
Init.cond=c(R0,D0,C0,A0) 
param=c(mu,Bmax,beta,gamma,mu_p,emax,EC50,H,gamma_a)


# Execute
result <- lsoda(Init.cond, Time, bac.model, param)

# Name columns for convenience
colnames(result) <- c("Time", "R", "D", "C", "A")

# Quick look
head(result)

# Log-transform
result[,2:4] = log10(result[,2:4])

# Plot each line
par(mfrow=c(2,1)) # splits up the plotting area in two lines, one column
plot(Time,result[,"R"],type="l",col="red",xlab="Time",ylab="Bacteria",ylim=c(0,log10(Bmax)))
lines(Time,result[,"D"],col="blue")
lines(Time,result[,"C"],col="green")
legend("bottomright",legend=c("R","D","C"),col=c("red","blue","green"),lty=1,cex=0.5)

plot(Time,result[,"A"],type="l",col="black",xlab="Time",ylab="Antibiotic concentration")
par(mfrow=c(1,1)) # restores the plotting area to one line, one column


# BONUS QUESTION ##
# set a starting value to explore A0
A0 = 10

repeat{
  A0 = A0+0.1 # increase emax progressively
  Init.cond=c(R0,D0,C0,A0) 
  result <- lsoda(Init.cond, Time, bac.model, param)
  if(any(result[,2] < 1)) break # stop as soon as any value for R (column 2) is less than 1
}

A0

colnames(result) <- c("Time", "R", "D", "C", "A")
result[,2:4] = log10(result[,2:4])

par(mfrow=c(2,1)) # splits up the plotting area in two lines, one column
plot(Time,result[,"R"],type="l",col="red",xlab="Time",ylab="Bacteria",ylim=c(0,log10(Bmax)))
lines(Time,result[,"D"],col="blue")
lines(Time,result[,"C"],col="green")
legend("bottomright",legend=c("R","D","C"),col=c("red","blue","green"),lty=1,cex=0.5)
plot(Time,result[,"A"],type="l",col="black",xlab="Time",ylab="Antibiotic concentration")
par(mfrow=c(1,1)) # restores the plotting area to one line, one column


# BONUS QUESTION ##
concentrations = seq(0.1,32,0.1)  # concentrations to test

plot(concentrations, emax*(concentrations^H/(concentrations^H+10^H)), type="l",
     ylab="Antibiotic action rate", xlab="Antibiotic concentration")
