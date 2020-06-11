## Code to compare random walk with different process model structures
##
## Code skeleton and random walk with normal process from
## https://github.com/afredston/EFI_tick_project/blob/master/RandomWalkNymph.R
##
## Random walk with log-normal process from
## https://github.com/afredston/EFI_tick_project/blob/master/siegmundTickProjectNotes.Rmd

library(ecoforecastR)
library(tidyverse)

dat <- read.csv("project/groupProjectData/nymphDataWithMet.csv") # nymph data
sites <- unique(dat$Site)[1] # sites vector
tau.add <- ci <- list() # storage
tau.add2 <- ci2 <- list() # storage

# pick site = 1

i = 1

# JAGS Code ---------------------------------------------------------------

# * Random Walk with Normal Process Model -----------------------------------

RandomWalk <- "
    model{
    #### Data Model
    for(t in 1:n){
      y[t] ~ dpois(x[t])
    }
    #### Process Model
    for(t in 2:n){
      x[t] ~ dnorm(x[t-1],tau_add) T(0,)
    }
    #### Priors
    x[1] ~ dpois(x_ic)
    tau_add ~ dgamma(a_add,r_add)
  }"

# * Random Walk with Log Normal Process Model -----------------------------------

RandomWalkLogNormal = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }
  
  #### Process Model
  for(t in 2:n){
    mu[t-1] = log(max(0.0001, x[t-1]))
    x[t] ~ dlnorm(mu[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dlnorm(x_ic,tau_add[1])
  tau_add ~ dgamma(a_add,r_add)

}"

# Model fitting -----------------------------------------------------------

 # subset data
  data.train <- dat %>%
    filter(trainValidate == "Train") %>%
    filter(Site == sites[i])
  
  data.val <- dat %>% 
    filter(trainValidate == "Validate") %>%
    filter(Site == sites[i])
  
  na.add <- nrow(data.val) # number of validation data points (number of NAs to add)
  
  y <- data.train$totalNymph
  y <- c(data.train$totalNymph, rep(NA, na.add)) # add NAs
  n <- length(y)

# * Random Walk with Normal Process Model -----------------------------------

  data <- list(
    y = y,
    n = n,
    x_ic = round(runif(1, 1, 15)),
    a_add = 0.001,
    r_add = 0.001
  )
  
  j.model <- jags.model(file = textConnection(RandomWalk),
                        data = data,
                        # inits = init,
                        n.chains = 3)
  
  jags.out <- coda.samples(
    model = j.model,
    variable.names = c("x", "tau_add"),
    n.iter = 15000
  )
  
  ## check burnin
  GBR <- gelman.plot(jags.out, ask = FALSE)
  burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  cat("Burnin after", burnin, "iterations\n")
  out <- window(jags.out, start = burnin)
  out.mat <- as.matrix(out)
  out.predict <- out.mat[, -1]
  
  # save process error
  tau.add[[i]] <- out.tau_add <- out.mat[, "tau_add"]
  
  ## state CIs
  ci[[i]] <- apply(out.predict, 2, quantile, c(0.025, 0.5, 0.975))

# * Random Walk with Log-Normal Process Model -----------------------------------

  data <- list(
    y = y,
    n = n,
    x_ic = log(round(runif(1, 1, 15))),
    a_add = 0.001,
    r_add = 0.001
  )

  j.model <- jags.model(file = textConnection(RandomWalkLogNormal),
                        data = data,
                        # inits = init,
                        n.chains = 3)
  
  jags.out2 <- coda.samples(
    model = j.model,
    variable.names = c("x", "tau_add"),
    n.iter = 15000
  )
  
  ## check burnin
  GBR2 <- gelman.plot(jags.out, ask = FALSE)
  burnin2 <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  cat("Burnin after", burnin2, "iterations\n")
  out2 <- window(jags.out2, start = burnin)
  out.mat2 <- as.matrix(out2)
  out.predict2 <- out.mat2[, -1]
  
  # save process error
  tau.add2[[i]] <- out.tau_add2 <- out.mat2[, "tau_add"]
  
  ## state CIs
  ci2[[i]] <- apply(out.predict2, 2, quantile, c(0.025, 0.5, 0.975))
  

# Plot results -----------------------------------

par(mfrow=c(2,1))
  
  # subset data 
  data <- dat %>%
    filter(Site == sites[i])
  
  y <- data$totalNymph # all data, for plotting
  
  time.axis <- lubridate::ymd(data$Day) # convert axis to date
  
  plot(time.axis, y, 
       pch = "", 
       main = paste(sites[i], "Nymphs"), 
       ylab = "Individuals", 
       xlab = "", ylim=c(0,range(y)[2]+25))

  
  ciEnvelope(time.axis, ci[[i]][1,], ci[[i]][3,], col = "lightblue")
  lines(time.axis, ci[[i]][2,], col = 'blue')
  points(time.axis, y, pch = 16)

  ciEnvelope(time.axis, ci2[[i]][1,], ci2[[i]][3,], col = adjustcolor("red",alpha.f=0.25) )
  lines(time.axis, ci2[[i]][2,], col = 'red')
 
  legend.x <- lubridate::ymd(data$Day[1]) # convert axis to date
  
  legend(x=legend.x,y=150,
         c("Normal Process","Log-Normal Process"),
         col=c("blue","red"),lty=c(1,1),
         cex=.5)
  
  plot(time.axis, y, 
       pch = "", 
       main = paste(sites[i], "Nymphs"), 
       ylab = "Individuals", 
       xlab = "", ylim=c(0,2000))
  
  
  ciEnvelope(time.axis, ci[[i]][1,], ci[[i]][3,], col = "lightblue")
  lines(time.axis, ci[[i]][2,], col = 'blue')
  points(time.axis, y, pch = 16)
  
  ciEnvelope(time.axis, ci2[[i]][1,], ci2[[i]][3,], col = adjustcolor("red",alpha.f=0.25) )
  lines(time.axis, ci2[[i]][2,], col = 'red')
  
  legend.x <- lubridate::ymd(data$Day[1]) # convert axis to date
  
  legend(x=legend.x,y=1500,
         c("Normal Process","Log-Normal Process"),
         col=c("blue","red"),lty=c(1,1),
         cex=.5)
  
  