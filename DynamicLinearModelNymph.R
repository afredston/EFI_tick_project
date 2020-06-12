library(ecoforecastR)
library(tidyverse)


dat <- read.csv("nymphDataWithMet.csv") # nymph data
sites <- unique(dat$Site) # sites vector
tau.add <- ci <- list() # storage

#for (i in seq_along(sites)) {

  site = 1
  
  # subset data
  data.train <- dat %>%
    filter(trainValidate == "Train") %>%
    filter(Site == sites[i])
  
  data.val <- dat %>% 
    filter(trainValidate == "Validate") %>%
    filter(Site == sites[i])
  
  data.full <- dat %>% 
    filter(Site == sites[i])
  
  na.add <- nrow(data.val) # number of validation data points (number of NAs to add)
  
  y <- data.train$totalNymph
  y <- c(data.train$totalNymph, rep(NA, na.add)) # add NAs
  temp <- data.full$dailyMaxAirTemp
  temp_center <- temp-mean(temp,na.rm=TRUE)
  uncert <- data.full$dailyAvgExpUncert
  n <- length(y)
  
  DynamicLinearModel <- "
   model{

    #### Data Model
    for(t in 1:n){
      y[t] ~ dpois(x[t])
      temp[t] ~ dnorm(rho[t], tau_temp[t])
      rho[t] ~ dnorm(rho_mu, tau_rho)
    }

    #### Process Model
    for(t in 2:n){
      x[t] ~ dnorm(mu[t],tau_add) T(0,)
      log(mu[t]) <- x[t-1] + beta*rho[t-1]
    }

    #### Priors
    x[1] ~ dpois(x_ic)
    beta ~ dnorm(0, 0.01)
    tau_add ~ dgamma(a_add,r_add)
  }"
  
  data <- list(
    y = y,
    temp = temp_center,
    n = n,
    x_ic = round(runif(1, 1, 15)),
    a_add = 0.001,
    r_add = 0.001,
    tau_p = 0.0001,
    rho = rep(0,length(temp_center)),
    tau_temp = 1/(rep(mean(uncert,na.rm=TRUE),length(uncert))^2)
  )
  
  nchain = 3
  init <- list()
  for(i in 1:nchain){
    y.samp = sample(y,length(y),replace=TRUE)
    init[[i]] <- list(tau_add=1/var(diff(log(y.samp))), x = y)
  }
  
  
  j.model <- jags.model(file = textConnection(DynamicLinearModel),
                        data = data,
                        inits = init,
                        n.chains = 3)
  
  jags.out <- coda.samples(
    model = j.model,
    variable.names = c("tau_add", "beta"),
    n.iter = 15000
  )
  
  ## check burnin
  plot(jags.out)
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


# plot
par(mfrow = c(2,2))
for(i in seq_along(sites)){
  
  # subset data 
  data <- dat %>%
    filter(Site == sites[i])
  
  y <- data$totalNymph # all data, for plotting
  
  time.axis <- lubridate::ymd(data$Day) # convert axis to date
  
  plot(time.axis, y, 
       pch = "", 
       main = paste(sites[i], "Nymphs"), 
       ylab = "Individuals", 
       xlab = "")
  ciEnvelope(time.axis, ci[[i]][1,], ci[[i]][3,], col = "lightblue")
  lines(time.axis, ci[[i]][2,])
  points(time.axis, y, pch = 16)
}




