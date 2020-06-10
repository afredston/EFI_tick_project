library(ecoforecastR)
library(tidyverse)


dat <- read.csv("nymphDataWithMet.csv") # nymph data
sites <- unique(dat$Site) # sites vector
tau.add <- ci <- list() # storage

for (i in seq_along(sites)) {
  
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
  
}

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




