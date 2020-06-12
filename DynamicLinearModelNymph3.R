library(ecoforecastR)
library(tidyverse)


dat <- read.csv("nymphDataWithMet.csv") # nymph data
sites <- unique(dat$Site) # sites vector
tau.add <- ci <- list() # storage

#for (i in seq_along(sites)) {

i = 1

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
uncert[is.na(uncert)] <- mean(uncert, na.rm = T)
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
      x[t] ~ dnorm(mu[t],tau_add) #T(0,)
      log(mu[t]) <- x[t-1] + beta*rho[t-1]
    }

    #### Priors
    x[1] ~ dpois(x_ic)
    beta ~ dnorm(0, 0.01)
    tau_add ~ dgamma(a_add,r_add)
  }"

data <- list(
  y = y,
  n = n,
  temp = temp,
  tau_temp = 1/uncert^2,
  x_ic = round(runif(1, 1, 15)),
  a_add = 0.001,
  r_add = 0.001,
  rho_mu = mean(temp, na.rm = TRUE),
  tau_rho = 0.001
)

init <- list(x = y)

j.model <- jags.model(file = textConnection(DynamicLinearModel),
                      data = data,
                      inits = init,
                      n.chains = 3)

jags.out <- coda.samples(
  model = j.model,
  variable.names = c("tau_add","x","mu","beta"),
  n.iter = 15000
)

## split output
out <- list(params = NULL, predict = NULL)
mfit <- as.matrix(jags.out, chains = TRUE)
pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
chain.col <- which(colnames(mfit) == "CHAIN")
out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])

## check burnin
# plot(jags.out)
GBR <- gelman.plot(out$params, ask = FALSE)
burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
cat("Burnin after", burnin, "iterations\n")
out.burn <- window(jags.out)
out.mat <- as.matrix(out.burn)
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

##########

# define IC

IC <- as.matrix(out$predict[,"x[43]"]) # initial conditions (posterior of X) in final year of training data
# if not interested in partitioning initial condition uncertainty, just take the mean of IC instead of using the full posterior
Nmc <- 20
NT <- nrow(data.val)

# mu <- exp(as.matrix(mfit[,"mu[43]"]))

forecastN <- function(IC, beta, rho, Q=tau_add, n=Nmc, NT){
  N <- matrix(NA,n,NT+1)  ## storage: NT = number of time steps forecasted, n = Nmc = number of ensembles
  Nprev <- IC           ## initialize
  for(t in 2:(NT+1)){
    
    mu = (Nprev + beta*rho[t-1])
    N[,t] <- rnorm(n,(mu),Q)                         ## predict next step
    N[,t] <- pmax(N[,t], 0)
    Nprev <- (N[,t] )                                 ## update IC
  }
  return(N)
}

# ndet
N.deterministic <- forecastN(IC=mean(IC),
                             beta=mean(out.mat[,"beta"]),
                             rho=(data.val$dailyMaxAirTemp-mean(data.full$dailyMaxAirTemp, na.rm=TRUE)),
                             Q=0,
                             NT=NT,
                             n=1)

ylim = c(0,150)  ## set Y range on plot
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy
ci <- apply(as.matrix(out$predict),2,quantile,c(0.025,0.5,0.975))
time = 1:(nrow(data.train)+nrow(data.val))    ## total time
time1 = 1:nrow(data.train)       ## calibration period
time2 = (max(time1))+(1:nrow(data.val))   ## forecast period

plot.run <- function(){
  plot(time,time,type='n',ylim=ylim,ylab="N")
  ecoforecastR::ciEnvelope(time1,ci[1,(time1)],ci[3,time1],col=col.alpha("lightBlue",0.6))
  lines(time1,ci[2,time1],col="blue")
  points(time1,y[time1])
}

# ### settings
# # Nmc = 1000         ## set number of Monte Carlo draws
#  ylim = c(0,100)  ## set Y range on plot
# N.cols <- c("black","red","green","blue","orange") ## set colors
# trans <- 0.8       ## set transparancy
# time = 1:(NT*2)    ## total time
# time1 = 1:NT       ## calibration period
# time2 = time1+NT   ## forecast period


## Plot run
plot.run()
lines(time2,N.deterministic[-1],col="purple",lwd=3)

## IC conditions
## sample parameter rows from previous analysis
params = out.mat[,"beta"]

prow = sample(params,Nmc,replace=TRUE)

N.I <- forecastN(IC=IC[prow],
                             beta=mean(out.mat[,"beta"]),
                             rho=(data.val$dailyMaxAirTemp-mean(data.full$dailyMaxAirTemp, na.rm=TRUE)),
                             Q=0,
                             NT=NT,
                             n=Nmc)

## Plot run
plot.run()
N.I.ci = apply(N.I[,-1],2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

#points(time2,data.val$totalNymph)

# Parameter uncertainty
N.IP <- forecastN(IC=IC[prow],  ## sample IC
                  beta=out.mat[,"beta"][prow],
                  rho=(data.val$dailyMaxAirTemp-mean(data.full$dailyMaxAirTemp, na.rm=TRUE)),
                  Q=0,
                  NT=NT,
                  n=Nmc)


## Plot run
par(mfrow=c(1,1))
plot.run()
N.IP.ci = apply(N.IP[,-1],2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

## Process uncertainty
## process error samples
Qmc <- 1/sqrt(out.mat[,"tau_add"][prow])  ## convert from precision to standard deviation

N.IPDE <- forecastN(IC=IC[prow],  ## sample IC
                    beta=out.mat[,"beta"][prow],
                    rho=(data.val$dailyMaxAirTemp-mean(data.full$dailyMaxAirTemp, na.rm=TRUE)),
                    Q=Qmc,
                    NT=NT,
                    n=Nmc)

ylim=c(0,400)
## Plot run
plot.run()
N.IPDE.ci = apply(N.IPDE[,-1],2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IPDE.ci[1,],N.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
#ecoforecastR::ciEnvelope(time2,N.IPD.ci[1,],N.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)
