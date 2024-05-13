## -----------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
theme_set(theme_light())

tr <- seq(15, 85, length = 10)
traps <- data.frame(X = rep(tr, each = length(tr)),
                    Y = rep(tr, times = length(tr))) # 100 coord. traps

set.seed(10)
xlim <- c(0, 100)
ylim <- c(0, 100) # area 100 * 100 = 1e4

N <- 50 # generate population

s <- data.frame(s.x = runif(N, xlim[1], xlim[2]), 
                s.y = runif(N, ylim[1], ylim[2]))

sigma <- 5
lambda0 <- 0.4
J <- nrow(traps) # nb of traps
K <- 5 # nb capture occasions
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((traps$X[j] - s$s.x)^2 + (traps$Y[j] - s$s.y)^2)
  lambda <- lambda0 * exp(-dist^2 / (2 * sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}
n <- apply(yy, c(1,2), sum)


## -----------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
theme_set(theme_light())
viz_traps <- traps %>% 
  ggplot(aes(x = X, y = Y)) +
  geom_point(pch = 3, size=4) + 
  xlim(0, 100) +
  ylim(0, 100)
viz_traps

## ----------------------------------------------------------------------------------------------------------------------------------

viz_traps_ac <- viz_traps + 
  geom_point(data = s, aes(x = s.x, y = s.y), pch = 16, size =4, color = "red")
viz_traps_ac

## -----------------------------------------------------------------------------------------------------------------------------------
n 

nobs <- n[apply(n,1,sum)!=0,] # remove undetected individuals
nobs


## -----------------------------------------------------------------------------------------------------------------------------------
tot <- apply(n, 1, sum)
dat <- data.frame(traps, tot = tot)

viz_traps_ac +
  geom_point(data = dat, aes(x = X, y = Y, size = tot), alpha = 0.3) +
  scale_size(range = c(0, 20)) +
  labs(x = "",
       y = "",
       size = "# detections")


## -----------------------------------------------------------------------------------------------------------------------------------
 library(nimble)


## -----------------------------------------------------------------------------------------------------------------------------------
scr <- nimbleCode({
 
   # priors
   sigma ~ dunif(0, 10000)
   p0 ~ dunif(0, 1)
   psi ~ dunif(0, 1)
 
   for(i in 1:M) { # for each possible individuak
 
     z[i] ~ dbern(psi) # decide if it is real or ghost
 
     # assign ACs to individual i
     s[i,1] ~ dunif(xlim[1], xlim[2])
     s[i,2] ~ dunif(ylim[1], ylim[2])
 
     # calculate distance btw AC and each trap
     dist[i,1:J] <- sqrt((s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2)
     # calculate detection probability as an half-normal function
     p[i,1:J] <- p0*exp(-dist[i,1:J]^2 / (2 * sigma^2)) * z[i]
 
   # likelihood
   for(j in 1:J){ # for each trap
       n[i,j] ~ dpois(p[i,j]*K)
   }
   }
   # calculate population size as a derived parameter
   N <- sum(z[1:M])
 })

## -----------------------------------------------------------------------------------------------------------------------------------
## # data augmentation parameters
 M <- 200
 
 # a matrix of 0 for potential individuals
 n0 <- matrix(0,nrow = M-nrow(nobs), ncol = ncol(nobs))
 ndata <- rbind(nobs, n0)
 
 # constants
 constants <- list(M = M,
                   K = K,
                   J = J)
 # data
 data <- list(n = ndata,
              X = traps,
              xlim = xlim,
              ylim = ylim)
 
 # initial values
 s <- cbind(runif(M, xlim[1], xlim[2]),
            runif(M, ylim[1], ylim[2]))
 z <- rep(1, M)
 
 inits <- list(sigma = 1000,
               p0 = 0.6,
               s = s,
               z = z,
               psi = 0.5)
 


## ----All in one -----------------------------------------------------------------------------------------------------
## samples <- nimbleMCMC(code = scr,
##            data = data,
##            constants = constants,
##            inits = inits,
##            monitors =  c("N", "sigma", "p0"),
##            niter = 5000,
##            nburnin = 1000,
##            nchains = 2)
## 


## -------Step by step -----------------------------------------------------------------------------------

 # Approximate Running Time (ART) 1min
 Rmodel <- nimbleModel(code = scr,
                      constants = constants,
                      data = data,
                      inits = inits)

Rmodel$initializeInfo()
Rmodel$calculate() # - 61769.03

Cmodel <- compileNimble(Rmodel)

calculate(Cmodel)

conf <- configureMCMC(Rmodel,
                      monitors = c("N", "sigma", "p0"))

Rmcmc <- buildMCMC(conf)

Cmcmc <- compileNimble(Rmcmc, project = Cmodel)


# Approximate Running Time 5 min
niter = 5000
nburnin = 1000
nchains = 2

samples <- runMCMC(Cmcmc,
                   niter = niter,
                   nburnin = nburnin,
                   nchains = nchains)



## -----------------------------------------------------------------------------------------------------------------------------------

# load results
load("scr_simulation.rdata")


library(MCMCvis)
MCMCtrace(samples, pdf = F)

MCMCsummary(samples, round = 2)



## -----------------------------------------------------------------------------------------------------------------------------------
samplesBind <- rbind(samples$chain1[1001:(niter- nburnin),],
                     samples$chain2[1001:(niter- nburnin),])


# do the burnin and plot
samplesBind %>% 
  as_tibble() %>% 
  mutate(chain = rep(c("chain1","chain2"), each = length(1001:(niter- nburnin)))) %>% 
  ggplot() + 
  geom_density(aes(x = N, fill = chain), alpha = 0.4)+
  labs(title= "Posterior distribution of N, the population size",
       subtitle = "SCR model")+
  theme(text = element_text(face = "bold", color = "black"),
        axis.text.y = element_blank(),
        axis.text.x =  element_text(face = "bold", color = "black", size = 14))


## -------Reconfig to save activity centers----------------------------------------------------------------------------------------------
conf2 <- configureMCMC(Rmodel,
                       monitors = c("N", "sigma", "p0", "s", "z"))

Rmcmc2 <- buildMCMC(conf2)

Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel)


 # ART 5 min

samples2 <- runMCMC(Cmcmc2,
                        niter = 5000,
                        nburnin = 1000,
                        nchains = 2)



## -----------------------------------------------------------------------------------------------------------------------------------
str(samples2)

