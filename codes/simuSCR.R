## -----Simulaciones de SCR data ----------------------------------------------------------------------------------------

library(tidyverse)
library(sf)
theme_set(theme_light())


set.seed(10)
xlim <- c(0, 100)
ylim <- c(0, 100) # area 100 * 100 = 1e4

# two points to design a covariate named "distance to bar"
bar <- tibble(x = c(50, 60), y = c(90, 40)) %>% 
  st_as_sf(coords = c("x","y"))

set.seed(10)
xlim <- c(0, 100)
ylim <- c(0, 100) # area 100 * 100 = 1e4

# parameter that direct the intensity regression mu = m0 + m1 * dbar
m0 <- - 1
m1 <- 0.8

# make grid of the study area
ggrid <- sf::st_make_grid(cellsize = 10, offset = c(0,0), n = c(20)) %>% 
  st_as_sf() 

# create the covariate
ggrid <- ggrid%>% 
  mutate(dbar = st_distance(ggrid ,bar) %>% apply(1,min),
         dbar.sc = scale(dbar)[,1])

# calculate the intensity
mu <- ggrid %>% 
  mutate(lam = exp(m0 + m1 * dbar.sc))

# Randomly select 50 sampled sites among all grid-cells
nsamp <- sample(1:nrow(ggrid), 50, replace = F)

# traps
traps <- ggrid[nsamp,] %>%
  st_centroid() %>% 
  st_coordinates()

# sites
sites <- ggrid %>%
  st_centroid() %>% 
  st_coordinates()

# check the map
ggplot() +  geom_sf(data = ggrid, aes(fill = dbar)) + 
  geom_sf(data = ggrid[nsamp,], color = "yellow3", lwd = 2, alpha=0)+
  geom_sf(data = bar, color = "red", size = 4)

# define population size as the sum of the density in each grid-cell
N <- rpois(1,sum(mu$lam))

# define some 
idd <- rmultinom(N, 1, mu$lam)
id <- rep(0,N)
for(i in 1:N){
  id[i] <- which(idd[,i]==1)
}
id

s <- data.frame(s.x = sites[id,1], 
                s.y = sites[id,2])

sigma <- 5
lambda0 <- 0.4
J <- length(nsamp) # nb of traps
K <- 5 # nb capture occasions
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((traps[j,"X"] - s$s.x)^2 + (traps[j,"Y"] - s$s.y)^2)
  lambda <- lambda0 * exp(-dist^2 / (2 * sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}
n <- apply(yy, c(1,2), sum)


## ---------plot traps ----------------------------------------------------------------------
library(tidyverse)
theme_set(theme_light())
viz_traps <- ggrid[nsamp,] %>%
  st_centroid() %>% 
  ggplot() +
  geom_sf(color = "darkgreen", size = 3) 
viz_traps

## ---------plot traps and activity centers ---------------------------------------------------------------------------------
ac <- ggrid[id,] %>% st_centroid()
viz_traps_ac <- viz_traps + 
  geom_sf(data =ac, color = "red")
viz_traps_ac

## -----------------------------------------------------------------------------------------------------------------------------------
n 

nobs <- n[apply(n,1,sum)!=0,] # remove undetected individuals
nobs


## ----------------- plot with detections ------------------------------------------------------------------
tot <- apply(n, 2, sum)
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
 M <- 500
 
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
 samples <- nimbleMCMC(code = scr,
            data = data,
            constants = constants,
            inits = inits,
            monitors =  c("N", "sigma", "p0"),
            niter = 3000,
            nburnin = 1000,
            nchains = 2)
 
 ## -----------------------------------------------------------------------------------------------------------------------------------
 
 # load results if you dont want to run it
 # load("scr_simulation.rdata")
 # results are here https://github.com/valentinlauret/diberticurso/tree/main/codes
 
 
 library(MCMCvis)
 MCMCtrace(samples, pdf = F)
 
 MCMCsummary(samples, round = 2)
 
 
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
 
 

# ----------- SCR WITH DENSITY ----------------

scr.ipp <-  nimbleCode({
  
  # priors
  b0 ~ dnorm(0,1) # Intercept of lambda-density regression
  b1 ~ dnorm(0,1) # slope of bathymetry
  
  # priors for the sigma 
  sigma ~ dunif(10,10^5)
  sig2 <- 2*sigma*sigma
  
  # priors for the p0
  p0 ~ dunif(0,1) 
  
  # Intensity process
  log(mu[1:nsites]) <- b0 + b1 * dbar[1:nsites]  # abundance model via IPP 
  
  # Population size
  EN <- sum(mu[1:nsites]) # expected number of individual for the study area
  
  #psi <-   EN/M #
  psi ~ dunif(0,1) # parameter for data augmentation

  # cell prob
  probs[1:nsites] <- mu[1:nsites]/EN
  
  # data augmentation
  for(i in 1:M) {
    
    z[i] ~ dbern(psi)  # is i alive or not ?
    
    # induce the uniform distribution of activity centers over the discrete set 1:nsites
    id[i] ~ dunif(0.5, nsites+0.5) 
    
    # location of AC
    SX[i] <- sites[round(id[i]),1]
    SY[i] <- sites[round(id[i]),2]
    
    #zeros's trick
    nll[i] <- - log(probs[round(id[i])])
    zeros[i] ~ dpois(nll[i])
    
    dist[i, 1:ntrap] <- sqrt((traps[1:ntrap,1] - SX[i])^2 + (traps[1:ntrap,2] - SY[i])^2)
    
    p[i, 1:ntrap] <- p0 * exp( - dist[i, 1:ntrap]^2 / sig2) * z[i]
    
    for(j in 1:ntrap){
      y.scr[i, j] ~ dpois(p[i,j]*K)
    }
    
  }
  N <- sum(z[1:M]) 
})


## # data augmentation parameters
M <- 500

# a matrix of 0 for potential individuals
n0 <- matrix(0,nrow = M-nrow(nobs), ncol = ncol(nobs))
ndata <- rbind(nobs, n0)

nsites = nrow(ggrid)
ntrap = length(nsamp)

# constants
constants <- list(M = M,
                  K = K,
                  ntrap = ntrap,
                  nsites = nsites,
                  traps = traps,
                  dbar = ggrid$dbar.sc)
# data
data <- list(y.scr = ndata,
             sites = ggrid %>% st_centroid() %>%  st_coordinates(),
             zeros = rep(0,M))


# initial values
id <- round(runif(M, 0,nsites))


z <- apply(ndata, 1, max)
z[z>0] <- 1 

inits <- list(b0 = rnorm(1, 0, 1),
              b1 = rnorm(1, 0, 1),
              psi = runif(1,0,1),
              sigma = 1000,
              p0 = 0.6,
              id = id,
              z = z,
              psi = 0.5)


# Approximate Running Time (ART) 1min
Rmodel <- nimbleModel(code = scr.ipp,
                      constants = constants,
                      data = data,
                      inits = inits)

Rmodel$initializeInfo()
Rmodel$calculate() # 

Cmodel <- compileNimble(Rmodel)

calculate(Cmodel)

conf <- configureMCMC(Rmodel,
                      monitors = c("N", "sigma", "p0"))

Rmcmc <- buildMCMC(conf)

Cmcmc <- compileNimble(Rmcmc, project = Cmodel)


# Approximate Running Time 5 min
niter = 3000
nburnin = 1000
nchains = 2

samples2 <- runMCMC(Cmcmc,
                   niter = niter,
                   nburnin = nburnin,
                   nchains = nchains)

# check SCR density
library(MCMCvis)
MCMCtrace(samples2, pdf = F)

MCMCsummary(samples2, round = 2)


samplesBind <- rbind(samples$chain1,
                     samples$chain2)

