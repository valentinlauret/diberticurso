# Marios y Albertos del Womad 

# Simulaciones


tr <- seq(20, 80, length = 10)
traps <- data.frame(X = rep(tr, each = length(tr)),
                    Y = rep(tr, times = length(tr))) # 100 coord. traps

bar <- tibble(x = c(50, 60), y = c(90, 40)) %>% 
  st_as_sf(coords = c("x","y"))

set.seed(10)
xlim <- c(0, 100)
ylim <- c(0, 100) # area 100 * 100 = 1e4

m0 <- - 1
m1 <- 0.02
a0 <- - 0.5
a1 <- -0.9

# make grid
ggrid <- sf::st_make_grid(cellsize = 10, offset = c(0,0), n = c(10)) %>% 
  st_as_sf() 

ggrid <- ggrid%>% 
  mutate(dbar = st_distance(ggrid ,bar) %>% apply(1,min),
         dbar.sc = scale(dbar)[,1])

lam_mario <- ggrid %>% 
  mutate(lam = exp(m0 + m1 * dbar.sc))

lam_alberto <- ggrid %>% 
  mutate(lam = exp(a0 + a1 * dbar.sc))

# Nmario <- 162
# Nalberto <- 300 # generate population

# sampled sites
 nsamp <- sample(1:100, 10, replace = F)
 
 ggplot() +  geom_sf(data = ggrid, aes(fill = dbar)) + 
   geom_sf(data = ggrid[nsamp,], color = "yellow3", lwd = 2, alpha=0)+
   geom_sf(data = bar, color = "red", size = 4)+ 
   xlim(0, 100) +
   ylim(0, 100)
 
mobs <- lam_mario[nsamp,] %>% 
  mutate(nmario = rpois(length(nsamp), lam))
aobs <- lam_alberto[nsamp,] %>% 
  mutate(nalberto = rpois(length(nsamp), lam))
 

## Bayesian model
mariosyalbertos <- nimbleCode({
  
  # priors
  m0 ~ dnorm(0,1)
  m1 ~ dnorm(0,1)
  a0 ~ dnorm(0,1)
  a1 ~ dnorm(0,1)
  
  # likelihood
  for(i in 1:nsites){
    
    # lambda
    log(lamario[i]) <-  m0 + m1 * distToBar[i]
    log(lalberto[i]) <-  a0 + a1 * distToBar[i]
    
    
    marios[i] ~ dpois(lamario[i])
    albertos[i]  ~ dpois(lalberto[i])
    
  }
 
  
})

data <- list(marios = mobs$nmario,
             albertos = aobs$nalberto)

constants <- list(nsites = length(nsamp),
                  distToBar = mobs$dbar.sc)

inits <- list(m0 = 1,
              m1 = 0,
              a0 = 1,
              a1 = 0)

samples <- nimbleMCMC(code = mariosyalbertos,
                      data = data,
                      constants = constants,
                      inits = inits,
                      monitors =  c("a0", "a1", "m0", "m1"),
                      niter = 1000,
                      nburnin = 100,
                      nchains = 2)

library(MCMCvis)

MCMCtrace(samples, pdf = F)
MCMCsummary(samples, round =2)

str(samples)

samplesbind <- rbind(samples$chain1, samples$chain2)

pred.dbar <- ggrid$dbar.sc

pred.alberto <- pred.mario <- matrix(NA, nrow = nrow(samplesbind), ncol = length(pred.dbar))
for(i in 1:nrow(pred.alberto)){
  pred.alberto[i,] <- exp(samplesbind[i,"a0"] + samplesbind[i,"a1"] * pred.dbar)
  pred.mario[i,]<- exp(samplesbind[i,"m0"] + samplesbind[i,"m1"] * pred.dbar)
}

palberto <- apply(pred.alberto,2,mean)
pmario <- apply(pred.mario,2,mean)

ggplot() + geom_sf(data = ggrid, aes(fill = palberto))
ggplot() + geom_sf(data = ggrid, aes(fill = pmario)) +
  geom_sf(data = bar, size = 4, color = 'red')

