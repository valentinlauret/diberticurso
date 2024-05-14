

## -----------------------------------------------------------------------------------------------------------------------
# load data
load("data_womad.rdata")

# load packages
library(nimble)
library(tidyverse)

# conteos de marios 
nmarios

# conteos de albertos
nalbertos

# grid of the womad
ggrid















####### REPUESTAS #####

## -----------------------------------------------------------------------------------------------------------------------
library(nimble)

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


## -----------------------------------------------------------------------------------------------------------------------
data <- list(marios = nmarios$nmario,
             albertos = nalbertos$nalberto)

constants <- list(nsites = length(nsamp),
                  distToBar = nmarios$dbar.sc)

# initial values
inits <- list(m0 = 1,
              m1 = 0,
              a0 = 1,
              a1 = 0)


## -----------------------------------------------------------------------------------------------------------------------

# lets run Nimble
samples <- nimbleMCMC(code = mariosyalbertos,
                      data = data,
                      constants = constants,
                      inits = inits,
                      monitors =  c("a0", "a1", "m0", "m1"),
                      niter = 1000,
                      nburnin = 100,
                      nchains = 2)


## -----------------------------------------------------------------------------------------------------------------------
# check convergence
library(MCMCvis)

MCMCtrace(samples, pdf = F)


## -----------------------------------------------------------------------------------------------------------------------
MCMCsummary(samples, round =2)


## -----------------------------------------------------------------------------------------------------------------------
str(samples)


## -----------------------------------------------------------------------------------------------------------------------
# juntar las dos cadenas
samplesbind <- rbind(samples$chain1, samples$chain2)

# extraer la valores para predecir
pred.dbar <- ggrid$dbar.sc

# preallocate memory
pred.alberto <- pred.mario <- matrix(NA, nrow = nrow(samplesbind), ncol = length(pred.dbar))

# por cada iteracion i
for(i in 1:nrow(pred.alberto)){
  
  # predecimos la densidad lambda en cada celda pred.bar
  pred.alberto[i,] <- exp(samplesbind[i,"a0"] + samplesbind[i,"a1"] * pred.dbar)
  pred.mario[i,]<- exp(samplesbind[i,"m0"] + samplesbind[i,"m1"] * pred.dbar)
  
}

# calculamos la media de cada celda por cada iteracion
palberto <- apply(pred.alberto,2,mean)
pmario <- apply(pred.mario,2,mean)


# Calculamos el tamano de poblacion de Marios y Albertos sumando la densidad de todas las celdas
Nmar <- sum(pmario)
Nalb <- sum(palberto)

# Plot of density maps of Marios and Albertos
plot_alb <- ggplot() + geom_sf(data = ggrid, aes(fill = palberto))+
  geom_sf(data = bar, size = 4, color = 'red')+
   labs(title = "Densidad de Albertos",
       fill = "Densidad") + 
  theme(plot.title = element_text(family = "Comic Sans MS", face = "bold", color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") + 
  paletteer::scale_fill_paletteer_c(palette = "grDevices::RdYlBu", direction = -1)+
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5,
                                barwidth = unit(10, 'lines'), barheight = unit(.5, 'lines')),
         colour = "none") 

plot_mar <- ggplot() + geom_sf(data = ggrid, aes(fill = pmario)) +
  geom_sf(data = bar, size = 4, color = 'red') + 
   labs(title = "Densidad de Marios",
       fill = "Densidad") + 
  theme(plot.title = element_text(family = "Comic Sans MS", face = "bold", color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") + 
  paletteer::scale_fill_paletteer_c(palette = "grDevices::RdYlBu", direction = -1)+
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5,
                                barwidth = unit(10, 'lines'), barheight = unit(.5, 'lines')),
         colour = "none") 

library(cowplot)

plot_grid(plot_alb, plot_mar)



## -----------------------------------------------------------------------------------------------------------------------
Nmar


## -----------------------------------------------------------------------------------------------------------------------
Nalb


### Simulaciones de datos

library(tidyverse)
library(sf)

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
nsamp <- sample(1:100, 20, replace = F)


nmarios <- lam_mario[nsamp,] %>% 
  mutate(nmario = rpois(length(nsamp), lam))%>% 
  st_drop_geometry() %>% 
  select(dbar, dbar.sc, nmario)
nalbertos <- lam_alberto[nsamp,] %>% 
  mutate(nalberto = rpois(length(nsamp), lam)) %>% 
  st_drop_geometry() %>% 
  select(dbar, dbar.sc, nalberto)

# save(nmarios, nalbertos, ggrid, file= "data_womad.rdata")


## ----Map of Womad study area"----------------------------------------------------------------------------------


ggplot() +  geom_sf(data = ggrid, aes(fill = dbar)) + 
  geom_sf(data = ggrid[nsamp,], color = "yellow3", lwd = 2, alpha=0)+
  geom_sf(data = bar, color = "red", size = 4)+ 
  xlim(0, 100) +
  ylim(0, 100)+ 
  labs(title = "Zona de estudio del Womad Caceres",
       subtitle = "Puntos rojos representan los bares \n Celdas amarillas han estado visitado para el muestreo",
       fill = "Distancia hasta \n el bar mas cerca") + 
  theme(plot.title = element_text(family = "Comic Sans MS", face = "bold", color = "black"),
        legend.position = "bottom",
        plot.title.position = "plot") + 
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = .5,
                                barwidth = unit(10, 'lines'), barheight = unit(.5, 'lines')),
         colour = "none") 

