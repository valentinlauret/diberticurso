## ---------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(nimble)


## ---------------------------------------------------------------------------------------------------------------------
# Generate count data
set.seed(1)
temp <- round(runif(100,-2,10),1)
lam <- exp(1 + 0.5*temp)
conteo <- rpois(100, lam)
datos <- data.frame(conteo = conteo, temp = temp)
datos

## ---------------------------------------------------------------------------------------------------------------------
# Fit frequentist glm
mod <- glm(conteo ~ temp, data = datos, family = poisson(link = "log"))



## ---------------------------------------------------------------------------------------------------------------------
 
# nimble model
glm.wb <- nimbleCode({
 
   # priors
   b0 ~ dnorm(0, 1)
   b1 ~ dnorm(0, 1)
 
   # likelihood
   for(i in 1:nsites){
 
     log(lambda[i]) <- b0 + b1 * temp[i]
 
     conteo[i] ~ dpois(lambda[i])
   }
 
 })
 
# read in data
 my.data <- list(conteo = datos$conteo)
 my.constants <- list(temp = datos$temp,
                      nsites = nrow(datos))
 
# initial values 
 initial.values <- list(b0 = rnorm(1,0,1),
                        b1 = rnorm(1, 0, 1))
 
# parameters to save 
 parameters.to.save <- c("b0", "b1")
 
# MCMC settings
 n.iter <- 5000
 n.burnin <- 1000
 n.chains <- 2
 n.thin <- 1


# ---------------------------------------------------------------------------------------------------------------------
# Run with Nimble
  mcmc.output <- nimbleMCMC(code = glm.wb,
                           data = my.data,
                           constants = my.constants,
                           inits = initial.values,
                           monitors = parameters.to.save,
                           thin = n.thin,
                           niter = n.iter,
                           nburnin = n.burnin,
                           nchains = n.chains)


## ------- View results and check convergence -------------------------------------------------------------------------------------
str(mcmc.output)

head(mcmc.output$chain1)

library(MCMCvis)
MCMCsummary(mcmc.output, round = 2)

# compare with glm frequentist
summary(mod)

# check convergence
MCMCtrace(mcmc.output,  pdf = F) 



## ----- Prediction --------------------------------------------------------------------------------------

mcmc.bind <- rbind(mcmc.output$chain1, mcmc.output$chain2)


mcmc.bind

# the gradient of temperature we want to predict on
pred.temp <- seq(from = -5, to = 9, by = 0.1)

pred.lambda <- matrix(NA, nrow = nrow(mcmc.bind), ncol = length(pred.temp))

for(i in 1:nrow(pred.lambda)){
  pred.lambda[i,] <- exp(mcmc.bind[i,"b0"] + mcmc.bind[i,"b1"] * pred.temp)
}


mean.lambda <- apply(pred.lambda, 2, mean)
sd.lambda <- apply(pred.lambda, 2, sd)
ci <- apply(pred.lambda, 2, quantile, c(0.1, 0.9))

pred.results <- tibble(temp = pred.temp, 
                       mean = mean.lambda,
                       sd = sd.lambda,
                       cinf = ci[1,],
                       csup = ci[2,])


mean.lambda <- apply(pred.lambda, 2, mean) # extract mean
sd.lambda <- apply(pred.lambda, 2, sd) # extract sd
ci <- apply(pred.lambda, 2, quantile, c(0.1, 0.9)) # 80% CI

pred.results <- tibble(temp = pred.temp, 
                       mean = mean.lambda,
                       sd = sd.lambda,
                       cinf = ci[1,],
                       csup = ci[2,])

pred.results %>% 
  ggplot(aes(x = temp, y = mean), color = "#000A39") +
  geom_line() + 
  geom_line(aes(x = temp, y= cinf), linetype = 2, color = "#D2A34E" )+
  geom_line(aes(x = temp, y= csup), linetype = 2, color = "#D2A34E" )+
  labs(title = "Number of wild boar poop (WBP) as a function of temperature",
       subtitle = paste0("b0 = ", round(mean(mcmc.bind[,"b0"]),2),
                         " & b1 = ",  round(mean(mcmc.bind[,"b1"]),2)),
       x = "Temperature",
       y = "Nb of WBP") + 
  theme_minimal()+
  theme(text = element_text(face = "bold", family = "Comic Sans MS", color ="#840032"),
        plot.subtitle =  element_text(family = "Times"))

