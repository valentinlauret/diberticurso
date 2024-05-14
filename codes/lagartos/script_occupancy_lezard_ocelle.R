

library(AICcmodavg)	# charger le package AICcmodavg
library(boot)


# ------------------------------------- #
# 	 Site occupancy         		#		----
# ------------------------------------- #


data <- read.csv2('codes/ocelle_singleseason.csv', h = T)	# abrir el dataframe

y <- data[, 2:4]				# extraer conteos
yp <- y					# transformarlos como presencia ausencia
yp[yp>0] <- 1
covsite <- data[, 5:16]			#  covariables de sitios
temp <- data[, 17:19]			# rcovariables de occasiones
vent <- data[, 20:22]
mois <- data[, 23:25]

# Con nimble ----
load(nimble)

# por ejemplo un modelo con conejos
occu <- nimbleCode({
  
  # priors
  a0 ~ dnorm(0,1)
  a1 ~ dnorm(0,1)
  a2 ~ dnorm(0,1)
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  
  
  # likelihood
  for(i in 1:nsites){
    
    z[i] ~ dbern(psi[i])
    
    logit(psi[i]) <- b0 + b1 * lapin[i]
    
    for(j in 1:nocc){
      
    logit(p[i,j]) <- a0 + a1 * temp[i,j] + a2 * vent[i,j]
    
    y[i,j] ~ dbern(z[i] * p[i,j])
      
    }
  }
  
})

data <- list(y = yp)
constants <- list(nsites = nrow(yp), 
                  nocc = ncol(yp),
                  lapin = covsite$terlapin,
                  vent = vent,
                  temp = temp)

inits <- list(a0 = rnorm(1,0,1),
              a1 = rnorm(1,0,1),
              a2 = rnorm(1,0,1),
              b0 = rnorm(1,0,1),
              b1 = rnorm(1,0,1),
              z = apply(yp,1,max))

samples <- nimbleMCMC(code = occu,
                      data = data,
                      constants = constants,
                      inits = inits,
                      monitors =  c("a0", "a1", "a2", "b0", "b1"),
                      niter = 1000,
                      nburnin = 100,
                      nchains = 2)

library(MCMCvis)
MCMCsummary(samples, round = 2)
MCMCtrace(samples, pdf =F)


# Correr un modelo con {unmarked} ----
library(unmarked)	# load  package unmarked

obscov <- list(temp = temp, vent = vent, mois = mois)	# hacer una lista con covariables de occasion

counto <- unmarkedFrameOccu(yp, siteCovs = covsite, obsCovs = obscov)	# funcion unmarked para occupancy
str(counto)
?unmarkedFrameOccu

summary(counto) # ver que hay dentro del objecto

res_null <- occu(~1 ~1, counto)		# ajustar un modelo nulo
res_null

psi <- backTransform(res_null, type = "state") # recuperar la valor de psi haciendo la tansformacion inverse-logit
psi
confint(psi)				   # et ses IC ? 95%	

p <- backTransform(res_null, type = "det") 	# ecuperar la valor de p haciendo la tansformacion inverse-logit
p
confint(p)				#  IC 95%



# probar otros modelos y calculando AIC
res_mois <- occu(~mois ~1, counto)

res_temp <- occu(~temp ~1, counto)

res_null@AIC
res_mois@AIC
res_temp@AIC


# predecir las valores
newData1 <- data.frame(mois = c('p1', 'p2', 'p3'))			
pred <- predict(res_mois, type = "det", newdata = newData1, appendData = TRUE) # hacer la prediccion


res_mois2 <- occu(~mois+temp ~1, counto)
newData2 <- data.frame(mois = rep('p2', 10), temp = seq(20, 25, length = 10))
pred <- predict(res_mois2, type = "det", newdata = newData2, appendData = TRUE) # prediccion

# seleccion de modelos 
res_mois_terrier <- occu(~mois ~terlapin, counto)
res_mois_terrier2 <- occu(~mois+terlapin ~1, counto)
res_mois_terrier3 <- occu(~mois+terlapin ~terlapin, counto)

res_mois_habitat <- occu(~mois ~habitat, counto)
res_mois_habitat2 <- occu(~mois+habitat ~1, counto)
res_mois_habitat3 <- occu(~mois+habitat ~habitat, counto)

newData2 <- data.frame(habitat = c("dune", "clairiere"))
E.p <- predict(res_mois_habitat, type = "state", newdata = newData2, appendData = TRUE)


# dataframe con todos los modelos
fms <- fitList("psi(.)p(.)" = res_null, "psi(.)p(mois)" = res_mois,
"psi(lapin)p(mois)" = res_mois_terrier, "psi(.)p(mois+lapin)" = res_mois_terrier2, "psi(lapin)p(mois+lapin)" = res_mois_terrier3,
"psi(habitat)p(mois)" = res_mois_habitat, "psi(.)p(mois+habitat)" = res_mois_habitat2, "psi(habitat)p(mois+habitat)" = res_mois_habitat3
)

# comparison de los modelos
ms <- modSel(fms)
ms

# predecir variables continuas
newData2 <- data.frame(terlapin = seq(0, 20, by = 1))
E.p <- predict(res_mois_terrier, type = "state", newdata = newData2, appendData = TRUE)
plot(Predicted ~ terlapin, E.p, type = "l", ylim = c(0,1),col = "darkred",lwd = 3,
xlab = "nb terriers de lapin",
ylab = "Expected occupancy probability")
lines(lower ~ terlapin, E.p, type = "l", col = gray(0.5))
lines(upper ~ terlapin, E.p, type = "l", col = gray(0.5))



# ------------------------------------- #


# 	MODELOS de N-mixture		#	----	
# ------------------------------------- #
data <- read.csv2('ocelle_singleseason.csv', h = T)	# load data

y <- data[, 2:4]				# recuperar conteos
covsite <- data[, 5:16]			# recuperar covariables de sitios
temp <- data[, 17:19]			  # recuperar covariables de ocasiones
vent <- data[, 20:22]
mois <- data[, 23:25]

# Bayesian Nmix ---- -Con nimble-----

load(nimble)

# por ejemplo un modelo con conejos
nmix <- nimbleCode({
  
  # priors
  a0 ~ dnorm(0,1)
  a1 ~ dnorm(0,1)
  a2 ~ dnorm(0,1)
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  
  
  # likelihood
  for(i in 1:nsites){
    
    n[i] ~ dpois(lam[i])
    
    log(lam[i]) <- b0 + b1 * lapin[i]
    
    for(j in 1:nocc){
      
      logit(p[i,j]) <- a0 + a1 * temp[i,j] + a2 * vent[i,j]
      
      y[i,j] ~ dbinom(p[i,j],n[i])
      
    }
  }
  
  Npop <- sum(n[1:nsites])
  
})

data <- list(y = y)
constants <- list(nsites = nrow(y), 
                  nocc = ncol(y),
                  lapin = covsite$terlapin,
                  vent = vent,
                  temp = temp)

inits <- list(a0 = rnorm(1,0,1),
              a1 = rnorm(1,0,1),
              a2 = rnorm(1,0,1),
              b0 = rnorm(1,0,1),
              b1 = rnorm(1,0,1),
              n = apply(y,1,max))

samplesNmix <- nimbleMCMC(code = nmix,
                      data = data,
                      constants = constants,
                      inits = inits,
                      monitors =  c("a0", "a1", "a2", "b0", "b1", "Npop"),
                      niter = 1000,
                      nburnin = 100,
                      nchains = 2)

library(MCMCvis)
MCMCsummary(samplesNmix, round = 2)
MCMCtrace(samplesNmix, pdf =F)

# # con package unmakred ----

load(unmarked)

obscov <- list(temp = temp, vent = vent, mois = mois)	# lista con  covariables de ocasiones

count <- unmarkedFramePCount(y, siteCovs = covsite, obsCovs = obscov) # objecto para el Nmix

res_nullP <- pcount(~1 ~1, count, mixture = c("P"))
res_nullNB <- pcount(~1 ~1, count, mixture = c("NB"))
res_nullZIP <- pcount(~1 ~1, count, mixture = c("ZIP"))
res_nullP@AIC
res_nullNB@AIC
res_nullZIP@AIC

res_nullPk105 <- pcount(~1 ~1, count, mixture = c("P"))
res_nullPk200 <- pcount(~1 ~1, count, K = 200, mixture = c("P"))

res_mois <- pcount(~mois ~1, count, mixture = c("P"))

newData1 <- data.frame(mois = c('p1', 'p2', 'p3'))
pred <- predict(res_mois, type = "det", newdata = newData1, appendData = TRUE)

lambda <- backTransform(res_mois, type = "state") 
lambda
confint(lambda)

res_mois_terrier <- pcount(~mois ~terlapin, count, mixture = c("P"))
res_mois_terrier2 <- pcount(~mois+terlapin ~1, count, mixture = c("P"))
res_mois_terrier3 <- pcount(~mois+terlapin ~terlapin, count, mixture = c("P"))
res_mois_habitat <- pcount(~mois ~habitat, count, mixture = c("P"))
res_mois_habitat2 <- pcount(~mois+habitat ~1, count, mixture = c("P"))
res_mois_habitat3 <- pcount(~mois+habitat ~habitat, count, mixture = c("P"))

newData2 <- data.frame(terlapin = seq(0, 20, by = 1))
E.p <- predict(res_mois_terrier, type = "state", newdata = newData2, appendData = TRUE)
plot(Predicted ~ terlapin, E.p, type = "l", ylim = c(0, 10),
     xlab = "nb terriers de lapin",
     ylab = "Expected abundance")
lines(lower ~ terlapin, E.p, type = "l", col = gray(0.5))
lines(upper ~ terlapin, E.p, type = "l", col = gray(0.5))
