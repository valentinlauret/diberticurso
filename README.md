# (d)IBERticurso

**What:** This is a workshop on introduction to Bayesian inference. Hopefully, you will learn about how to make Bayesian models for population ecology. Our hope is to provide you with what you need to go your own path in Bayesian modelling. We will present how to run occupancy model, N-mixture model, and Spatial Capture-Recapture models.

**For whom:** This is a workshop for ecologists. No previous experience with Nimble or Bayesian statistics is assumed, and few knowledge of R is required.

**How:** Through a combination of lectures and live demonstrations, you will get acquainted with the Bayesian approach, MCMC methods and the R Nimble package to fit single-site, multi-site, multi-state and multi-event models to capture-recapture data within the hidden Markov modeling (HMM) framework.

**Who:** Valentin Lauret, Javi Fernandez-Lopez

**When:** May 13-15, 2024

**Where:** IREC, Royal City.

## Requirements

  * Install `R` and `RStudio`.

  * Install `Nimble` following these guidelines. Then run the following code in R. If that runs without error, you’re all set. If not, please get in touch with us.

```
library(nimble)
  code <- nimbleCode({
  y ~ dnorm(0,1)
  })
  model <- nimbleModel(code)
  cModel <- compileNimble(model)
  
```

  * Install the following R packages: `tidyverse`, `mcmcplots`, `coda`. You can install them all at once by running the following code in R:

```
install.packages(c("tidyverse", "mcmcplots", "coda"))
```


