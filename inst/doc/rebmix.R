### R code from vignette source 'rebmix.Rnw'
### Encoding: CP1250

###################################################
### code chunk number 1: rebmix-code
###################################################
##############################################
## R sources for reproducing the results in ##
##   Marko Nagode:                          ##
##   rebmix: The Rebmix Package             ##
##############################################

options(prompt = "R> ", continue = "+  ", width = 80,
  useFancyQuotes = FALSE, digits = 3)


###################################################
### code chunk number 2: rebmix-code
###################################################
###################
## Preliminaries ##
###################

## load package and set prompt before starting new page to TRUE.

library("rebmix")
devAskNewPage(ask = TRUE)


###################################################
### code chunk number 3: rebmix-code
###################################################
######################
##  Gamma datasets  ##
######################

## Generate gamma datasets.

n <- c(100, 100, 100, 100)

Theta <- rbind(pdf = "gamma",
  theta1 = c(1/100, 1/100, 1/100, 1/100),
  theta2 = c(200, 400, 600, 800))

gamma1 <- RNGMIX(Dataset = "gamma1", n = n, Theta = Theta)

n <- c(40, 360)

Theta <- rbind(pdf = "gamma",
  theta1 = c(1/27, 1/270),
  theta2 = c(9, 90))

gamma2 <- RNGMIX(Dataset = "gamma2", n = n, Theta = Theta)

n <- c(80, 240, 80)

Theta <- rbind(pdf = "gamma",
  theta1 = c(1/20, 1, 1/20),
  theta2 = c(40, 6, 200))

gamma3 <- RNGMIX(Dataset = "gamma3 ", n = n, Theta = Theta)


###################################################
### code chunk number 4: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

gamma1est <- REBMIX(Dataset = gamma1$Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = c("AIC", "BIC"),
  Variables = "continuous",
  pdf = "gamma",
  K = 30:80)

gamma2est <- REBMIX(Dataset = gamma2$Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  Variables = "continuous",
  pdf = "gamma",
  K = 30:80)

gamma3est <- REBMIX(Dataset = gamma3$Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  Variables = "continuous",
  pdf = "gamma",
  K = 30:80)


###################################################
### code chunk number 5: rebmix-code
###################################################
summary(gamma1est)


###################################################
### code chunk number 6: gamma2-fig
###################################################
plot(gamma2est, pos = 1, what = c("den", "dis"), ncol = 2, npts = 1000)


###################################################
### code chunk number 7: rebmix-code
###################################################
coef(gamma2est)


###################################################
### code chunk number 8: rebmix-code
###################################################
## Bootstrap finite mixture.
gamma3boot <- boot.REBMIX(x = gamma3est, pos = 1, Bootstrap = "p", B = 10, n = NULL, replace = TRUE, prob = NULL)


###################################################
### code chunk number 9: rebmix-code
###################################################
gamma3boot


###################################################
### code chunk number 10: rebmix-code
###################################################
summary(gamma3boot)


###################################################
### code chunk number 11: rebmix-code
###################################################
#########################
##   Poisson dataset   ##
#########################

## Generate the Poisson dataset.

n <- c(200, 200, 200)

Theta <- rbind(rep("Poisson", 3), c(3, 9, 15), rep("Poisson", 3), c(2, 10, 16))

poisson <- RNGMIX(Dataset = paste("Poisson_", 1:100, sep = ""), n = n, Theta = Theta)


###################################################
### code chunk number 12: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

poissonest <- REBMIX(Dataset = poisson$Dataset,
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "MDL5",
  Variables = rep("discrete", 2),
  pdf = rep("Poisson", 2),
  K = 1)

c <- as.numeric(poissonest$summary$c)
IC <- as.numeric(poissonest$summary$IC)


###################################################
### code chunk number 13: rebmix-code
###################################################
## Visualize results.

summary(c)
summary(IC, digits = 5)


###################################################
### code chunk number 14: poisson-fig
###################################################
plot(poissonest, pos = 58, what = c("dens", "marg", "IC", "D", "logL"), nrow = 2, ncol = 3, npts = 1000)


###################################################
### code chunk number 15: rebmix-code
###################################################
rm(list = ls())


