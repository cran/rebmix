### R code from vignette source 'rebmix.Rnw'

###################################################
### code chunk number 1: rebmix-code
###################################################
##############################################
## R sources for reproducing the results in ##
##   Marko Nagode:                          ##
##   rebmix: An R Package for Continuous    ##
##   and Discrete Finite Mixture Models     ##
##############################################

options(prompt = "R> ", continue = "+  ", width = 70,
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
####################
## Galaxy dataset ##
####################

## Load galaxy dataset.

data("galaxy")

galaxyest <- list(normal = NULL, lognormal = NULL, Weibull = NULL)

## Estimate number of components, component weights and component parameters.

pdf <- c("normal", "lognormal", "Weibull")

for (i in 1:3) {
  galaxyest[[i]] <- REBMIX(Dataset = list(galaxy = galaxy),
    Preprocessing = c("histogram", "Parzen window"),
    cmax = 8,
    Criterion = c("AIC", "BIC"),
    Variables = "continuous",
    pdf = pdf[i],
    K = 7:20)
}


###################################################
### code chunk number 4: rebmix-code
###################################################
summary(galaxyest$normal)

summary(galaxyest$lognormal)

summary(galaxyest$Weibull)


###################################################
### code chunk number 5: galaxy-fig
###################################################
## Visualize Figure 1.

plot(galaxyest$Weibull, pos = 1, what = c("den", "dis"), ncol = 2, npts = 1000)


###################################################
### code chunk number 6: rebmix-code
###################################################
plot(galaxyest$Weibull, pos = 1, what = c("den", "dis"), ncol = 2, npts = 1000)


###################################################
### code chunk number 7: rebmix-code
###################################################
coef(galaxyest$Weibull, pos = 1)


###################################################
### code chunk number 8: rebmix-code
###################################################
##################
## Iris dataset ##
##################

data("iris")

## Remove the Species column from the dataset.

iris <- iris[, !(colnames(iris) %in% "Species")]


###################################################
### code chunk number 9: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

irisest <- REBMIX(Dataset = list(iris = iris),
  Preprocessing = c("histogram", "Parzen window", "k-nearest neighbour"),
  Criterion = c("AIC", "AWE", "BIC", "CLC", "ICL", "ICL-BIC"),
  Variables = rep("continuous", 4),
  pdf = rep("normal", 4),
  K = list(6:25, 6:25, 3:13))


###################################################
### code chunk number 10: rebmix-code
###################################################
summary(irisest)


###################################################
### code chunk number 11: rebmix-code
###################################################
plot(irisest, pos = 5, what = c("den", "IC", "logL", "D"), nrow = 3, ncol = 3, npts = 1000)


###################################################
### code chunk number 12: iris-fig
###################################################
## Visualize Figure 2.

plot(irisest, pos = 5, what = c("den", "IC", "logL", "D"), nrow = 3, ncol = 3, npts = 1000)


###################################################
### code chunk number 13: rebmix-code
###################################################
##################
## Wine dataset ##
##################

data("wine")

## Remove the Cultivar column from the dataset.

wine <- wine[, !(colnames(wine) %in% "Cultivar")]

## Estimate number of components, component weights and component parameters.

wineest <- REBMIX(Dataset = list(wine = wine),
  Preprocessing = c("histogram", "Parzen window"),
  Criterion = c("AIC", "AWE", "BIC", "CLC", "ICL", "ICL-BIC"),
  Variables = rep("continuous", 13),
  pdf = rep("normal", 13),
  K = 8:27)


###################################################
### code chunk number 14: rebmix-code
###################################################
summary(wineest)


###################################################
### code chunk number 15: rebmix-code
###################################################
#######################
## Complex 1 dataset ##
#######################

## Generate the complex 1 dataset.

n <- c(998, 263, 1086, 487, 213, 1076, 232,
  784, 840, 461, 773, 24, 811, 1091, 861)

Theta <- rbind(pdf = "normal",
  theta1 = c(688.4, 265.1, 30.8, 934, 561.6, 854.9, 883.7,
  758.3, 189.3, 919.3, 98, 143, 202.5, 628, 977),
  theta2 = c(12.4, 14.6, 14.8, 8.4, 11.7, 9.2, 6.3, 10.2,
  9.5, 8.1, 14.7, 11.7, 7.4, 10.1, 14.6))

complex1 <- RNGMIX(Dataset = "complex1",
  n = n,
  Theta = Theta)


###################################################
### code chunk number 16: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

time <- system.time(
complex1est1 <- REBMIX(Dataset = complex1$Dataset,
  Preprocessing = "histogram",
  cmax = 20,
  Criterion = "BIC",
  Variables = "continuous",
  pdf = "normal",
  K = seq(14, 200, 4))
)


###################################################
### code chunk number 17: rebmix-code
###################################################
plot(complex1est1, npts = 1000)


###################################################
### code chunk number 18: complex1-fig
###################################################
## Visualize Figure 3.

plot(complex1est1, npts = 1000)


###################################################
### code chunk number 19: rebmix-code
###################################################
#######################
## Complex 2 dataset ##
#######################

## Generate the complex 2 dataset.

n <- c(390, 110, 300, 70, 130)

Theta <- rbind(pdf1 = rep("lognormal", 5),
  theta1.1 = c(0.8, 1.3, 3.4, 2.7, 4.3),
  theta2.1 = c(0.5, 0.7, 0.2, 0.4, 0.1),
  pdf2 = rep("Poisson", 5),
  theta1.2 = c(10.0, 7.3, 1.7, 3.3, 5.0),
  pdf3 = rep("binomial", 5),
  theta1.3 = c(10, 10, 10, 10, 10),
  theta2.3 = c(0.9, 0.7, 0.5, 0.3, 0.1),
  pdf4 = rep("Weibull", 5),
  theta1.4 = c(20, 45, 60, 90, 120),
  theta2.4 = c(2.0, 3.1, 6.3, 2.5, 7.0))

complex2 <- RNGMIX(Dataset = "complex2",
  n = n,
  Theta = Theta)


###################################################
### code chunk number 20: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.
time <- system.time(
complex2est <- REBMIX(Dataset = complex2$Dataset,
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  Variables = c("continuous", "discrete", "discrete", "continuous"),
  pdf = c("lognormal", "Poisson", "binomial", "Weibull"),
  Theta1 = c(NA, NA, 10, NA),
  K = seq(10, 64, 1))
)


###################################################
### code chunk number 21: rebmix-code
###################################################
plot(complex2est, what = c("dens", "marg", "IC", "D"), nrow = 4, ncol = 3, npts = 1000)


###################################################
### code chunk number 22: complex2-fig
###################################################
## Visualize Figure 4.

plot(complex2est, what = c("dens", "marg", "IC", "D"), nrow = 4, ncol = 3, npts = 1000)


###################################################
### code chunk number 23: rebmix-code
###################################################
## Bootstrap finite mixture.

complex2boot <- boot.REBMIX(x = complex2est, pos = 1, Bootstrap = "p", B = 100, n = NULL, replace = TRUE, prob = NULL)


###################################################
### code chunk number 24: rebmix-code
###################################################
complex2boot


###################################################
### code chunk number 25: rebmix-code
###################################################
#########################
## Simulated 1 dataset ##
#########################

## Generate the simulated 1 dataset.

n <- c(75, 100, 125, 150, 175)

Theta <- rbind(rep("normal", 5),
  c(10, 8.5, 12, 13, 7),
  c(1, 1, 1, 2, 3),
  rep("normal", 5),
  c(12, 10.5, 14, 15, 9),
  c(1, 1, 1, 2, 3),
  rep("normal", 5),
  c(10, 8.5, 12, 7, 13),
  c(1, 1, 1, 2, 3),
  rep("normal", 5),
  c(12, 10.5, 14, 9, 15),
  c(1, 1, 1, 2, 3))

simulated1 <- RNGMIX(Dataset = paste("Simulated1_", 1:100, sep = ""),
  n = n,
  Theta = Theta)


###################################################
### code chunk number 26: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

time <- system.time(
simulated1est1 <- REBMIX(simulated1$Dataset,
  Preprocessing = "histogram",
  Criterion = "BIC",
  Variables = rep("continuous", 4),
  pdf = rep("normal", 4),
  K = seq(10, 28, 2))
)

c <- as.numeric(simulated1est1$summary$c)
IC <- as.numeric(simulated1est1$summary$IC)


###################################################
### code chunk number 27: rebmix-code
###################################################
## Visualize results.

summary(c)
summary(IC, digits = 5)


###################################################
### code chunk number 28: rebmix-code
###################################################
rm(list = ls())


