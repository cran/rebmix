### R code from vignette source 'rebmix.Rnw'

###################################################
### code chunk number 1: rebmix-code
###################################################
##############################################
## R sources for reproducing the results in ##
##              rebmix package              ##
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

library(rebmix)
devAskNewPage(ask = TRUE)


###################################################
### code chunk number 3: rebmix-code
###################################################
######################
##  Gamma datasets  ##
######################

## Generate gamma datasets.

n <- c(100, 100, 100, 100)

Theta <- new("RNGMIX.Theta", c = 4, pdf = "gamma")

a.theta1(Theta) <- rep(1/100, 4)
a.theta2(Theta) <- c(200, 400, 600, 800)

gamma1 <- RNGMIX(Dataset.name = "gamma1", n = n, Theta = a.Theta(Theta))

n <- c(40, 360)

Theta <- new("RNGMIX.Theta", c = 2, pdf = "gamma")

a.theta1(Theta) <- c(1/27, 1 / 270)
a.theta2(Theta) <- c(9, 90)

gamma2 <- RNGMIX(Dataset.name = "gamma2", n = n, Theta = a.Theta(Theta))

n <- c(80, 240, 80)

Theta <- new("RNGMIX.Theta", c = 3, pdf = "gamma")

a.theta1(Theta) <- c(1/20, 1, 1/20)
a.theta2(Theta) <- c(40, 6, 200)

gamma3 <- RNGMIX(Dataset.name = "gamma3", rseed = -4, n = n, Theta = a.Theta(Theta))


###################################################
### code chunk number 4: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

gamma1est <- REBMIX(Dataset = a.Dataset(gamma1),
  Preprocessing = "kernel density estimation",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma")

gamma2est <- REBMIX(Dataset = a.Dataset(gamma2),
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma")

gamma3est <- REBMIX(Dataset = a.Dataset(gamma3),
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  K = 23:27)


###################################################
### code chunk number 5: gamma3-fig
###################################################
plot(gamma3est, pos = 1, what = c("pdf", "marginal cdf"), ncol = 2, npts = 1000)


###################################################
### code chunk number 6: rebmix-code
###################################################
summary(gamma2est)

a.theta1.all(gamma1est, pos = 1)

a.theta2.all(gamma1est, pos = 1)


###################################################
### code chunk number 7: rebmix-code
###################################################
## Bootstrap finite mixture.

gamma3boot <- boot(x = gamma3est, pos = 1, Bootstrap = "p", B = 10)

gamma3boot

summary(gamma3boot)


###################################################
### code chunk number 8: rebmix-code
###################################################
## EM.control object creation.

EM <- new("EM.Control",
  strategy = "best",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

gamma1est.em <- REBMIX(Dataset = a.Dataset(gamma1),
  Preprocessing = "kernel density estimation",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  EMcontrol = EM)

gamma2est.em <- REBMIX(Dataset = a.Dataset(gamma2),
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  EMcontrol = EM)

gamma3est.em <- REBMIX(Dataset = a.Dataset(gamma3),
  Preprocessing = "histogram",
  cmax = 8,
  Criterion = "BIC",
  pdf = "gamma",
  K = 23:27,
  EMcontrol = EM)

summary(gamma1est.em)

summary(gamma2est.em)

summary(gamma3est.em)


###################################################
### code chunk number 9: rebmix-code
###################################################
#########################
##   Poisson dataset   ##
#########################

## Generate the Poisson dataset.

n <- c(200, 200, 200)

Theta <- new("RNGMIX.Theta", c = 3, pdf = rep("Poisson", 2))

a.theta1(Theta, 1) <- c(3, 2)
a.theta1(Theta, 2) <- c(9, 10)
a.theta1(Theta, 3) <- c(15, 16)

poisson <- RNGMIX(Dataset.name = paste("Poisson_", 1:10, sep = ""), n = n, Theta = a.Theta(Theta))


###################################################
### code chunk number 10: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

poissonest <- REBMIX(Dataset = a.Dataset(poisson),
  Preprocessing = "histogram",
  cmax = 10,
  Criterion = "MDL5",
  pdf = rep("Poisson", 2),
  K = 1)


###################################################
### code chunk number 11: poisson-fig
###################################################
plot(poissonest, pos = 1, what = c("pdf", "marginal pdf", "IC", "D", "logL"), nrow = 2, ncol = 3, npts = 1000)


###################################################
### code chunk number 12: poisson-clu-fig
###################################################
poissonclu <- RCLRMIX(x = poissonest, pos = 1, Zt = a.Zt(poisson))

plot(poissonclu)


###################################################
### code chunk number 13: rebmix-code
###################################################
## Visualize results.

summary(poissonest)

a.theta1.all(poissonest, pos = 1)

a.theta2.all(poissonest, pos = 1)


###################################################
### code chunk number 14: rebmix-code
###################################################
## EM.control object creation.

EM <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

poissonest.em <- REBMIX(Dataset = a.Dataset(poisson),
  Preprocessing = "histogram",
  cmax = 10,
  Criterion = "MDL5",
  pdf = rep("Poisson", 2),
  K = 1,
  EMcontrol  = EM)

summary(poissonest.em)


###################################################
### code chunk number 15: rebmix-code
###################################################
###################################
##  Multivariate normal dataset  ##
###################################

## Generate normal dataset.

n <- c(50, 50, 50, 50, 50)

Theta <- new("RNGMVNORM.Theta", c = 5, d = 2)

a.theta1(Theta, 1) <- c(2.7, 3.7)

a.theta1(Theta, 2) <- c(5.7, 9.1)

a.theta1(Theta, 3) <- c(2.0, 9.0)

a.theta1(Theta, 4) <- c(9.5, 6.6)

a.theta1(Theta, 5) <- c(6.3, 0.6)

a.theta2(Theta, 1) <- c(0.9, -0.1, -0.1, 0.4)

a.theta2(Theta, 2) <- c(2.8, -1.3, -1.3, 1.5)

a.theta2(Theta, 3) <- c(0.1, 0.0, 0.0, 0.3)

a.theta2(Theta, 4) <- c(1.3, -0.4, -0.4, 0.4)

a.theta2(Theta, 5) <- c(0.5, 0.3, 0.3, 2.5)

mvnorm.simulated <- RNGMIX(model = "RNGMVNORM",
  Dataset.name = "mvnormdataset",
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))


###################################################
### code chunk number 16: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters.

mvnormest <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 20,
  Criterion = "BIC")


###################################################
### code chunk number 17: mvnorm-fig
###################################################
plot(mvnormest)


###################################################
### code chunk number 18: mvnorm-clu-fig
###################################################
mvnormclu <- RCLRMIX(model = "RCLRMVNORM", x = mvnormest)

plot(mvnormclu)


###################################################
### code chunk number 19: rebmix-code
###################################################
summary(mvnormest)


###################################################
### code chunk number 20: rebmix-code
###################################################
summary(mvnormclu)


###################################################
### code chunk number 21: rebmix-code
###################################################
## EM.control object creation.

EM <- new("EM.Control",
  strategy = "single",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Optimal K estimation.

K <- optbins(Dataset = a.Dataset(mvnorm.simulated), Rule = "Knuth equal", kmin = 2, kmax = 100)

## Finite mixture estimation

mvnormest.em <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 20,
  K = K,
  Criterion = "BIC",
  EMcontrol = EM)

summary(mvnormest.em)


###################################################
### code chunk number 22: rebmix-code
###################################################
## EM.control object creation.

CEM <- new("EM.Control",
  strategy = "exhaustive",
  variant = "ECM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Estimate number of components, component weights and component parameters.

mvnormest.cem <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 10,
  Criterion = "ICL",
  EMcontrol = CEM)

mvnorm.clu <- RCLRMIX(model = "RCLRMVNORM", x = mvnormest.cem)


###################################################
### code chunk number 23: mvnorm-clu-fig-cem
###################################################
plot(mvnorm.clu)


###################################################
### code chunk number 24: rebmix-code
###################################################
## Create the EM.Control object to utilize one of the REBMIX&EM strategies.

EM.normal <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Estimate number of components, component weights and component parameters.

mvnormestest.em.normal <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC",
  EMcontrol = EM.normal)

cat("Total number of EM algorithm iterations: ",
  a.summary.EM(mvnormestest.em.normal, pos = 1, col.name = "total.iterations.nbr"),
  ". Value of BIC: ", a.summary(mvnormestest.em.normal, pos = 1, col.name = "IC"))


###################################################
### code chunk number 25: rebmix-code
###################################################
## Create the EM.Control object to utilize one of the REBMIX&EM strategies.

EM.fixed1.5 <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.5,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Estimate number of components, component weights and component parameters.

mvnormest.em.fixed1.5 <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC",
  EMcontrol = EM.fixed1.5)

cat("Total number of EM algorithm iterations: ",
  a.summary.EM(mvnormest.em.fixed1.5, pos = 1, col.name = "total.iterations.nbr"),
  ". Value of BIC: ", a.summary(mvnormest.em.fixed1.5, pos = 1, col.name = "IC"))


###################################################
### code chunk number 26: rebmix-code
###################################################
## Create the EM.Control object to utilize one of the REBMIX&EM strategies.

EM.line <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "line",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Estimate number of components, component weights and component parameters.

mvnormest.em.line <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC",
  EMcontrol = EM.line)

cat("Total number of EM algorithm iterations: ",
  a.summary.EM(mvnormest.em.line, pos = 1, col.name = "total.iterations.nbr"),
  ". Value of BIC: ", a.summary(mvnormest.em.line, pos = 1, col.name = "IC"))


###################################################
### code chunk number 27: rebmix-code
###################################################
## Create the EM.Control object to utilize one of the REBMIX&EM strategies.

EM.golden <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "golden",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 0)

## Estimate number of components, component weights and component parameters.

mvnormest.em.golden <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(mvnorm.simulated),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC",
  EMcontrol = EM.golden)

cat("Total number of EM algorithm iterations: ",
  a.summary.EM(mvnormest.em.golden, pos = 1, col.name = "total.iterations.nbr"),
  ". Value of BIC: ", a.summary(mvnormest.em.golden, pos = 1, col.name = "IC"))


###################################################
### code chunk number 28: rebmix-code
###################################################
data(sensorlessdrive)

## Split dataset into train (75%) and test (25%) subsets.

set.seed(5)

Drive <- split(p = 0.75, Dataset = sensorlessdrive, class = 4)


###################################################
### code chunk number 29: rebmix-code
###################################################
## Estimate number of components, component weights and component
## parameters for train subsets.

driveest <- REBMIX(model = "REBMVNORM",
  Dataset = a.train(Drive),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC")


###################################################
### code chunk number 30: rebmix-code
###################################################
## Selected features.

drivecla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(driveest),
  Dataset = a.test(Drive),
  Zt = a.Zt(Drive))


###################################################
### code chunk number 31: rebmix-code
###################################################
drivecla

summary(drivecla)


###################################################
### code chunk number 32: drive-cla-fig
###################################################
# Plot selected features.

plot(drivecla, nrow = 3, ncol = 2)


###################################################
### code chunk number 33: rebmix-code
###################################################
## EM.control object creation.

EM <- new("EM.Control",
  strategy = "exhaustive",
  variant = "EM",
  acceleration = "fixed",
  acceleration.multiplier = 1.0,
  tolerance = 1e-4,
  maximum.iterations = 1000,
  K = 300)

## Estimate number of components, component weights and component
## parameters for train subsets.

driveest <- REBMIX(model = "REBMVNORM",
  Dataset = a.train(Drive),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC",
  EMcontrol = EM)

drivecla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(driveest),
  Dataset = a.test(Drive),
  Zt = a.Zt(Drive))

summary(drivecla)


###################################################
### code chunk number 34: rebmix-code
###################################################
data(adult)

## Find complete cases.

adult <- adult[complete.cases(adult),]

## Replace levels with numbers.

adult <- as.data.frame(data.matrix(adult))


###################################################
### code chunk number 35: rebmix-code
###################################################
## Find numbers of levels.

cmax <- unlist(lapply(apply(adult[, c(-1, -16)], 2, unique), length))

cmax


###################################################
### code chunk number 36: rebmix-code
###################################################
## Split adult dataset into train and test subsets for two Incomes
## and remove Type and Income columns.

Adult <- split(p = list(type = 1, train = 2, test = 1),
  Dataset = adult, class = 16)


###################################################
### code chunk number 37: rebmix-code
###################################################
## Estimate number of components, component weights and component parameters
## for the set of chunks 1:14.

adultest <- list()

for (i in 1:14) {
  adultest[[i]] <- REBMIX(Dataset = a.train(chunk(Adult, i)),
    Preprocessing = "histogram",
    cmax = min(120, cmax[i]),
    Criterion = "BIC",
    pdf = "Dirac",
    K = 1)
}


###################################################
### code chunk number 38: rebmix-code
###################################################
## Class membership prediction based upon the best first search algorithm.

adultcla <- BFSMIX(x = adultest,
  Dataset = a.test(Adult),
  Zt = a.Zt(Adult))


###################################################
### code chunk number 39: rebmix-code
###################################################
adultcla

summary(adultcla)


###################################################
### code chunk number 40: adult-cla-fig
###################################################
## Plot selected chunks.

plot(adultcla, nrow = 5, ncol = 2)


###################################################
### code chunk number 41: rebmix-code
###################################################
rm(list = ls())


