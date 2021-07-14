##############################################
## R sources for reproducing the results in ##
##   Marko Nagode, Branislav Panic,         ##
##   Jernej Klemenc, Simon Oman:            ##
##   Fault Detection and Classification     ##
##   with the rebmix R Package              ##
##############################################

options(prompt = "> ", continue = "+ ", width = 70,
  useFancyQuotes = FALSE, digits = 3)

library("rebmix")
library("e1071")
library("FNN")
library("MASS")

data("bearings")
data("steel.plates")
data("sensorless.drive")

bearings <- bearings[, c(1, 5, 6, 9, 10, 11, 14)]
steel.plates <- steel.plates[, c(2, 9, 17, 19, 24, 28)]

# Data normalization.

range <- apply(bearings, 2, range)

for (i in 1:(ncol(range) - 1)) {
  bearings[, i] <- (bearings[, i] - range[1, i]) / (range[2, i] - range[1, i])
}

range <- apply(steel.plates, 2, range)

for (i in 1:(ncol(range) - 1)) {
  steel.plates[, i] <- (steel.plates[, i] - range[1, i]) / (range[2, i] - range[1, i])
}

range <- apply(sensorless.drive, 2, range)

for (i in 1:(ncol(range) - 1)) {
  sensorless.drive[, i] <- (sensorless.drive[, i] - range[1, i]) / (range[2, i] - range[1, i])
}

# Data split into train and test datasets.

set.seed(1)

Bearings <- split(p = 0.6, Dataset = bearings, class = 7)
Steel.plates <- split(p = 0.6, Dataset = steel.plates, class = 6)
Sensorless.drive <- split(p = 0.6, Dataset = sensorless.drive, class = 4)

########## rebmix ##########

# Bearings dataset.

EM <- new("EM.Control", strategy = "exhaustive", variant = "ECM")

system.time({

bearings.model <- REBMIX(model = "REBMVNORM", 
  Dataset = Bearings@train, 
  Preprocessing  = "k-nearest neighbour",
  Criterion = "AIC",
  EMcontrol = EM)

bearings.class <- RCLSMIX(model = "RCLSMVNORM", 
  x = list(bearings.model), 
  Dataset = Bearings@test, 
  Zt = Bearings@Zt)

})

bearings.class

plot(bearings.class, nrow = 5, ncol = 3)

# Steel.plates dataset.

a.strategy(EM) <- "best"
a.variant(EM) <- "EM"

system.time({

steel.plates.model <- REBMIX(model = "REBMVNORM", 
  Dataset = Steel.plates@train, 
  Preprocessing  = "histogram",
  Criterion = "BIC",
  K = 2:100,
  EMcontrol = EM)

steel.plates.class <- RCLSMIX(model = "RCLSMVNORM", 
  x = list(steel.plates.model), 
  Dataset = Steel.plates@test, 
  Zt = Steel.plates@Zt)

})

steel.plates.class

plot(steel.plates.class, nrow = 2, ncol = 5)

# Sensorless.drive dataset.

a.strategy(EM) <- "single"

system.time({

bins <- optbins(Dataset = Sensorless.drive@train, Rule = "Knuth equal", kmin = 2, kmax = 100)

sensorless.drive.model <- REBMIX(model = "REBMIX", 
  Dataset = Sensorless.drive@train, 
  Preprocessing  = "histogram",
  pdf = rep("normal", 3), 
  K = bins, 
  EMcontrol = EM)

sensorless.drive.class <- RCLSMIX(model = "RCLSMIX", 
  x = list(sensorless.drive.model), 
  Dataset = Sensorless.drive@test, 
  Zt = Sensorless.drive@Zt)
})

sensorless.drive.class

plot(sensorless.drive.class, nrow = 3, ncol = 1)

# Bearings data preparation.

train <- NULL; Zr <- NULL

for (i in 1:length(Bearings@ntrain)) {
  Zr <- c(Zr, Bearings@Zr[[i]])
  train <- rbind(train, Bearings@train[[i]])
}

rownames(train) <- NULL

Zr <- as.factor(Zr)

test <- Bearings@test; Zt <- as.numeric(Bearings@Zt)

rownames(test) <- NULL

########## svm ##########

system.time({

model <- svm(x = train, y = Zr)

Zp <- predict(model, test)

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error

########## knn ##########

system.time({

knn <- knn(train = train, test = test, cl = Zr, k = 10)

Zp <- knn[1:nrow(test)]

})

Zp <- as.numeric(as.vector(Zp))

Error <- 1.0 - sum(Zt == Zp) / length(Zp)

Error

########## lda ##########

system.time({

lda <- lda(train, Zr)

Zp <- predict(lda, test)$class

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error

# Steel.plates data preparation.

train <- NULL; Zr <- NULL

for (i in 1:length(Steel.plates@ntrain)) {
  Zr <- c(Zr, Steel.plates@Zr[[i]])
  train <- rbind(train, Steel.plates@train[[i]])
}

rownames(train) <- NULL

Zr <- as.factor(Zr)

test <- Steel.plates@test; Zt <- as.numeric(Steel.plates@Zt)

rownames(test) <- NULL

########## svm ##########

system.time({

model <- svm(x = train, y = Zr)

Zp <- predict(model, test)

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error

########## knn ##########

system.time({

knn <- knn(train = train, test = test, cl = Zr, k = 10)

Zp <- knn[1:nrow(test)]

})

Zp <- as.numeric(as.vector(Zp))

Error <- 1.0 - sum(Zt == Zp) / length(Zp)

Error

########## lda ##########

system.time({

lda <- lda(train, Zr)

Zp <- predict(lda, test)$class

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error

# Sensorless.drive data preparation.

train <- NULL; Zr <- NULL

for (i in 1:length(Sensorless.drive@ntrain)) {
  Zr <- c(Zr, Sensorless.drive@Zr[[i]])
  train <- rbind(train, Sensorless.drive@train[[i]])
}

rownames(train) <- NULL

Zr <- as.factor(Zr)

test <- Sensorless.drive@test; Zt <- as.numeric(Sensorless.drive@Zt)

rownames(test) <- NULL

########## svm ##########

system.time({

model <- svm(x = train, y = Zr)

Zp <- predict(model, test)

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error

########## knn ##########

system.time({

knn <- knn(train = train, test = test, cl = Zr, k = 10)

Zp <- knn[1:nrow(test)]

})

Zp <- as.numeric(as.vector(Zp))

Error <- 1.0 - sum(Zt == Zp) / length(Zp)

Error

########## lda ##########

system.time({

lda <- lda(train, Zr)

Zp <- predict(lda, test)$class

})

Zp <- as.numeric(Zp)

Error <- 1.0 - sum(Zt == Zp) / length(Zt)

Error
