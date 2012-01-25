\name{print.RNGMIX}
\alias{print.RNGMIX}
\title{
Prints Univariate or Multivariate RNGMIX Output
}
\description{
Returns the \code{\link[base]{print}} output for class RNGMIX.
}
\usage{
\method{print}{RNGMIX}(x, ...)
}
\arguments{
\item{x}{
an object of class \code{RNGMIX}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\value{
\code{print.RNGMIX} returns (invisibly) the full value of x itself.
}
\examples{
## Generate and print the simulated 5 dataset.

Theta <- rbind(pdf = rep("Poisson", 2),
  theta1 = c(1.0, 10.0))

simulated5 <- RNGMIX(Dataset = paste("simulated5_", 1:5, sep = ""),
  rseed = -1,
  n = c(100, 50),
  Theta = Theta)
  
print(simulated5, pos = 3)

## Estimate number of components, component weights and component parameters.

simulated5est <- REBMIX(Dataset = simulated5$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 4, 
  Criterion = "AIC", 
  Variables = "discrete",
  pdf = "Poisson", 
  K = 1, 
  b = 1.0)
  
## Print and plot the finite mixture.    

print(simulated5est, pos = 2)

plot(simulated5est, pos = 3)
}
\keyword{print}

