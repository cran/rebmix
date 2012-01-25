\name{summary.REBMIX}
\alias{summary.REBMIX}
\title{
Prints Univariate or Multivariate REBMIX Summary
}
\description{
Returns the \code{summary} output for class REBMIX.
}
\usage{
\method{summary}{REBMIX}(object, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{\dots}{
further arguments to \code{\link[base]{summary}}.
}
}
\examples{
## Generate and print the simulated 6 dataset.

Theta <- rbind(pdf = rep("Dirac", 4),
  theta1 = c(1, 3, 4, 8))

simulated6 <- RNGMIX(Dataset = "simulated6",
  rseed = -1,
  n = c(40, 20, 50, 10),
  Theta = Theta)
  
simulated6

## Estimate number of components, component weights and component parameters.

simulated6est <- REBMIX(Dataset = simulated6$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 8, 
  Criterion = "D", 
  Variables = "discrete",
  pdf = "Dirac", 
  K = 1, 
  b = 0.0)
  
## Print summary and plot the finite mixture.  

summary(simulated6est)

plot(simulated6est)
}
\keyword{print}

