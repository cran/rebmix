\name{coef.REBMIX}
\alias{coef.REBMIX}
\title{
Prints Univariate or Multivariate REBMIX Coefficients
}
\description{
Returns the \code{w} and \code{Theta} output at \code{pos} for class REBMIX.
}
\usage{
\method{coef}{REBMIX}(object, pos = 1, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} to be printed. The default value is \code{1}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\examples{
## Generate and print the simulated 7 dataset.

Theta <- rbind(pdf = rep("Weibull", 3),
  theta1 = c(2.0, 10.0, 30),
  theta2 = c(2.0, 4.0, 7.0))

simulated7 <- RNGMIX(Dataset = "simulated7",
  rseed = -1,
  n = c(40, 60, 50),
  Theta = Theta)
  
simulated7

## Estimate number of components, component weights and component parameters.

simulated7est <- REBMIX(Dataset = simulated7$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 6, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = 8:25, 
  b = 1.0)
  
## Print coefficients and plot the finite mixture.  

coef(simulated7est)

plot(simulated7est)
}
\keyword{print}

