\name{RCLRMIX-methods}
\docType{methods}
\alias{RCLRMIX}
\alias{RCLRMIX-methods}
\alias{RCLRMIX,RCLRMIX-method}
\alias{RCLRMIX,RCLRMVNORM-method}
\alias{show,RCLRMIX-method}
\alias{show,RCLRMVNORM-method}
\alias{summary,RCLRMIX-method}
\alias{summary,RCLRMVNORM-method}
\title{
Predicts Cluster Membership Based Upon a Model Trained by REBMIX
}
\description{
Returns as default the RCLRMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities. If \code{model} equals \code{"RCLRMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{RCLRMIX}{RCLRMIX}(model = "RCLRMIX", x = NULL, pos = 1, Zt = factor(), ...)
## ... and for other signatures
\S4method{summary}{RCLRMIX}(object, ...)
## ... and for other signatures     
}
\arguments{
\item{model}{
see Methods section below.
}
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{x@summary} to be clustered. The default value is \code{1}.
}
\item{Zt}{
a factor of true cluster membership. The default value is \code{factor()}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used; additional arguments affecting the summary produced.
}
}
\value{
Returns an object of class \code{RCLRMIX} or \code{RCLRMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLRMIX")}}{a character giving the default class name \code{"RCLRMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.}
\item{\code{signature(model = "RCLRMVNORM")}}{a character giving the class name \code{"RCLRMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "RCLRMIX")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(object = "RCLRMVNORM")}}{an object of class \code{RCLRMVNORM}.}
}
}
\author{Marko Nagode}
\examples{
devAskNewPage(ask = TRUE)

# Generate Poisson dataset.

n <- c(500, 200, 400)

Theta <- list(pdf1 = "Poisson",
  theta1.1 = 3,
  theta2.1 = NA,
  pdf2 = "Poisson",
  theta1.2 = 12,
  theta2.2 = NA,
  pdf3 = "Poisson",
  theta1.3 = 36,
  theta2.3 = NA)

poisson <- RNGMIX(Dataset.name = "Poisson_1", n = n, Theta = Theta)

# Estimate number of components, component weights and component parameters.

poissonest <- REBMIX(Dataset = poisson@Dataset,
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC",
  pdf = rep("Poisson", 1),
  K = 1)

summary(poissonest)

# Plot finite mixture.

plot(poissonest)

# Cluster dataset.

poissonclu <- RCLRMIX(x = poissonest, Zt = poisson@Zt)

# Plot clusters.

plot(poissonclu)
}
\keyword{clustering}