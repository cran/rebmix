\name{RNGMIX}
\alias{RNGMIX}
\title{
Random Univariate or Multivariate Finite Mixture Generation
}
\description{
Returns univariate or multivariate random datasets for mixtures of conditionally independent normal, lognormal or Weibull component densities.
}
\usage{
RNGMIX(Dataset = NULL, rseed = -1, n = NULL, Theta = NULL)
}
\arguments{
\item{Dataset}{
a character vector containing the names of the tab delimited ASCII files, e.g. \code{"Simulated1_1.txt"} that \emph{d}-dimensional datasets are written in.
If they do not contain an absolute path, the file names are relative to the current working directory, see also \code{\link[base]{getwd}}.
}
\item{rseed}{
set the random seed to any negative integer value to initialize the sequence. The first file in \code{Dataset} corresponds to it.
For each next file the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{n}{
a vector containing total number of observations in classes \eqn{n_{l}}, where total number of observations \eqn{n = \sum_{l = 1}^{c} n_{l}}.
}
\item{Theta}{
a \eqn{3 \times d \times c} matrix. The first, fourth, etc. rows contain \eqn{c} Parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"} or \code{"Weibull"}.
The second, fifth, etc. rows contain \eqn{c} component parameters \code{theta1.i}. One of \eqn{\mu_{il}} for normal and lognormal distributions or \eqn{\theta_{il}} for Weibull distribution.
The third, sixth, etc. rows contain \eqn{c} component parameters \code{theta2.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions or \eqn{\beta_{il}} for Weibull distribution.
}
}
\details{
RNGMIX is based on the "Minimal" random number generator \code{ran1} of Park and Miller with the Bays-Durham shuffle and added safeguards that returns a uniform random deviate between 0.0 and 1.0 
(exclusive of the endpoint values).
}
\value{
\item{w}{
a data frame containing \eqn{c} component weights \eqn{w_{l}} summing to 1.
}
\item{Theta}{
a \eqn{3\times d \times c} data frame. The first, fourth, etc. rows contain \eqn{c} Parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"} or \code{"Weibull"}.
The second, fifth, etc. rows contain \eqn{c} component parameters \code{theta1.i}. One of \eqn{\mu_{il}} for normal and lognormal distributions or \eqn{\theta_{il}} for Weibull distribution.
The third, sixth, etc. rows contain \eqn{c} component parameters \code{theta2.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions or \eqn{\beta_{il}} for Weibull distribution.
}
}
\references{
Press WH, Teukolsky SA, Vetterling VT, Flannery BP (1992). \emph{Numerical Recipes in C: The Art
of Scientific Computing}. Cambridge University Press, Cambridge.
}
\examples{
## Simulated1 dataset consists of n = 625 four-dimensional observations  
## obtained by generating samples separately from each of the five normal 
## distributions.

Dataset <- paste("Simulated1_", 1:25, ".txt", sep = "")

rseed <- -1

n <- c(75, 100, 125, 150, 175)

Theta <- rbind(pdf1 = rep("normal", 5),
  theta1.1 = c(10, 8.5, 12, 13, 7),
  theta2.1 = c(1, 1, 1, 2, 3),
  pdf2 = rep("normal", 5),
  theta1.2 = c(12, 10, 14, 15, 9),
  theta2.2 = c(1, 1, 1, 2, 3),
  pdf3 = rep("normal", 5),
  theta1.3 = c(10, 8.5, 12, 7, 13),
  theta2.3 = c(1, 1, 1, 2, 3),
  pdf4 = rep("normal", 5),
  theta1.4 = c(12, 10.5, 14, 9, 15),
  theta2.4 = c(1, 1, 1, 2, 3))

RNGMIX <- RNGMIX(Dataset = Dataset,
  rseed = rseed,
  n = n,
  Theta = Theta)

RNGMIX
}
\keyword{finite mixtures}
\keyword{random number generation}
\keyword{probability distributions}
