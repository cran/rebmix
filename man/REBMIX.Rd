\name{REBMIX}
\alias{REBMIX}
\title{
REBMIX Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns the REBMIX algorithm output for mixtures of conditionally independent normal, lognormal or Weibull component densities.
}
\usage{
REBMIX(Dataset = NULL, Preprocessing = NULL, D = 0.025,
       cmax = 15, InformationCriterion = "AIC", pdf = NULL, 
       K = NULL, Rmin = 0.001, ar = 0.1, Restraints = "loose")
}
\arguments{
\item{Dataset}{
a character vector containing the names of the tab delimited ASCII files, e.g. \code{"Simulated1_1.txt"} that the \emph{d}-dimensional datasets are read from.
If they do not contain an absolute path, the file names are relative to the current working directory, see also \code{\link[base]{getwd}}. Each of the \eqn{d} columns represents
one random variable. Total number of observations \eqn{n} equals the number of rows in a dataset. Use \code{sep = "\t"} and \code{eol = "\n"} in \code{write.table}.
}
\item{Preprocessing}{
a character string, giving the Preprocessing type. One of \code{"histogram"}, \code{"Parzen window"} or \code{"k-nearest neighbour"}.
}
\item{D}{
a total of positive relative deviations standing for the maximum acceptable measure of
distance between predictive and empirical densities. It satisfies the relation \eqn{0 < D \leq 1}.
The default value is \code{0.025}. However, if components with a low probability of occurrence are expected, it has to
decrease.
}
\item{cmax}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{InformationCriterion}{
a character string giving the Infromation criterion type. One of default \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"}, \code{"AICc"}, \code{"BIC"},
\code{"CAIC"}, \code{"HQC"}, \code{"MDL2"}, \code{"MDL5"}, \code{"AWE"}, \code{"CLC"}, \code{"ICL"}, \code{"PC"} or \code{"ICL-BIC"}.
}
\item{pdf}{
a character vector of length \eqn{d} containing Parametric family types. One of \code{"normal"}, \code{"lognormal"} or \code{"Weibull"}.
}
\item{K}{
a vector containing total numbers of bins \eqn{k} for the histogram and the Parzen window or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{k}. Consequently,
the REBMIX algorithm identifies it from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{k = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{k = 10 \mathrm{log_{10}}(n)} or RootN 
rule \eqn{k = 2 \sqrt{n}} can be applied to estimate the limiting total numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours.
}
\item{Rmin}{
minimum radius \eqn{0 < R_{\mathrm{min}} \leq 1}. It is advisable to keep it very close to zero. The
default value is \code{0.001}.
}
\item{ar}{
acceleration rate \eqn{0 < a_{\mathrm{r}} \leq 1}. The default value is \code{0.1} and in most cases does not have to be altered.
}
\item{Restraints}{
a character string giving the Restraints type. One of \code{"rigid"} or default \code{"loose"}.
The rigid restraints are obsolete and applicable for well separated components only. 
}
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
\item{summary}{
a data frame with additional information about Dataset, Preprocessing, \eqn{D}, \eqn{c_{\mathrm{max}}}, Information criterion type, 
\eqn{R_{\mathrm{min}}}, \eqn{a_{\mathrm{r}}}, Restraints type, optimal \eqn{c}, optimal \eqn{k}, \eqn{y_{i0}}, optimal \eqn{h_{i}},
calculation time \eqn{t_{\mathrm{c}}} in ms, information criterion \eqn{\mathrm{IC}} and log likelihood \eqn{\mathrm{log}\, L}.
}
}
\references{
Sturges HA (1926). The choice of a class interval. \emph{Journal of American Statistical Association},
21, 65-66.\cr\cr
Nagode M, Fajdiga M (1998). A General Multi-Modal Probability Density Function Suitable for 
the Rainflow Ranges of Stationary Random Processes. \emph{International Journal of Fatigue}, 20, 211-223.\cr\cr
Nagode M, Fajdiga M (2000). An Improved Algorithm for Parameter Estimation Suitable
for Mixed Weibull Distributions. \emph{International Journal of Fatigue}, 22, 75-80.\cr\cr
Nagode M, Klemenc J, Fajdiga M (2001). Parametric Modelling and Scatter Prediction of
Rainflow Matrices. \emph{International Journal of Fatigue}, 23, 525-532.\cr\cr
Nagode M, Fajdiga M (2006). An Alternative Perspective on the Mixture Estimation Problem. 
\emph{Reliability Engineering & System Safety}, 91, 388-397.\cr\cr
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Univariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}, 40(5), 876-892.\cr\cr
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Multivariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}. 40(11), 2022-2034.
}
\examples{
## Generate the complex 1 dataset.

Dataset <- "complex1.txt"

rseed <- -1

n <- c(998, 263, 1086, 487, 213, 1076, 232, 
  784, 840, 461, 773, 24, 811, 1091, 861)

Theta <- rbind(pdf = "normal",
  theta1 = c(688.4, 265.1, 30.8, 934, 561.6, 854.9, 883.7, 
  758.3, 189.3, 919.3, 98, 143, 202.5, 628, 977),
  theta2 = c(12.4, 14.6, 14.8, 8.4, 11.7, 9.2, 6.3, 10.2,
  9.5, 8.1, 14.7, 11.7, 7.4, 10.1, 14.6))

RNGMIX <- RNGMIX(Dataset = Dataset,
  rseed = rseed,
  n = n,
  Theta = Theta)

## Estimate number of components, component weights and component parameters. 

Dataset <- "complex1.txt"

Preprocessing <- "histogram"

D <- 0.0025

cmax <- 30

InformationCriterion <- "BIC"

pdf <- "normal"

k <- c(as.integer(1 + log2(sum(n))), ## Minimum k follows the Sturges rule.
  as.integer(2 * sum(n)^0.5)) ## Maximum k follows the RootN rule.

## Number of classes or nearest neighbours to be processed.

N <- as.integer(log(k[2] / (k[1] + 1)) / log(1 + 1 / k[1]))

K <- c(k[1], as.integer((k[1] + 1) * (1 + 1 / k[1])^(0:N)))

Rmin <- 0.001

ar <- 0.1

Restraints <- "loose"

REBMIX <- REBMIX(Dataset = Dataset, 
  Preprocessing = Preprocessing, 
  D = D, 
  cmax = cmax, 
  InformationCriterion = InformationCriterion, 
  pdf = pdf, 
  K = K, 
  Rmin = Rmin, 
  ar = ar, 
  Restraints = Restraints)
                 
REBMIX$w
REBMIX$Theta
REBMIX$summary                 

## Plot the finite mixture.

plot(REBMIX, npts = 1000)
}
\keyword{finite mixtures}
\keyword{parameter estimation}
\keyword{probability distributions}

