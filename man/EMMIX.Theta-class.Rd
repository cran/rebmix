\name{EMMIX.Theta-class}
\docType{class}
\alias{EMMIX.Theta-class}
\alias{EMMVNORM.Theta-class}
\alias{show,EMMIX.Theta-method}
\alias{show,EMMVNORM.Theta-method}
\title{Class \code{"EMMIX.Theta"}}
\description{Object of class \code{EMMIX.Theta}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EMMIX.Theta", ...)}. Accessor methods for the slots are \code{a.c(x = NULL)}, \code{a.d(x = NULL)},
\code{a.pdf(x = NULL)} and \code{a.Theta(x = NULL)}, where \code{x} stands for an object of class \code{EMMIX.Theta}. Setter methods
\code{a.theta1(x = NULL, l = numeric())}, \code{a.theta2(x = NULL, l = numeric())}, \code{a.theta3(x = NULL, l = numeric())}, 
\code{a.theta1.all(x = NULL)}, \code{a.theta2.all(x = NULL)}, \code{a.theta3.all(x = NULL)} and \code{a.w(x = NULL)} 
are provided to write to \code{Theta} slot, where \eqn{l = 1, \ldots, c}.
}
\section{Slots}{
\describe{
\item{\code{c}:}{
number of components \eqn{c > 0}. The default value is \code{1}.
}
\item{\code{d}:}{
number of dimensions.
}
\item{\code{pdf}:}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"} or \code{"vonMises"}.
}
\item{\code{Theta}:}{
a list containing \eqn{c} parametric family types \code{pdfl}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"} or circular \code{"vonMises"} defined for \eqn{0 \leq y_{i} \leq 2 \pi}.
Component parameters \code{theta1.l} follow the parametric family types. One of \eqn{\mu_{il}} for normal, lognormal, Gumbel and von Mises distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.l} follow \code{theta1.l}. One of \eqn{\sigma_{il}} for normal, lognormal and Gumbel distributions, \eqn{\beta_{il}} for Weibull and gamma distributions, \eqn{p_{il}} for binomial distribution, \eqn{\kappa_{il}} for von Mises distribution.
Component parameters \code{theta3.l} follow \code{theta2.l}. One of \eqn{\xi_{il} \in \{-1, 1\}} for Gumbel distribution.
}
\item{\code{w}:}{
a vector of length \eqn{c} containing component weights \eqn{w_{l}} summing to 1.
}
}
}
\author{Branislav Panic}
\examples{
Theta <- new("EMMIX.Theta", c = 2, pdf = c("normal", "Gumbel"))  

a.w(Theta) <- c(0.4, 0.6)

a.theta1(Theta, l = 1) <- c(2, 10)
a.theta2(Theta, l = 1) <- c(0.5, 2.3)
a.theta3(Theta, l = 1) <- c(NA, 1.0)
a.theta1(Theta, l = 2) <- c(20, 50)
a.theta2(Theta, l = 2) <- c(3, 4.2)
a.theta3(Theta, l = 2) <- c(NA, -1.0)

Theta

Theta <- new("EMMIX.Theta", c = 2, pdf = c("normal", "Gumbel", "Poisson"))  

a.w(Theta) <- c(0.4, 0.6)

a.theta1.all(Theta) <- c(2, 10, 30, 20, 50, 60)
a.theta2.all(Theta) <- c(0.5, 2.3, NA, 3, 4.2, NA)
a.theta3.all(Theta) <- c(NA, 1.0, NA, NA, -1.0, NA)

Theta

Theta <- new("EMMVNORM.Theta", c = 2, d = 3)

a.w(Theta) <- c(0.4, 0.6)

a.theta1(Theta, l = 1) <- c(2, 10, -20)
a.theta2(Theta, l = 1) <- c(9, 0, 0, 0, 4, 0, 0, 0, 1)
a.theta1(Theta, l = 2) <- c(-2.4, -15.1, 30)
a.theta2(Theta, l = 2) <- c(4, -3.2, -0.2, -3.2, 4, 0, -0.2, 0, 1)

Theta

Theta <- new("EMMVNORM.Theta", c = 2, d = 3)

a.w(Theta) <- c(0.4, 0.6)

a.theta1.all(Theta) <- c(2, 10, -20, -2.4, -15.1, 30)

a.theta2.all(Theta) <- c(9, 0, 0, 0, 4, 0, 0, 0, 1, 
  4, -3.2, -0.2, -3.2, 4, 0, -0.2, 0, 1)

Theta
}
\keyword{classes}
