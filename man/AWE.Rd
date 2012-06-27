\name{AWE}
\alias{AWE}
\alias{AWE.default}
\alias{AWE.REBMIX}
\title{
Approximate Weight of Evidence Criterion
}
\description{
Returns the approximate weight of evidence criterion at \code{pos} for class REBMIX.
}
\usage{
\method{AWE}{REBMIX}(object, pos = 1, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\references{
Banfield JD, Raftery AE (1993). Model-Based Gaussian and Non-Gaussian Clustering.
\emph{Biometrics}, 49, 803-821.
}
\keyword{information criterion}