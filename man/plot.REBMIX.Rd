\name{plot.REBMIX}
\alias{plot.REBMIX}
\title{
Plots Univariate or Multivariate REBMIX Algorithm Output
}
\description{
Takes an object of class \code{REBMIX} and plots empirical and predictive univariate \eqn{d = 1} mixture densities
or pairs of empirical and predictive bivariate marginal mixture densities for \eqn{d > 1}. The empirical
densities follow the histogram, Parzen window or \emph{k}-nearest neighbour approach.
}
\usage{
\method{plot}{REBMIX}(x, pos = 1, nrow = 2, ncol = 2, npts = 200, cex = 0.8,
     fg = "black", lty = "solid", lwd = 1, pty = "m", tcl = 0.5, 
     plot.cex = 0.8, plot.pch = 19, plot.type = "p", 
     contour.drawlabels = FALSE, contour.labcex = 0.8, 
     contour.method = "flattest", contour.nlevels = 12, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{x$summary} to be ploted. The default value is \code{1}.
}
\item{nrow}{
a desired number of rows in which the empirical and predictive densities are to be ploted. The default value is \code{2}.
}
\item{ncol}{
a desired number of columns in which the empirical and predictive densities are to be ploted. The default value is \code{2}.
}
\item{npts}{
a number of points at which the predictive densities are to be ploted. The default value is \code{200}.
}
\item{cex}{
a numerical value giving the amount by which the plotting text and symbols should be magnified
relative to the default, see also \code{\link[graphics]{par}}. The default value is \code{0.8}.
}
\item{fg}{
a colour used for things like axes and boxes around plots, see also \code{\link[graphics]{par}}. 
The default value is \code{"black"}.
}
\item{lty}{
a line type, see also \code{\link[graphics]{par}}. The default value is \code{"solid"}.
}
\item{lwd}{
a line width, see also \code{\link[graphics]{par}}. The default value is \code{1}.
}
\item{pty}{
a character specifying the type of the plot region to be used. One of \code{"s"} generating a square plotting
region or \code{"m"} generating the maximal plotting region. The default value is \code{"m"}.
}
\item{tcl}{
a length of tick marks as a fraction of the height of a line of the text, see also \code{\link[graphics]{par}}. 
The default value is \code{0.5}.
}
\item{plot.cex}{
a numerical vector giving the amount by which plotting characters and symbols should be
scaled relative to the default. It works as a multiple of \code{\link[graphics]{par}("cex")}. \code{NULL} and \code{NA} are
equivalent to \code{1.0}. Note that this does not affect annotation, see also \code{\link[graphics]{plot.default}}.
The default value is \code{0.8}.
}
\item{plot.pch}{
a vector of plotting characters or symbols, see also \code{\link[graphics]{points}}. The default value is \code{19}.
}
\item{plot.type}{
a type of plot to be drawn, see also \code{\link[graphics]{plot}}. The default value is \code{"p"}.
}
\item{contour.drawlabels}{
logical. The contours are labelled if \code{TRUE}. The default value is \code{FALSE}.
}
\item{contour.labcex}{
\code{cex} for contour labelling. The default value is \code{0.8}. This is an absolute size, not a multiple of
\code{\link[graphics]{par}("cex")}.
}
\item{contour.method}{
a character string specifying where the labels will be located. The possible values
are \code{"simple"}, \code{"edge"} and default \code{"flattest"}, see also \code{\link[graphics]{contour}}.
}
\item{contour.nlevels}{
a number of desired contour levels. The default value is \code{12}.
}
\item{\dots}{
further arguments to \code{\link[graphics]{par}}.
}
}
\references{
Bishop CM (1995). \emph{Neural Networks for Pattern Recognition}. Clarendon Press, Oxford.
}
\examples{
data("wine")

colnames(wine)

## Remove the Cultivar column from the wine dataset.

winecolnames <- !(colnames(wine) \%in\% "Cultivar")

wine <- wine[ , winecolnames]

## Determine number of dimensions d and wine dataset size n.

d <- ncol(wine) 
n <- nrow(wine)

## Write the dataset without the Cultivar information into the 
## tab delimited ASCII file.

write.table(wine, file = "wine.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)

## Estimate number of components, component weights and component parameters.

Dataset <- "wine.txt"

Preprocessing <- "Parzen window"

D <- 0.025

cmax <- 15

InformationCriterion <- "ICL-BIC"

pdf <- rep("normal", d)

k <- c(as.integer(1 + log2(n)), ## Minimum k follows the Sturges rule.
  as.integer(2 * n^0.5)) ## Maximum k follows the RootN rule.

## Number of classes or nearest neighbours to be processed.

N <- as.integer(log(k[2] / (k[1] + 1)) / log(1 + 1 / k[1]))

K <- c(k[1], as.integer((k[1] + 1) * (1 + 1 / k[1])^(0:N)), k[2])

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

## Plot the finite mixture.

plot(REBMIX, nrow = 4, ncol = 5, npts = 500, pty = "s")
}
\keyword{finite mixtures}
\keyword{probability distributions}
