\name{wine}
\docType{data}
\alias{wine}
\title{Wine Recognition Data}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three
different cultivars (1-3). The analysis determined the quantities of 13 constituents: alcohol, malic acid, ash, alcalinity of ash,
magnesium, total phenols, flavanoids, nonflavanoid phenols, proanthocyanins, colour intensity, hue, OD280/OD315 of diluted wines,
and proline found in each of the three types of the wines. The number of instances in classes 1 to 3 is 59, 71 and 48, respectively.
}
\usage{
wine
}
\format{
\code{wine} is a data frame with 178 cases (rows) and 14 variables (columns) named \code{Alcohol}, \code{Malic.Acid}, \code{Ash}, \code{Alcalinity.of.Ash},
\code{Magnesium}, \code{Total.Phenols}, \code{Flavanoids}, \code{Nonflavanoid.Phenols}, \code{Proanthocyanins}, \code{Color.Intensity}, \code{Hue},
\code{OD280.OD315.of.Diluted.Wines}, \code{Proline} and \code{Cultivar}.
}
\source{
Asuncion A, Newman DJ (2007). UCI Machine Learning Repository. \url{http://www.
ics.uci.edu/~mlearn/MLRepository.html}.
}
\references{
Roberts SJ, Everson R, Rezek I (2000). Maximum Certainty Data Partitioning. \emph{Pattern
Recognition}, 33, 833-839.
}
\examples{
data("wine")

colnames(wine)

## Split the wine dataset into three subsets for the three Cultivars
## and remove the Cultivar column from the datasets.

winecolnames <- !(colnames(wine) \%in\% "Cultivar")

wine1 <- wine[wine$Cultivar == 1, winecolnames]
wine2 <- wine[wine$Cultivar == 2, winecolnames]
wine3 <- wine[wine$Cultivar == 3, winecolnames]

wine <- wine[ , winecolnames]

## Write the datasets without the Cultivar information into the tab 
## delimited ASCII files.

write.table(wine, file = "wine.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine1, file = "wine1.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine2, file = "wine2.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine3, file = "wine3.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
}
\keyword{datasets}
