\name{overallStatistic-methods}
\docType{methods}
\alias{overallStatistic}
\alias{overallStatistic-methods}
\alias{overallStatistic,CAnDResult-method}
\title{ ~~ Methods for Function \code{overallStatistic}  ~~}
\description{
This function returns the calculated CAnD test statistic from applying the CAnD test to a set of ancestral proportions.
}
\usage{
overallStatistic(object)
}
\arguments{
    \item{object}{An object of type \code{CAnDResult} }
}
\value{
The CAnD test statistic stored in the object.
}
\section{Methods}{
\describe{Defined methods include:
\item{
    \code{signature(object = "CAnD")}
    }{
    Returns the statistic calculated from performing the CAnD test on a set of ancestral proportions for a \code{CAnD} object
    }
}
}
\seealso{
    \code{\linkS4class{CAnDResult}}-class
}
\examples{
    data(ancestries)
    overallStatistic( CAnD(ancestries[,c(2:5)]) ) 
}
\author{
Caitlin McHugh
}
\keyword{methods}
