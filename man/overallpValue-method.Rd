\name{overallpValue-methods}
\docType{methods}
\alias{overallpValue}
\alias{overallpValue-methods}
\alias{overallpValue,CAnDResult-method}
\title{ ~~ Methods for Function \code{overallpValue}  ~~}
\description{
This function returns the CAnD test p-value from performing a CAnD test on a set of ancestral proportions.
}
\usage{
overallpValue(object)
}
\arguments{
    \item{object}{An object of type \code{CAnDResult} }
}
\value{
The CAnD p-value.
}
\section{Methods}{
\describe{ Defined methods include:
\item{
    \code{signature(object = "CAnD")}
    }{
    Returns the CAnD test p-value from applying the CAnD test to a set of ancestral proportions for a \code{CAnD} object
    }
}
}
\seealso{
    \code{\linkS4class{CAnDResult}}-class
}
\examples{
    data(ancestries)
    overallpValue( nonParam_CAnD(ancestries[,c(2:5)]) ) 
    overallpValue( nonParam_CAnD(ancestries[,c(3:6)]) ) 
}
\author{
Caitlin McHugh
}
\keyword{methods}
