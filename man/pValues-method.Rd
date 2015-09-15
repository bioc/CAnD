\name{pValues-methods}
\docType{methods}
\alias{pValues}
\alias{pValues-methods}
\alias{pValues,CAnDResult-method}
\title{ ~~ Methods for Function \code{pValues}  ~~}
\description{
This function returns each p-value calculated from the CAnD test on all chromosomes/chromosomal segments.
}
\usage{
pValues(object)
}
\arguments{
    \item{object}{An object of type \code{CAnDResult} }
}
\value{
A vector of p-values from the CAnD test.
}
\section{Methods}{
\describe{ Defined methods include:
\item{
    \code{signature(object = "CAnD")}
    }{
    Returns the p-values for each chromosome/chromosomal segment calculated using the CAnD test for a \code{CAnD} object
    }
}
}
\seealso{
    \code{\linkS4class{CAnDResult}}-class
}
\examples{
    data(ancestries)
    pValues( CAnD(ancestries[,c(5:12)]) ) 
}
\author{
Caitlin McHugh
}
\keyword{methods}
