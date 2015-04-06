\name{test-methods}
\docType{methods}
\alias{test}
\alias{test-methods}
\alias{test,CAnDResult-method}
\title{ ~~ Methods for Function \code{test}  ~~}
\description{
This function returns the type of CAnD test applied to a set of ancestral proportions.
}
\usage{
test(object)
}
\arguments{
    \item{object}{An object of type \code{CAnDResult} }
}
\value{
Either `non-parametric' or `parametric,' depending on which test was run to produce the object.
}
\section{Methods}{
\describe{ Defined methods include:
 \item{
    \code{signature(object = "CAnD")}
    }{
    Returns the type of CAnD test applied to calculate the p-values for a \code{CAnD} object
    }
}
}
\seealso{
    \code{\linkS4class{CAnDResult}}-class
}
\examples{
    data(ancestries)
    test( nonParam_CAnD(ancestries[,c(2:5)]) ) 
    test( CAnD(ancestries[,c(2:5)]) ) 
}
\author{
Caitlin McHugh
}
\keyword{methods}
