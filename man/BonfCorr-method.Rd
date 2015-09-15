\name{BonfCorr-methods}
\docType{methods}
\alias{BonfCorr}
\alias{BonfCorr-methods}
\alias{BonfCorr,CAnDResult-method}
\title{ ~~ Methods for Function \code{BonfCorr}  ~~}
\description{
This function returns whether the Bonferroni multiple 
testing correction was applied in calculating the final 
p-values by chromosome/chromosomal segment.
}
\usage{
BonfCorr(object)
}
\arguments{
    \item{object}{An object of type \code{CAnDResult} }
}
\value{
A logical indicator of whether the Bonferroni multiple testing correction was used.
}
\section{Methods}{
\describe{Defined methods include:
\item{
    \code{signature(object = "CAnD")}
    }{
    Returns the whether the Bonferroni multiple testing 
correction was applied in calculating the final 
p-values for a \code{CAnD} object
    }
}
}
\seealso{
    \code{\linkS4class{CAnDResult}}-class
}
\examples{
    data(ancestries)
    BonfCorr( CAnD(ancestries[,c(2:5)]) ) 
}
\author{
Caitlin McHugh
}
\keyword{methods}
