\name{CAnDResult-class}
\Rdversion{1.1}
\docType{class}
\alias{CAnDResult-class}

\title{Class \code{"CAnDResult"}}
\description{
Objects of this class store results from running CAnD methods. 

}
\section{Objects from the Class}{
Objects can be created by calls of \code{nonParam_CAnD} or \code{CAnD}.
}
\section{Slots}{
  \describe{
  \item{\code{test}:}{Object of class \code{"character"}, the type of CAnD test performed }
  \item{\code{pValues}:}{Object of class \code{"numeric"}, the p-values for each chromosome/chromosomal segment tested }
  \item{\code{overallStatistic}:}{Object of class \code{"numeric"}, the overall CAnD test statistic for the set }
  \item{\code{overallpValue}:}{Object of class \code{"numeric"}, the overall p-value for the set }
  \item{\code{BonfCorr}:}{Object of class \code{"logical"}, whether Bonferroni multiple testing correction was applied to the p-values for each chromosome/chromosomal segment tested }
  }
}
\section{Methods}{
No methods defined with class "CAnDResult" in the signature.
}

\author{
Caitlin McHugh
}
\examples{
showClass("CAnDResult")
}
\keyword{classes}
