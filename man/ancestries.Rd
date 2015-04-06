\name{ancestries}
\docType{data}
\alias{ancestries}
\title{Sample Chromosomal Proportions of Ancestral Subpopulations to Use for CAnD Testing}
\description{
	This dataset provides simulated ancestral proportions by chromosomes for 50 samples. The samples are assumed to have three different ancestral proportions, denoted by \code{Euro_}, \code{Afr_} and \code{Asian_}. Furthermore, these proportions are included as averages across each chromosome 1-22 and the X chromosome.
}
\usage{ancestries}
\format{A data.frame containing 50 rows and 70 columns:
\tabular{ll}{
	IID    \tab unique sample id\cr
	Euro_AA  \tab simulated proportion European ancestry on each chromosome, where \code{AA} ranges from 1-22, X\cr
	Afr_AA \tab simulated proportion African ancestry on each chromosome, where \code{AA} ranges from 1-22, X\cr
	Asian_AA \tab simulated proportion Asian ancestry on each chromosome, where \code{AA} ranges from 1-22, X
}
}
\value{
A character value of the name of the dataset.
}
\keyword{datasets}