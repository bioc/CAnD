
##' Plots CAnD p-values for each chromosome/chromosomal region
##'
##' Creates a plot of all p-values for each chromosome or 
##' chromosomal region.
##' @title Create a Plot of P-Values for Each Chromosome or 
##' Chromosomal Region
##' @param set An object of class \code{CAnDResult}. 
##' @param title A character string containing the title of the plot.
#' Default is "", a blank title.
##' @param xlab A character vector with the label for the x-axis on
#'  the plot. Default is \code{Chromosome}.
##' @param ylab A character vector holding the label for the y-axis 
#' on the plot. Default is \code{-log10(Bonferroni PValue)} 
#' or \code{-log10(PValue)}, 
#' depending on whether Bonferroni correction was used.
##' @param ... Further arguments to be passed to the plotting methods,
#' such as graphical parameters.
##' @return Creates a plot.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @examples
#' data(ancestries)
#' euroEsts <- ancestries[,c(seq(from=2,to=24))]
#' res <- CAnD(euroEsts)
#' #plotPvals(res,main="CAnD P-Values")
##' @export
plotPvals <- function(set, title="", xlab = "Chromosome", 
                      ylab = "-log10(PValue)", ...)
{
  if(!is(set,"CAnDResult")){
    stop("set must be a CAnDResult object.")}
  
  if(BonfCorr(set)){ylab="-log10(Bonferroni PValue)"}
  
  pvalRes <- data.frame("pval"=as.numeric(pValues(set)))
  pvalRes$numPs <- 1:nrow(pvalRes)
  
  ggplot(pvalRes, aes(x=numPs,y=-log10(pval), ...)) +
    geom_point() + xlab(xlab) + ylab(ylab) + ggtitle(title)
}