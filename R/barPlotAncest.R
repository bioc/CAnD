
##' Plots ancestry proportion estimates for each sample
##'
##' Creates a barplot of ancestry proportions for each sample
##' for a given chromosome or chromosomal region.
##' 
##' @title Create a Barplot of Ancestry Proportion Estimates for 
##' Every Sample and a Given Chromosome or Chromosomal Region
##' @param set A \code{data.frame} with columns of the proportion 
##' ancestry for a given chromosome or chromosomal region, and one
##' row per sample (bar). 
##' @param order A \code{logical} argument determining whether the 
##' samples should be ordered in increasing proportion of the first
##' ancestry. Default is TRUE.
##' @param title A character string containing the title of the plot.
#' Default is "", a blank title.
##' @param xlab A character vector with the label for the x-axis on
#'  the plot. Default is \code{Sample}.
##' @param ylab A character vector holding the label for the y-axis 
#' on the plot. Default is \code{Ancestry Proportion}.
##' @param ... Further arguments to be passed to the plotting methods,
#' such as graphical parameters.
##' @return Creates a plot.
##' @author Caitlin McHugh \email{mchughc@@uw.edu}
##' @examples
#' data(ancestries)
#' chr1 <- ancestries[,c("Euro_1","Afr_1","Asian_1")]
#' #barPlotAncest(chr1,title="Chr 1 Ancestry Proportions")
##' @export
barPlotAncest <- function(set, order=TRUE, title="", xlab = "Sample", 
                      ylab = "Ancestry Proportion", ...)
{
  if(!is.data.frame(set)){
    stop("set must be a data.frame object.")}
  
  if(order){  set <- set[order(set[,1]),] }
  
  set$id <- 1:nrow(set)
  set <- melt(set,id="id")
  
  ggplot(set,aes(id,value,fill=variable)) + geom_bar(stat="identity") +
    xlab(xlab) + ylab(ylab) + ggtitle(title)
}
