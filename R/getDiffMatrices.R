##' A helper function to calculate the mean ancestry proportion for a 
##' given subpopulation, excluding each chromosome/chromosomal segment
##' in turn.
##'
##' This function calculates the mean ancestry proportion of a given
##' subpopulation excluding each chromosome in turn.
##' @title Calculate the Mean Ancestry Proportion Excluding Each 
##' Chromosome/Chromosomal Segment in Turn
##' @param chrAncest A data.frame holding the ancestral proportions; 
#' each row corresponds to a sample and each column corresponds to a 
#' chromosomal/chromosomal segment ancestry proportion. 
#' Note: only include the proportions for one ancestral population at a time.
##' @param diff A logical argment indicating whether the difference 
##' between the pooled mean and the chromosomal mean should be returned,
##' or whether simply the pooled mean should be returned.
##' @return A matrix of chromosomal ancestry differences.
##' @author Caitlin McHugh
##' \email{mchughc@@uw.edu}
##' @export
getDiffMatrices <- function(chrAncest,diff=TRUE){
  
  numChrs <- ncol(chrAncest)
  n <- nrow(chrAncest)
  
  rowmns <- function(x){
    rowMeans(chrAncest[-c(x+1)])
  }
  res <- vapply(seq_len(numChrs),rowmns,rep(0,n))
  #  res holds the mean ancestry prop excluding each of the chrs in turn
  
  if(diff){diff_means <- res-chrAncest
  }else{ diff_means <- res }
  # diff_means holds the mean diff between each chr and an 
  # average of all the others
  
  return(diff_means)
  
}