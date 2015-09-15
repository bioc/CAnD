
##' Perform the CAnD test on a set of ancestry proportions 
#' estimated for a particular ancestral subpopulation of interest
#'
##' @title Peform the CAnD Test
##' @param chrAncest A data.frame holding the ancestral proportions; 
#' each row corresponds to a sample and each column corresponds to a 
#' chromosomal/chromosomal segment ancestry proportion. 
#' Note: only include the proportions for one ancestral population at a time.
##' @param bonfCorr A logical argument indicating whether the p-value 
#' should be corrected for multiple testing using Bonferroni 
#' correction. The default is \code{TRUE}.
##' @return A \code{CAnDResult} object holding the p-value for each 
#' chromosome/chromosomal segment, the overall CAnD p-value, the
#' CAnD statistic and whether the Bonferroni multiple testing 
#' correction was used.
##' @references  McHugh, C., Brown, L., Thornton, T.
#' Detecting heterogeneity in population structure across chromosomes in admixed populations.
#' Manuscript in Preparation.
##' @author Caitlin McHugh
#' \email{mchughc@@uw.edu}
##' @examples
#' data(ancestries)
#' euroCols <- grep("Euro",colnames(ancestries))
#' euro <- ancestries[,euroCols]
#' res <- CAnD(euro)
#' res
##' @export 

CAnD <- function(chrAncest,bonfCorr=TRUE){
  
  if(!(is(chrAncest,c("matrix","data.frame")))){
    stop("chrAncest must be a matrix or data.frame.")
  }
  
  if(!is.logical(bonfCorr)){
    stop("bonfCorr must be a logical argument.")
  }
    
  numChrs <- ncol(chrAncest)
  
  if(ncol(chrAncest)<=1){
    stop("chrAncest must have at least two columns of 
         ancestry proportions for testing.")
  }
  
  if(!all(vapply(chrAncest,is.numeric,logical(1)))){
    stop("Ancestry proportions must be numeric values.") }
  
  if(!sum(is.na(chrAncest))==0){
    warning("NA values will be excluded from the analysis.") }
    
  diff_means <- getDiffMatrices(chrAncest,diff=FALSE)
  
  pairedTtest <- function(x) 
  {
    t.test(chrAncest[,x], diff_means[,x],paired=TRUE)$p.value
  }
  pval <- vapply(seq_len(numChrs), pairedTtest, 0)
  
  # calculate correlation between cand statistics
  combinedRes <- calc_combP(chrAncest)
  
  if(bonfCorr){ pval <- pval*numChrs }
  
  pval <- ifelse(pval>1,1,pval)
  names(pval) <- colnames(chrAncest)

  return( new("CAnDResult",
              test = "parametric",
              pValues = pval,
              overallStatistic = as.numeric(combinedRes["statistic"]),
              overallpValue = as.numeric(combinedRes["pvalue"]),
              BonfCorr = bonfCorr) )  
}

