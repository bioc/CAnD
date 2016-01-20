
##' Calculate the combined CAnD test statistic p-value on a set of 
#' ancestry proportions estimated for a particular ancestral 
#' subpopulation of interest
#'
##' @title Calculate the Combined CAnD Test Statistic P-value
##' @param chrAncest A data.frame holding the ancestral proportions; 
#' each row corresponds to a sample and each column corresponds to a 
#' chromosomal/chromosomal segment ancestry proportion. 
#' Note: only include the proportions for one ancestral population at a time.
##' @return A \code{vector} of length two where `statistic' is the
##' combined CAnD statistic and `pvalue' is it's corresponding p-value, 
#' where the combined statistic is combined over all 
#' chromosomes/chromosomal segments included in \code{chrAncest}.
##' @references  McHugh, C., Brown, L., Thornton, T.
#' Detecting heterogeneity in population structure across chromosomes in admixed populations.
#' Manuscript in Preparation.
##' @author Caitlin McHugh
#' \email{mchughc@@uw.edu}
##' @examples
#' data(ancestries)
#' euroCols <- grep("Euro",colnames(ancestries))
#' euro <- ancestries[,euroCols]
#' res <- calc_combP(euro)
#' res
##' @export 

calc_combP <- function(chrAncest){
  
  if(!(is(chrAncest,c("matrix","data.frame")))){
    stop("chrAncest must be a matrix or data.frame.")
  }
  
  if(ncol(chrAncest)<=1){
    stop("chrAncest must have at least two columns of 
         ancestry proportions for testing.")
  }
  
  if(!all(vapply(chrAncest,is.numeric,logical(1)))){
    stop("Ancestry proportions must be numeric values.") }
  
  if(!sum(is.na(chrAncest))==0){
    warning("NA values will be excluded from the analysis.") }
  
  # get the w_cc',i for each individual i
  numChrs <- ncol(chrAncest)
  n <- nrow(chrAncest)
  
  pairChrs <- expand.grid(1:numChrs,1:numChrs) # all pairwise combos of chrs
  # remove rows where chrs are =
  cc <- pairChrs[,1]==pairChrs[,2]
  pairChrs <- pairChrs[!cc,]
  
  chrAncest$abar <- rowMeans(chrAncest)
  
  paircorr <- function(x){
    tmp <- cbind(x[pairChrs[,1]]-x["abar"],x[pairChrs[,2]]-x["abar"])
    mean(tmp[,1]*tmp[,2])
  }
  pairw <- apply(chrAncest,1,paircorr)
  # pairw_{cc',i} should be zero under the null, since we are adj for the mean ancestry w/in an individ
  # this is relative to an individ, so there will be no correlation
  # if this is relative to the population itself, there will be correlation
  
  chrAncest$abar <- NULL 
  
  # now, calculate sig2_c for each pair of chromosomes
  getTstat <- function(x,n){
    exclChr <- chrAncest[,-x]
    d <- chrAncest[,x] - rowMeans(exclChr)
    dbar <- mean(d)
    
    sig2 <- sum((d-dbar)^2)/(n*(n-1))
    tstat <- dbar/sqrt(sig2)
    
    return(c(sig2,tstat))
  }
  res <- vapply(seq_len(numChrs),getTstat,c(sig2=0,tstat=0),n=n) 
  # first row of this is chr level variance
  # second row is t statistic
  
  sig2.i <- apply(chrAncest,1,var) # this is individual level variance
  
  pairChrs <- expand.grid(1:numChrs,1:numChrs) # get all pairwise combos of chrs again
  denom <- n^2*sqrt(res["sig2",pairChrs[,1]])*sqrt(res["sig2",pairChrs[,2]])
  sig.matrix.values <- (1/denom)*sum((numChrs/(numChrs-1)^2)*(pairw-sig2.i))
  sig.matrix <- matrix(sig.matrix.values,byrow=TRUE,nrow=numChrs)
  diag(sig.matrix) <- 1
  
  # calculate new stat
  newstat <- res["tstat",]%*%solve(sig.matrix)%*%res["tstat",]
  pval <- pchisq(newstat,df=(numChrs-1),lower.tail=FALSE)
  return(c(statistic=newstat,pvalue=as.numeric(pval)))
}
