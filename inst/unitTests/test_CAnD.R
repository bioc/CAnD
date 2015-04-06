
# unit tests for the CAnD package

test_nonParamCAnD <- function(){
  data(ancestries)
  asianCols <- grep("Asian_",colnames(ancestries))
  asian <- ancestries[,c(1,asianCols)]
  checkTrue(class(nonParam_CAnD(asian))=="CAnDResult")
  checkException(nonParam_CAnD(asian[,2]))
  
  res <- nonParam_CAnD(asian)
  checkTrue(BonfCorr(res)==TRUE)
  checkTrue(test(res)=="non-parametric")
  checkTrue(overallpValue(res)<=1.0)
  checkTrue(all(pValues(res)<=1.0))
  
  res <- nonParam_CAnD(asian,FALSE)
  checkTrue(BonfCorr(res)==FALSE)
  checkTrue(test(res)=="non-parametric")
  
  asian[34,3] <- NA
  res <- tryCatch(nonParam_CAnD(asian),warning=conditionMessage)
  checkIdentical("NA values will be excluded from the analysis.", res)
  
  asian[34,3] <- "J"
  checkException(nonParam_CAnD(asian))
  
  asian[,3] <- "J"
  checkException(nonParam_CAnD(asian))
}

test_CAnD <- function(){
  data(ancestries)
  asianCols <- grep("Asian_",colnames(ancestries))
  asian <- ancestries[,c(1,asianCols)]
  checkTrue(class(CAnD(asian))=="CAnDResult")
  checkException(CAnD(asian[,2]))
  
  res <- CAnD(asian)
  checkTrue(BonfCorr(res)==TRUE)
  checkTrue(test(res)=="parametric")
  checkTrue(overallpValue(res)<=1.0)
  checkTrue(all(pValues(res)<=1.0))
  
  res <- CAnD(asian,FALSE)
  checkTrue(BonfCorr(res)==FALSE)
  checkTrue(test(res)=="parametric")
  
  res <- tryCatch(CAnD(asian[1:8,]),warning=conditionMessage)
  checkIdentical("The number of samples may be too small for assumptions to hold.
            Consider running 'nonParam_CAnD' instead.", res)
  
  asian[34,3] <- NA
  res <- tryCatch(CAnD(asian),warning=conditionMessage)
  checkIdentical("NA values will be excluded from the analysis.", res)
  
  asian[34,3] <- "J"
  checkException(CAnD(asian))
  
  asian[,3] <- "J"
  checkException(CAnD(asian))
}



