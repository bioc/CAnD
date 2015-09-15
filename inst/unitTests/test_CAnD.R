
# unit tests for the CAnD package

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
   
  asian[34,3] <- NA
  res <- tryCatch(CAnD(asian),warning=conditionMessage)
  checkIdentical("NA values will be excluded from the analysis.", res)
  
  asian[34,3] <- "J"
  checkException(CAnD(asian))
  
  asian[,3] <- "J"
  checkException(CAnD(asian))
}



