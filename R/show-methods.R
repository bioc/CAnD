## show
setMethod("show", "CAnDResult",
function(object){
    cat(paste0("CAnD results for ", object@test, " test\n"))
    cat("Bonferroni correction was ")
    if(!object@BonfCorr){cat("not ")}
    cat("used\n")
    cat(paste0("p-values = ", signif(object@pValues,3), "\n") )
    cat(paste0("observed CAnD statistic = ", 
        signif(object@overallStatistic,3), "\n") )
    cat(paste0("calculated CAnD p-value = ", 
        signif(object@overallpValue,3), "\n"))
}
)