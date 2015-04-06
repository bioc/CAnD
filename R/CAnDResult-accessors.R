setMethod("test", "CAnDResult",
          function(object){object@test}
          )

setMethod("pValues", "CAnDResult",
          function(object){object@pValues}
          )

setMethod("overallStatistic", "CAnDResult",
          function(object){object@overallStatistic}
          )

setMethod("overallpValue", "CAnDResult",
          function(object){object@overallpValue}
          )

setMethod("BonfCorr", "CAnDResult",
          function(object){object@BonfCorr}
          )
