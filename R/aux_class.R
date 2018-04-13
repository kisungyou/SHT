#' @keywords internal
#' @noRd
hypothesis <- function(method, statistic, alpha, pvalue, 
                       alternative, conclusion){
  structure(
    list(method       = method,
         statistic    = statistic,
         p.value      = pvalue,
         significance = alpha,
         alternative  = alternative,
         conclusion   = conclusion
         ),
    class = c("hypothesis")
  )
}