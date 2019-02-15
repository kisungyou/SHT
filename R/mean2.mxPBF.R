#' mxPBF
#' 
#' 
#' @export
mean2.mxPBF <- function(x){
  y = lgamma(x)
  z = testcpp_lgamma(x)
  
  output = list()
  output$R = y
  output$C = z
  return(output)
}