#' # S3 class \emph{hypothesis}
#' # 
#' 
#' #'
#' #' @keywords internal
#' #' @noRd
#' hypothesis <- function(method, statistic, alpha, pvalue,
#'                        alternative, conclusion){
#'   structure(
#'     list(method       = method,
#'          statistic    = statistic,
#'          p.value      = pvalue,
#'          significance = alpha,
#'          alternative  = alternative,
#'          conclusion   = conclusion
#'          ),
#'     class = c("hypothesis")
#'   )
#' }

# generic functions -------------------------------------------------------
#' class S3 method addition
#' print <- function(x){
#'  UseMethod("print")
#' }
#' #' S3 class 'hypothesis' print function.
#' print.hypothesis <- function(x){
#'   width = getOption("width")
#'   
#'   mname = x$method
#'   
#'   vert = "\n"
#'   ws     = rep(" ", max(floor((width - nchar(mname))/2),1))
#'   wshalf = rep(" ", min(abs(floor((width - nchar(mname))/4)),1))
#'   
#'   cat(vert)
#'   # cat(ws, mname, vert, sep="") # method title
#'   cat(wshalf,wshalf, mname, vert, sep="") # method title
#'   cat(vert)
#'   cat(wshalf,"* test statistic : ", x$statistic,vert,sep="")
#'   cat(wshalf,"* alternative    : ", x$alternative,vert,sep="")
#'   cat(wshalf,"* significance   : ", x$significance,vert,sep="")
#'   cat(wshalf,"* p-value        : ", x$p.value,vert,sep="")
#'   cat(wshalf,"* conclusion     : ", x$conclusion,vert,sep="")
#' }