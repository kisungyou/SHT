#' my version of edist : energy is faster than mine, so bad.
#' 
#' 
#' @keywords internal
#' @noRd
test_edist <- function(X,Y,alpha=1.0,nthreads=1){
  return(energy_distance(X,Y,alpha,nthreads))
}
# 
# data(iris)
# X = as.matrix(iris[1:50,1:4])
# Y = as.matrix(iris[51:100,1:4])
# 
# library(energy)
# ori = as.double(edist(rbind(X,Y), c(50,50)))
# my1 = test_edist(X,Y)
# my2 = test_edist(X,Y,nthreads = 2)
# my4 = test_edist(X,Y,nthreads = 4)
# my8 = test_edist(X,Y,nthreads = 8)
# c(ori,my1,my2,my4,my8)
# 
# library(microbenchmark)
# microbenchmark(list=alist(ori = as.double(edist(rbind(X,Y), c(50,50))),
#                           my1 = test_edist(X,Y),
#                           my2 = test_edist(X,Y,nthreads = 2),
#                           my4 = test_edist(X,Y,nthreads = 4),
#                           my8 = test_edist(X,Y,nthreads = 8)), times=10)
