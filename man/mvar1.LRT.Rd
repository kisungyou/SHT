% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvar1_LRT.R
\name{mvar1.LRT}
\alias{mvar1.LRT}
\title{One-sample Simultaneous Likelihood Ratio Test of Mean and Variance}
\usage{
mvar1.LRT(x, mu0 = 0, var0 = 1)
}
\arguments{
\item{x}{a length-\eqn{n} data vector.}

\item{mu0}{hypothesized mean \eqn{\mu_0}.}

\item{var0}{hypothesized variance \eqn{\sigma_0^2}.}
}
\value{
a (list) object of \code{S3} class \code{htest} containing: \describe{
\item{statistic}{a test statistic.}
\item{p.value}{\eqn{p}-value under \eqn{H_0}.}
\item{alternative}{alternative hypothesis.}
\item{method}{name of the test.}
\item{data.name}{name(s) of provided sample data.}
}
}
\description{
Given two univariate samples \eqn{x} and \eqn{y}, it tests
\deqn{H_0 : \mu_x = \mu_0, \sigma_x^2 = \sigma_0^2 \quad vs \quad H_1 : \textrm{ not } H_0}
using likelihood ratio test.
}
\examples{
## CRAN-purpose small example
mvar1.LRT(rnorm(10))

\dontrun{
## empirical Type 1 error 
niter   = 1000
counter = rep(0,niter)  # record p-values
for (i in 1:niter){
  x = rnorm(100)  # sample x from N(0,1)
  
  counter[i] = ifelse(mvar1.LRT(x)$p.value < 0.05, 1, 0)
}

## print the result
cat(paste("\n* Example for 'mvar1.LRT'\n","*\n",
"* number of rejections   : ", sum(counter),"\n",
"* total number of trials : ", niter,"\n",
"* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
}

}
\concept{mvar}
