% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean1_ttest.R
\name{mean1.ttest}
\alias{mean1.ttest}
\title{One-sample Student's t-test for Univariate Mean}
\usage{
mean1.ttest(x, mu0 = 0, alternative = c("two.sided", "less", "greater"))
}
\arguments{
\item{x}{a length-\eqn{n} data vector.}

\item{mu0}{hypothesized mean \eqn{\mu_0}.}

\item{alternative}{specifying the alternative hypothesis.}
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
Given an univariate sample \eqn{x}, it tests
\deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
using the procedure by Student (1908).
}
\examples{
## empirical Type 1 error 
niter   = 1000
counter = rep(0,niter)  # record p-values
for (i in 1:niter){
  x = rnorm(10)         # sample from N(0,1)
  counter[i] = ifelse(mean1.ttest(x)$p.value < 0.05, 1, 0)
}

## print the result
cat(paste("\n* Example for 'mean1.ttest'\n","*\n",
"* number of rejections   : ", sum(counter),"\n",
"* total number of trials : ", niter,"\n",
"* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))

}
\references{
\insertRef{student_probable_1908}{SHT}

\insertRef{student_probable_1908a}{SHT}
}
\concept{mean_univariate}
