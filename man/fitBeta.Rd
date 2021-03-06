% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Beta.R
\name{fitBeta}
\alias{fitBeta}
\title{Fit Beta Mixture Model
The function fits a two-component Beta mixture model.
This function directly maximize the log likelihood function through optimization.}
\usage{
fitBeta(y, inits = NULL, model = "V")
}
\arguments{
\item{y}{A vector representing the RNAseq raw count.}

\item{inits}{Initial value to fit the mixture model. A vector with elements alpha1, beta1, alpha2, beta2 and pi1.}

\item{model}{Character specifying E or V model. E model fits the mixture model with equal dispersion phi while V model doesn't put any constraint. Currently only V model is implemented.}
}
\value{
A vector consisting parameter estimates of alpha1, beta1, alpha2, beta2 and pi1, logLik and BIC.
}
\description{
Fit Beta Mixture Model
The function fits a two-component Beta mixture model.
This function directly maximize the log likelihood function through optimization.
}
\examples{
set.seed(100)
y=c(rbeta(100,1,4),rbeta(200,4,1))
fitBeta(y=y)
SIBER(y, model='Beta')
}

