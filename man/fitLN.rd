\name{fitLN}
\alias{fitLN}
\title{Fit Log Normal Mixture Model}
\description{
The function fits a two-component log normal mixture model. 
}
\usage{
fitLN(y, base=10, eps=10, d=NULL, model='E', zeroPercentThr=0.2, logLikToLN=TRUE)
}

\arguments{  
\item{y}{A vector representing the RNAseq raw count. }
\item{base}{The logarithm base defining the parameter estimates in the logarithm scale. This is also the base of log transformation applied to the data. } 
\item{eps}{A scalar to be added to the count data to avoid taking logarithm of zero. } 
\item{d}{A vector of the same length as y representing the normalization constant to be applied to the data. For the LN model, the original data would be devided by this vector. }
\item{model}{Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. } 
\item{zeroPercentThr}{A scalar specifying the minimum percent of zero counts needed when fitting a zero-inflated 
log normal model. This
parameter is used to deal with zero-inflation in RNAseq count data. When the percent of zero exceeds this threshold, 
1-comp mixture LN model is used to estimate mu and sigma from nonzero count.} 
\item{logLikToLN}{logical indicating if the log likelihod is defined on the transformed value or the orginal value from log noral distribution. } 
}

\details{
The parameter estimates from log normal mixture is obtained by taking logarithm and fit normal mixture. We use
mclust package to obtain parameter estimates of normal mixture model. In particular, \eqn{log_{base}(\frac{y+eps}{d})} is used to
fit to normal mixture model. 

With this function, three models can be fitted: (1) log normal mixture with equal variance (E model); 
(2) Generalized Poisson mixture with unequal variance (V model); (3) 0-inflated log normal model.
The 0-inflated log normal has the following density function:

\eqn{P(Y=y)=\pi D(y) + (1-\pi)LN(\mu, \sigma)} where D is the point mass at 0 while \eqn{LN(\mu, \sigma)} is the density
of log normal distribution with mean \eqn{\mu} and variance \eqn{\sigma^2}. 

The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
rule overrides the model argument (E or V) when observed percentae of zero count exceeds the threshold.

}

\value{
A vector consisting parameter estimates of mu1, mu2, sigma1, sigma2, pi1, logLik and BIC. 
For 0-inflated model, mu1=sigma1=0.

}

\references{
Wang, J.,Wen, S., Symmans,W., Pusztai, L., and Coombes, K. (2009). The bimodality
index: a criterion for discovering and ranking bimodal signatures from cancer gene
expression profiling data. Cancer informatics, 7, 199.

Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes
Using RNAseq Data. Bioinformatics, submitted.
}

\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (kcoombes@mdanderson.org)
}
\seealso{
\link{SIBER}
\link{fitNB}
\link{fitGP}
\link{fitNL}
}


\examples{
# artificial RNAseq data from negative binomial distribution
set.seed(1000)
dat <- rnbinom(100, mu=1000, size=1/0.2)
fitLN(y=dat)
}