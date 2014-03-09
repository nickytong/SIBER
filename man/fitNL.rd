\name{fitNL}
\alias{fitNL}
\title{Fit Negative Binomial Mixture Model}
\description{
The function fits a two-component Negative Binomial mixture model. 
}
\usage{
fitNL(y, d=NULL, model='E')
}

\arguments{  
\item{y}{A vector representing the transformed data that follows the normal mixture distribution.}
\item{d}{A vector of the same length as y representing the normalization constant to be applied to the data. }
\item{model}{Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. } 
}

\details{
This function calls the mclust package to fit the 2-component normal mixture. 
}

\value{
A vector consisting parameter estimates of mu1, mu2, phi1, phi2, pi1, logLik and BIC. 
}

\references{
Wang, J.,Wen, S., Symmans,W., Pusztai, L., and Coombes, K. (2009). The bimodality
index: a criterion for discovering and ranking bimodal signatures from cancer gene
expression profiling data. Cancer informatics, 7, 199.

Fraley, C. and Raftery, A. (2002). Model-based clustering, discriminant analysis, and
density estimation. Journal of the american statistical association, 97(458), 611:631.
Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes
Using RNAseq Data. Bioinformatics, submitted.
}

\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (kcoombes@mdanderson.org)
}
\seealso{
\link{SIBER}
\link{fitLN}
\link{fitNB}
\link{fitGP}
}


\examples{
# artificial microarray data from normal distribution
set.seed(1000)
dat <- rnorm(100, mean=5, sd=1)
fitNL(y=dat)
}