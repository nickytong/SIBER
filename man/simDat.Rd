\name{simDat}
\alias{simDat}
\alias{parList}
\alias{dataList}

\docType{data}

\title{Simulated Data From 2-component Mixture Models}

\description{Data from 2-component mixture models (NB, GP and LN) is simulated
with the true parameters given for testing and illustration purpose.}

\format{
  The data frame contains the following data objects:
  \describe{
    \item{parList}{A list of true parameters. There are three named elements 
	(NB, GP and LN) corresponding to
	the parameters used to simulate gene expression data from NB, GP and LN mixture models. Each
	element is a 6 by 5 matrix giving the true parameters generating the simulated data. }
    \item{dataList}{A list of matrices for simulated gene expression data. 
	There are three named elements (NB, GP and LN) corresponding to
	the simulate gene expression data from NB, GP and LN mixture models. Each
	element is a 6 by 200 matrix. That is, 6 genes (rows) are simulated with 200 samples (columns). 
	The first 3 genes in each matrix are from 2-component mixture model while the last 3 genes
	are from 0-inflated models. }
  }
}

\usage{
data(simDat)
}

\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (kcoombes@mdanderson.org)
}
\seealso{
\link{SIBER}
\link{fitNB}
\link{fitGP}
\link{fitLN}
\link{fitNL}
}


\examples{
library(SIBER)
data(simDat)
sapply(parList, dim)
sapply(dataList, dim)

}

\keyword{datasets}
