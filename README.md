# SIBER

SIBER (systematic identification of bimodally expressed genes using RNAseq data) is a methodology I developed to identify aberrant gene expression pattern which falls into the realm of outlier detection in statistics.  

The **SIBER** R package achieved this by calculating the Bimodality Index for both microarray and RNAseq data using
a mixture of Gaussian, Negative Bimomial and Generalized Poisson Models. 

Our approach compared favorably with alternative methods, including profile analysis using clustering and kurtosis (PACK) and cancer outlier profile analysis (COPA). Our method is robust, powerful, invariant to shifting and scaling, has no blind spots and has a sample-size-free interpretation.

Following is an example demonstrating on the pattern of measurements SIBER is good at detecting (Figure credit to PubMed Central):


![SIBER examples](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582265/bin/bts713f4p.jpg)

## Installation

The [**devtools** package](http://cran.r-project.org/web/packages/devtools/index.html) is used to install R packages hosted on Github. To install **SIBER**, type the following commands in the R console:

```r
    library(devtools)
    install_github("nickytong/SIBER")
```
## Cite SIBER
Tong, P., Chen, Y., Su, X., & Coombes, K. R. (2013). SIBER: systematic identification of bimodally expressed genes using RNAseq data. Bioinformatics, 29(5), 605-613. url: https://academic.oup.com/bioinformatics/article/29/5/605/247759 
