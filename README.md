# transCCA: A package for estimation transelliptical CCA using the transformed Kendall’s scatter matrix estimator

The `transCCA` package estimates the transelliptical CCA directions and
correlations, as well as their asymptotic variances using the
transformed Kendall’s scatter matrix estimator.

## Installation

``` install
devtools::install_github("blangworthy/transCCA")
library(transCCA)
```

\#\#Examples

``` r
library(transCCA)
library(mvtnorm)

sigmax <- diag(8)
sigmay <- diag(8)
sigmaxy <- diag(c(0.9,0.5,0.4,1/3,rep(0,4)))

sigma <- rbind(cbind(sigmax,sigmaxy),cbind(sigmaxy,sigmay))

z <- mvtnorm::rmvt(n=200,sigma = sigma,df = 5)

transccout <- transCCA(z[,1:8],z[,9:16])
```
