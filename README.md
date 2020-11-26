# mixedgcImp
This R package provides implemention to fit a full rank Gaussian copula model [1] or low rank Gaussian copula model [2], on continuous, ordinal and binary mixed data with missing values. The fitted model can be used for missing value imputation and correlation structure learning [1,2].
Online and mini-batch implementation of Gaussian copula [3] are now available in [Python](https://github.com/udellgroup/online_mixed_gc_imp). R implementation to come...

## install  and load the package
library(devtools)
install_github("udellgroup/mixedgcImp")
library(mixedgcImp)

## Vignette
Some examples can be found under "doc/example_sim.pdf".

## Reference
[1] Zhao, Y. and Udell, M. Missing value imputation for mixeddata via Gaussian copula, KDD 2020.

[2] Zhao, Y. and Udell, M. Matrix completion with quantified uncertainty through low rank Gaussian copula, NeurIPS 2020.

[3] Zhao, Y., Landgrebe, E., Shekhtman E., and Udell, M. Online missing value imputation and correlation change detection for mixed-type Data via Gaussian Copula, arXiv 2020.
