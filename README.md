# gcimputeR
This package provides a R implemention to **impute missing values** by fitting a **Gaussian copula model**, on incomplete dataset thay may contain continuous, ordinal and binary variables. The user could either fit a full rank Gaussian copula model [1] or a low rank Gaussian copula model [2]. The fitted model also provides latent correlation among all variables regardless of types. Future release will support mini-batching training and online data imputation [3], parallelization, and categorical variables. A [python implementation](https://github.com/udellgroup/gcimpute) is also available.



## Install  and load the package
The easiest way is to install using `install_github` from the package `devtools`:
```{r}
library(devtools)
install_github("udellgroup/gcimputeR")
library(gcimputeR)
```

## Examples
```{r}
library(gcimputeR)
set.seed(410)
# generate and mask 15-dim mixed data 
# 5 continuous, 5 ordinal (1-5) and 5 boolean
var_types = list('cont'=1:5, 'ord'=6:10, 'bin'=11:15) 
X = generate_GC(n = 2000, var_types = var_types)
X_mask = mask_MCAR(X, mask_fraction = 0.4)

# model fitting
fit = impute_GC(X_mask, verbose = TRUE)

# Evaluation: compute the scaled-MAE (SMAE) for each data type
# (scaled by MAE of median imputation) 
err_imp = cal_smae(xhat = fit$Ximp, xobs = X_obs, xtrue = X, reduce = FALSE) 
for (t in names(var_types)){
  err = round(mean(err_imp[var_types[[t]]]), 4)
  print(paste0('SMAE of ', var_types[t], ' : ', err))
}
```
More detailed examples are available under directory **doc**. 

## News

*06-13-2022 (Version 0.1.1)*

Our software was renamed from *mixedgcImp* to *gcimputeR* and went through substantial structural changes. We improved the code quality and the user interface. Moreover, this release achieves significant speed improvement when there are many ordinal variables. We strongly recommend an update if you are an early user!

## Reference
[1] Zhao, Y. and Udell, M. Missing value imputation for mixed data via Gaussian copula, KDD 2020.

[2] Zhao, Y. and Udell, M. Matrix completion with quantified uncertainty through low rank Gaussian copula, NeurIPS 2020.

[3] Zhao, Y., Landgrebe, E., Shekhtman E., and Udell, M. Online missing value imputation and change point detection using the Gaussian Copula, AAAI 2022.
