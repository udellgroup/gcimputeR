## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----simulate observation------------------------------------------------
library(mixedgcImp)
library(MASS)
set.seed(424)
n = 2000
p = 15
Sigma = matrix(rnorm(p*p), nrow = p)
Sigma = cov2cor(Sigma %*% t(Sigma)) # generate random correlation matrix
Z = mvrnorm(n, mu = rep(0,p), Sigma = Sigma) # generate copula observation 
X = Z # generate observation with specified marginals
# five continuous columns with exponential marginal
X[,1:5] = qexp(pnorm(Z[,1:5]), rate = 1/3) 
# five binary columns
for (i in 6:10) X[,i] = continuous2ordinal(Z[,i], k = 2) 
# five ordinal columns with five levels
for (i in 11:15) X[,i] = continuous2ordinal(Z[,i], k = 5) 
 # mask 50% of the original observation
X_obs = X
loc = sample(1:prod(n*p), size = floor(prod(n*p)*0.5))
X_obs[loc] = NA

## ----fit-----------------------------------------------------------------
fit = impute_mixedgc(X_obs)# takes around 20 secs

## ----evaluate imp--------------------------------------------------------
err_imp = cal_mae_scaled(xhat = fit$Ximp, xobs = X_obs, xtrue = X) # SMAE for each column
mean(err_imp[1:5]) # SMAE across continuous columns
mean(err_imp[6:10]) # SMAE across binary columns
mean(err_imp[11:15]) # SMAE across ordinal columns

## ----evaluate sigma------------------------------------------------------
# relative Frobeinus error of imputed correlation matrix
norm(fit$R - Sigma, type = 'F')/norm(fit$R, type = 'F') 

## ----monitor, fig.height = 4, fig.width = 6, fig.align = "center"--------
plot(fit$loglik, ylab = 'value', type = 'b', pch = 20,
     xlab = 'iterations', main = 'approxiamted likelihood values achieved during iterations')

