RIS Spatio-Temporal Probit
================
Philipp Hunziker
August 12, 2018

This package implements the Spatio-Temporal Autoregressive Probit model introduced in [Franzese, Hays and Cook (2016)](https://www.cambridge.org/core/journals/political-science-research-and-methods/article/spatial-and-spatiotemporalautoregressive-probit-models-of-interdependent-binary-outcomes/1F3614552196A3506BBA42C42C03C061) (FHC). Estimation is based on a Recursive Importance Sampling (RIS) Maximum Simulated Likelihood (MSL) procedure. The code is loosely based on Franzese et al's Matlab replication materials, with some speed-ups (C++) and convenience functions added.

NOTE: This package is alrgely untested - proceed with caution.

Usage
-----

Here's a short usage examples on simulated data similar to Experiment \#1 in FHC. Note that the model is implemented as a stateful R6 class.

``` r
library(risprobit)
set.seed(0)

## Set data specs
N <- 7*7
TT <- 5  # TT = 21 in FHC, but that takes forever

## Set true params
rho <- 0.1
gamma <- 0.3
beta <- c(-1.5, 3)

## Simulate some data
sim <- simulate_data(N, TT, rho, gamma, beta, X_params = c(-1, 2))
X <- sim$X
y <- sim$y
W_t <- sim$W_t

## Make model object, train
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT)
ris$train()
```

    ## [1] 59.0431

``` r
ris$coeftest()
```

    ## 
    ## z test of coefficients:
    ## 
    ##          Estimate Std. Error z value  Pr(>|z|)    
    ## beta_0 -1.0310565  0.1515921 -6.8015 1.035e-11 ***
    ## beta_1  2.2167956  0.2747631  8.0680 7.144e-16 ***
    ## rho     0.2359175  0.0043966 53.6586 < 2.2e-16 ***
    ## gamma   0.1845312  0.0034784 53.0507 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Notes on Identification and Estimation
--------------------------------------

I noticed some odd behavior in cases where the regressors aren't very predictive of the outcome, and thus identification of the dependence parameters has to occur solely via the errors. Check out the following example:

``` r
set.seed(0)

## Set true params
beta <- c(0, 0)  # X has no predictive value in true model

## Simulate some data
sim <- simulate_data(N, TT, rho, gamma, beta, X_params = c(-1, 2))
X <- sim$X
y <- sim$y
W_t <- sim$W_t

## Make model object, train
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT)
ris$train()
```

    ## [1] 135.5475

``` r
ris$coeftest()
```

    ## 
    ## z test of coefficients:
    ## 
    ##          Estimate Std. Error z value Pr(>|z|)
    ## beta_0  0.0835424  0.1083353  0.7711   0.4406
    ## beta_1 -0.0373322  0.1076553 -0.3468   0.7288
    ## rho    -0.0070533  0.1824547 -0.0387   0.9692
    ## gamma   0.0085199  0.1615840  0.0527   0.9579

I cannot rule out completely that the problem lies with my implementation, but after playing around for a while, I'm fairly certain the RIS-MSL estimation strategy is the culprit.

Another problem is that the RIS-based log-likelihood function does not seem to be well-behaved. Specifically, I frequently get varying results depending on the optimization algorithm employed (which shouldn't happen if the loss function was smooth and convex). Check out the following:

``` r
set.seed(0)

## Set data specs
N <- 7*7
TT <- 2

## Set true params
rho <- 0.1
gamma <- 0.3
beta <- c(0, 0)

## Simulate some data
sim <- simulate_data(N, TT, rho, gamma, beta, X_params = c(-1, 2))
X <- sim$X
y <- sim$y
W_t <- sim$W_t

## Train with optim Nelder-Mead (default)
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT)
cat('\nNelder-Mead:', fill = TRUE)
```

    ## 
    ## Nelder-Mead:

``` r
ll <- ris$train()
cat('ll:', round(ll, 3), fill = TRUE)
```

    ## ll: 31.585

``` r
ris$coeftest()
```

    ## 
    ## z test of coefficients:
    ## 
    ##          Estimate Std. Error z value Pr(>|z|)    
    ## beta_0 -0.3197200  0.1949379 -1.6401   0.1010    
    ## beta_1  0.0835082  0.2707839  0.3084   0.7578    
    ## rho     0.0097419  0.1096338  0.0889   0.9292    
    ## gamma   0.9399844  0.1096346  8.5738   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## Train with optim BFGS
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT)
cat('\nBFGS:', fill = TRUE)
```

    ## 
    ## BFGS:

``` r
ll <- ris$train(method = 'BFGS')
cat('ll:', round(ll, 3), fill = TRUE)
```

    ## ll: 31.609

``` r
ris$coeftest()
```

    ## 
    ## z test of coefficients:
    ## 
    ##         Estimate Std. Error z value Pr(>|z|)
    ## beta_0 -1.230994   1.808819 -0.6806   0.4962
    ## beta_1 -0.041990   0.195419 -0.2149   0.8299
    ## rho     0.013139   0.118592  0.1108   0.9118
    ## gamma  -0.625238   0.402005 -1.5553   0.1199

``` r
## Train with nloptr newuoa
ris <- RisSpatialProbit$new(X = X, y = y, W_t = W_t, N = N, TT = TT)
cat('\nNEWUOA:', fill = TRUE)
```

    ## 
    ## NEWUOA:

``` r
ll <- ris$train(method = 'newuoa')
```

    ## Loading required package: nloptr

``` r
cat('ll:', round(ll, 3), fill = TRUE)
```

    ## ll: 31.63

``` r
ris$coeftest()
```

    ## 
    ## z test of coefficients:
    ## 
    ##          Estimate Std. Error z value Pr(>|z|)
    ## beta_0 -0.3888009  0.3915575 -0.9930   0.3207
    ## beta_1 -0.0083487  0.2367581 -0.0353   0.9719
    ## rho     0.0033366  0.7147622  0.0047   0.9963
    ## gamma  -0.0023531  0.6158874 -0.0038   0.9970
