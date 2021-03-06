
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/clusteredinterference)](https://cran.r-project.org/package=clusteredinterference)
[![Travis-CI Build
Status](https://travis-ci.org/BarkleyBG/clusteredinterference.svg?branch=master)](https://travis-ci.org/BarkleyBG/clusteredinterference)
[![Coverage
Status](https://codecov.io/github/BarkleyBG/clusteredinterference/coverage.svg?branch=master)](https://codecov.io/github/BarkleyBG/clusteredinterference)
<!-- badges: end -->

# About the ‘clusteredinterference’ package

This package implements the estimators proposed in [Barkley et al.
(2017), *Causal Inference from Observational Studies with Clustered
Interference*](https://arxiv.org/abs/1711.04834) for estimating the
causal effects of different treatment policies in the presence of
partial or clustered interference. The package is available [on
CRAN](https://CRAN.R-project.org/package=clusteredinterference) with a
[companion website](https://barkleybg.github.io/clusteredinterference/)

## What is clustered interference?

### What is interference?

In causal inference, when one individual’s treatment may affect another
individual’s outcome, it’s often called **interference**. In most
applications, it is assumed that there is no interference whatsoever. In
some applications this must be relaxed - e.g., as in infectious disease
research.

### Clustered interference is partial interference

A relaxation of the assumption of “no interference” is to assume that
individuals may be partitioned into distinct clusters of individuals
(e.g., households, or classrooms, etc.) such that there may be
interference within the clusters, but not between the clusters.
Historically, this assumption has been referred to as **partial
interference** after Sobel (2006).

[Barkley et al. (2017)](https://arxiv.org/abs/1711.04834) introduces the
terminology **clustered interference** to refer to this same assumption.
This phrase may be sufficiently descriptive of the underlying
assumption, and perhaps clarifies the presumed restriction of
interference to clusters.

# About the method

[Barkley et al. (2017)](https://arxiv.org/abs/1711.04834) proposes new
causal estimands for defining treatment effects in the context of
observational studies when there may be interference or spillover
effects between units in the same cluster. The manuscript also
introduces IPTW estimators for thos estimands, which are implemented in
‘clusteredinterference’.

## The manuscript

A version of this manuscript is available [on arXiv
at 1711.04834](https://arxiv.org/abs/1711.04834):

Barkley, B. G., Hudgens, M. G., Clemens, J. D., Ali, M., and Emch, M. E.
(2017). Causal inference from observational studies with clustered
interference. *arXiv preprint arXiv:1711.04834*. URL
<https://arxiv.org/abs/1711.04834>.

# Using the ‘clusteredinterference’ package

## Install the package

This package is now [on
CRAN](https://cran.r-project.org/web/packages/clusteredinterference/index.html)\!

``` r
install.packages("clusteredinterference")
```

Or, visit the [GitHub
repo:](https://github.com/BarkleyBG/clusteredinterference)

``` r
# devtools::install_github("BarkleyBG/clusteredinterference")
```

## Load this package

``` r
library(clusteredinterference)
```

## A quick data example

``` r
data("toy_data")
head(toy_data)
#>   Outcome Treatment Cluster_ID      Age Distance
#> 1       1         0          1 37.62826 6.485258
#> 2       0         1          1 36.61508 6.928957
#> 3       1         1          1 31.74776 6.659470
#> 4       1         0          2 34.79259 8.138802
#> 5       0         0          2 48.05607 7.736209
#> 6       0         0          2 42.21215 8.023865
```

## Example

Estimation is carried out with one function:

``` r
set.seed(1113)
causal_fx <- policyFX(
  data = toy_data,
  formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
  alphas = c(.15, .25), 
  k_samps = 1
)
```

The estimates of causal estimands are printed in a tidy dataframe:

``` r
causal_fx
#> ------------- causal estimates --------------
#>       estimand estimate     se     LCI    UCI
#>       mu(0.15)   0.6985 0.0893  0.5234 0.8736
#>       mu(0.25)   0.6664 0.0702  0.5287 0.8041
#>      mu0(0.15)   0.7157 0.0917  0.5360 0.8954
#>      mu0(0.25)   0.6869 0.0775  0.5350 0.8388
#>      mu1(0.15)   0.1619 0.0429  0.0779 0.2460
#>      mu1(0.25)   0.2440 0.0536  0.1389 0.3491
#>  OE(0.25,0.15)  -0.0321 0.0275 -0.0861 0.0219
#>  OE(0.15,0.25)   0.0321 0.0275 -0.0219 0.0861
#>           ... and 4 more rows ... 
#> ---------------------------------------------
```

Use `summary()` for a little more information:

``` r
summary(causal_fx)
#> ------------- causal estimates --------------
#>       estimand estimate     se     LCI    UCI
#>       mu(0.15)   0.6985 0.0893  0.5234 0.8736
#>       mu(0.25)   0.6664 0.0702  0.5287 0.8041
#>      mu0(0.15)   0.7157 0.0917  0.5360 0.8954
#>      mu0(0.25)   0.6869 0.0775  0.5350 0.8388
#>      mu1(0.15)   0.1619 0.0429  0.0779 0.2460
#>      mu1(0.25)   0.2440 0.0536  0.1389 0.3491
#>  OE(0.25,0.15)  -0.0321 0.0275 -0.0861 0.0219
#>  OE(0.15,0.25)   0.0321 0.0275 -0.0219 0.0861
#> 
#>           ... and 4 more rows ... 
#> 
#> -------------- treatment model -------------
#> Generalized linear mixed model fit by maximum likelihood (Adaptive
#>   Gauss-Hermite Quadrature, nAGQ = 2) [glmerMod]
#>  Family: binomial  ( logit )
#> Formula: Treatment ~ Age + Distance + (1 | Cluster_ID)
#>    Data: data
#>      AIC      BIC   logLik deviance df.resid 
#> 137.0345 147.3743 -64.5172 129.0345       94 
#> Random effects:
#>  Groups     Name        Std.Dev.
#>  Cluster_ID (Intercept) 1.18    
#> Number of obs: 98, groups:  Cluster_ID, 30
#> Fixed Effects:
#> (Intercept)          Age     Distance  
#>    -1.44609     -0.00851      0.26097  
#> 
#> ------------- propensity scores -------------
#>      1      2      3      4      5      6      7      8      9     10 
#>  0.105  0.162  0.086  0.102  0.167  0.045  0.244 0.0934 0.0765  0.197 
#>     11     12     13     14     15     16     17     18     19     20 
#> 0.0653  0.281  0.104  0.365 0.0867  0.198  0.207  0.106 0.0847  0.134 
#>     21     22     23     24     25     26     27     28     29     30 
#>  0.103  0.111  0.105  0.302 0.0434 0.0943 0.0443 0.0512   0.13  0.263 
#> ---------------------------------------------
```

Note that `Treatment ~ Age + Distance + (1 | Cluster_ID)` in the the
middle of the `formula` argument is sent to `lme4::glmer()` to specify
the form of the (logit-link binomial) treatment model.

The `policyFX()` output list includes an element, `formula`, for the
`Formula` object:

``` r
causal_fx$formula
#> Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID
```

The output list also includes an element, `model`, which is the fitted
`glmerMod` S4 model object. Here we can see that the middle of `formula`
was passed into the `glmer()` logit-link binomial mixed model:

``` r
causal_fx$model@call
#> lme4::glmer(formula = Treatment ~ Age + Distance + (1 | Cluster_ID), 
#>     data = data, family = stats::binomial, nAGQ = nAGQ)
```

The fitted model estimates three fixed effects (intercept, a term for
`Age` and a term for `Distance`) and one random effect (for
`Cluster_ID`):

``` r
lme4::getME(causal_fx$model, c("beta", "theta"))
#> $beta
#> [1] -1.446087049 -0.008509771  0.260968952
#> 
#> $theta
#> Cluster_ID.(Intercept) 
#>               1.180325
```

## Vignette

The vignette provides more information on the formal arguments:

``` r
vignette("estimate-policyFX")
```

# News and version history

A changelog is found in the `NEWS.md` file. Version history is also
tracked by the [release
tags](https://github.com/BarkleyBG/clusteredinterference/releases) for
this GitHub repo.

# References and acknowledgments

  - The manuscript introducing the methods in ‘clusteredinterference’
    is:
      - Barkley, B. G., Hudgens, M. G., Clemens, J. D., Ali, M., and
        Emch, M. E. (2017). Causal inference from observational studies
        with clustered interference. *arXiv preprint arXiv:1711.04834*.
        URL <https://arxiv.org/abs/1711.04834>.
  - The terminology of **partial interference** is attributed to Sobel
    (2006):
      - Sobel, M. E. (2006). What do randomized studies of housing
        mobility demonstrate? Causal inference in the face of
        interference. *Journal of the American Statistical Association*,
        101(476), 1398-1407.
        [doi: 10.1198/016214506000000636](http://dx.doi.org/10.1198/016214506000000636)
  - Please see the
    [`inferference`](https://cran.r-project.org/package=inferference)
    package for related estimators from the following articles:
      - Perez‐Heydrich, C., Hudgens, M. G., Halloran, M. E., Clemens, J.
        D., Ali, M., & Emch, M. E. (2014). Assessing effects of cholera
        vaccination in the presence of interference. *Biometrics*,
        70(3), 731-741. [doi:
        10.1111/biom.12184](doi.wiley.com/10.1111/biom.12184)
      - Tchetgen, E. J. T., & VanderWeele, T. J. (2012). On causal
        inference in the presence of interference. *Statistical Methods
        in Medical Research*, 21(1), 55-75.
        [doi: 10.1177/0962280210386779](https://doi.org/10.1177/0962280210386779)  
  - An earlier version of the methods implemented in
    ‘clusteredinterference’ was implemented using the
    [`geex`](https://github.com/bsaul/geex) package for estimating
    equations.
  - Thanks to [Bradley Saul](https://github.com/bsaul) for
    [`inferference`](https://cran.r-project.org/package=inferference),
    [`geex`](https://github.com/bsaul/geex), and for comments and
    suggestions that were helpful in the creation of
    ‘clusteredinterference’.
