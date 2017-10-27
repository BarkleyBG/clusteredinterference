About the `clusteredinterference` package
=========================================

This package implements the estimator proposed in Barkley et al. (201X), *Causal Inference from Observational Studies with Clustered Interference*.

The manuscript is in submission.

About the method
----------------

Barkley et al. (201X) propose new causal estimands for defining treatment effects in the context of observational studies when there may be interference or spillover effects between units in the same cluster. The manuscript also introduces IPTW estimators for thos estimands, which are implemented in `clusteredinterference`.

Using the `clusteredinterference` package
=========================================

Install the package
-------------------

You'll need to use the `devtools` package to install from GitHub:

``` r
devtools::install_github("BarkleyBG/clusteredinterference")
```

Load this package
-----------------

``` r
library(clusteredinterference)
```

A quick data example
--------------------

``` r
library(inferference)
data("vaccinesim", package = "inferference") 
head(vaccinesim)
#>   Y        X1       X2 A B group
#> 1 1 5.3607405 1.715527 0 0     1
#> 2 0 0.1964597 1.730802 0 1     1
#> 3 0 0.4846243 1.769546 1 1     1
#> 4 0 0.8012977 1.715527 0 1     1
#> 5 0 2.1426629 1.772158 1 1     1
#> 6 0 1.2861017 1.715527 0 1     1
```

Estimation
----------

Estimation is implemented with one function:

``` r
causal_fx <- policyFX(
  data = vaccinesim[1:200,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  alphas = c(.35, .40), 
  k_samps = 1
)
```

``` r
knitr::kable(causal_fx$estimates, digits = 3)
```

|  alpha1\_num|  alpha2\_num|  trt|  alpha1|  alpha2| estimand | effect\_type | estVar |    k|  estimate|    var|     se|     LCI|    UCI|
|------------:|------------:|----:|-------:|-------:|:---------|:-------------|:-------|----:|---------:|------:|------:|-------:|------:|
|            1|           NA|   NA|    0.35|      NA| mu       | outcome      | TRUE   |    1|     0.203|  0.001|  0.032|   0.140|  0.265|
|            2|           NA|   NA|    0.40|      NA| mu       | outcome      | TRUE   |    1|     0.195|  0.001|  0.033|   0.130|  0.261|
|            1|           NA|    0|    0.35|      NA| mu0      | outcome      | TRUE   |    1|     0.243|  0.002|  0.040|   0.164|  0.321|
|            2|           NA|    0|    0.40|      NA| mu0      | outcome      | TRUE   |    1|     0.250|  0.002|  0.045|   0.162|  0.337|
|            1|           NA|    1|    0.35|      NA| mu1      | outcome      | TRUE   |    1|     0.140|  0.002|  0.039|   0.063|  0.217|
|            2|           NA|    1|    0.40|      NA| mu1      | outcome      | TRUE   |    1|     0.133|  0.001|  0.035|   0.064|  0.202|
|            1|            1|   NA|    0.35|    0.35| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|
|            2|            1|   NA|    0.40|    0.35| CE       | contrast     | TRUE   |    1|    -0.007|  0.000|  0.008|  -0.024|  0.009|
|            1|            2|   NA|    0.35|    0.40| CE       | contrast     | TRUE   |    1|     0.007|  0.000|  0.008|  -0.009|  0.024|
|            2|            2|   NA|    0.40|    0.40| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|
|            1|            1|    0|    0.35|    0.35| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|
|            2|            1|    0|    0.40|    0.35| QE0      | contrast     | TRUE   |    1|     0.007|  0.000|  0.010|  -0.013|  0.027|
|            1|            2|    0|    0.35|    0.40| QE0      | contrast     | TRUE   |    1|    -0.007|  0.000|  0.010|  -0.027|  0.013|
|            2|            2|    0|    0.40|    0.40| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|
|            1|            1|    1|    0.35|    0.35| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|
|            2|            1|    1|    0.40|    0.35| QE1      | contrast     | TRUE   |    1|    -0.008|  0.000|  0.011|  -0.028|  0.013|
|            1|            2|    1|    0.35|    0.40| QE1      | contrast     | TRUE   |    1|     0.008|  0.000|  0.011|  -0.013|  0.028|
|            2|            2|    1|    0.40|    0.40| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|  0.000|

Vignette
--------

The vignette provides more information on the formal arguments.

Acknowledgments
---------------

-   Please see the [`inferference`](https://cran.r-project.org/package=inferference) package for related estimators from the following articles:
    -   Perez-Heydrich, C., Hudgens, M. G., Halloran, M. E., Clemens, J. D., Ali, M. and Emch, M. E. (2014), *Assessing effects of cholera vaccination in the presence of interference*. Biometrics, 70: 731–741. [doi: 10.1111/biom.12184](doi.wiley.com/10.1111/biom.12184)
    -   Eric J Tchetgen Tchetgen, Tyler J VanderWeele (2014), *On causal inference in the presence of interference*. Statistical Methods in Medical Research. Vol 21, Issue 1, pp. 55 - 75 [doi: 10.1177/0962280210386779](https://doi.org/10.1177/0962280210386779)
-   An earlier version of the methods implemented in `clusteredinterference` was implemented using the [`geex`](https://github.com/bsaul/geex) package for estimating equations.
-   Thanks to [Bradley Saul](https://github.com/bsaul) for `inferference`, `geex`, and for comments and suggestions that were helpful in the creation of `clusteredinterference`.
