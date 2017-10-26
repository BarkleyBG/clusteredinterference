clusteredinterference
=====================

The goal of clusteredinterference is to estimate causal effects from observational study data in the presence of clustered interference.

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
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  alphas = c(.3, .5), 
  k_samps = 1,
  verbose = FALSE
)
```

``` r
knitr::kable(causal_fx$estimates, digits = 3)
```

|  alpha1\_num|  alpha2\_num|  trt|  alpha1|  alpha2| estimand | effect\_type | estVar |    k|  estimate|    var|     se|     LCI|     UCI|
|------------:|------------:|----:|-------:|-------:|:---------|:-------------|:-------|----:|---------:|------:|------:|-------:|-------:|
|            1|           NA|   NA|     0.3|      NA| mu       | outcome      | TRUE   |    1|     0.216|  0.001|  0.030|   0.158|   0.274|
|            2|           NA|   NA|     0.5|      NA| mu       | outcome      | TRUE   |    1|     0.160|  0.001|  0.026|   0.109|   0.210|
|            1|           NA|    0|     0.3|      NA| mu0      | outcome      | TRUE   |    1|     0.249|  0.001|  0.037|   0.177|   0.322|
|            2|           NA|    0|     0.5|      NA| mu0      | outcome      | TRUE   |    1|     0.209|  0.002|  0.045|   0.120|   0.298|
|            1|           NA|    1|     0.3|      NA| mu1      | outcome      | TRUE   |    1|     0.099|  0.001|  0.032|   0.036|   0.162|
|            2|           NA|    1|     0.5|      NA| mu1      | outcome      | TRUE   |    1|     0.087|  0.001|  0.023|   0.042|   0.133|
|            1|            1|   NA|     0.3|     0.3| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|   NA|     0.5|     0.3| CE       | contrast     | TRUE   |    1|    -0.056|  0.001|  0.025|  -0.106|  -0.007|
|            1|            2|   NA|     0.3|     0.5| CE       | contrast     | TRUE   |    1|     0.056|  0.001|  0.025|   0.007|   0.106|
|            2|            2|   NA|     0.5|     0.5| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            1|            1|    0|     0.3|     0.3| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|    0|     0.5|     0.3| QE0      | contrast     | TRUE   |    1|    -0.040|  0.001|  0.037|  -0.113|   0.032|
|            1|            2|    0|     0.3|     0.5| QE0      | contrast     | TRUE   |    1|     0.040|  0.001|  0.037|  -0.032|   0.113|
|            2|            2|    0|     0.5|     0.5| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            1|            1|    1|     0.3|     0.3| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|    1|     0.5|     0.3| QE1      | contrast     | TRUE   |    1|    -0.012|  0.001|  0.026|  -0.062|   0.039|
|            1|            2|    1|     0.3|     0.5| QE1      | contrast     | TRUE   |    1|     0.012|  0.001|  0.026|  -0.039|   0.062|
|            2|            2|    1|     0.5|     0.5| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|

Vignette
--------

The vignette provides more information on the formal arguments.

Acknowledgments
---------------

-   Please see the [inferference](https://cran.r-project.org/package=inferference) package for related estimators from the article [Tchetgen Tchetgen and VanderWeele (2012) article, "On causal inference in the presence of interference"](https://doi.org/10.1177/0962280210386779)
-   An earlier version of the methods implemented in `clusteredinterference` was implemented using the [`geex`](https://github.com/bsaul/geex) package for estimating equations.
-   Thanks to [Bradley Saul](https://github.com/bsaul) for `inferference`, `geex`, and for comments and suggestions that were helpful in the creation of `clusteredinterference`.
