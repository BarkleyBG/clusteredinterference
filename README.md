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
data("toy_data")
head(toy_data)
#>        Age Distance Treatment Outcome Cluster_ID
#> 1 37.62826 6.485258         0       1          1
#> 2 36.61508 6.928957         1       0          1
#> 3 31.74776 6.659470         1       1          1
#> 4 34.79259 8.138802         0       1          2
#> 5 48.05607 7.736209         0       0          2
#> 6 42.21215 8.023865         0       0          2
```

Estimation
----------

Estimation is implemented with one function:

``` r
causal_fx <- policyFX(
  data = toy_data,
  formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
  alphas = c(.35, .40), 
  k_samps = 1
)
```

``` r
knitr::kable(causal_fx$estimates, digits = 3)
```

|  alpha1\_num|  alpha2\_num|  trt|  alpha1|  alpha2| estimand | effect\_type | estVar |    k|  estimate|    var|     se|     LCI|     UCI|
|------------:|------------:|----:|-------:|-------:|:---------|:-------------|:-------|----:|---------:|------:|------:|-------:|-------:|
|            1|           NA|   NA|    0.35|      NA| mu       | outcome      | TRUE   |    1|     0.640|  0.004|  0.062|   0.519|   0.761|
|            2|           NA|   NA|    0.40|      NA| mu       | outcome      | TRUE   |    1|     0.628|  0.003|  0.058|   0.514|   0.742|
|            1|           NA|    0|    0.35|      NA| mu0      | outcome      | TRUE   |    1|     0.653|  0.006|  0.074|   0.508|   0.799|
|            2|           NA|    0|    0.40|      NA| mu0      | outcome      | TRUE   |    1|     0.633|  0.005|  0.073|   0.490|   0.777|
|            1|           NA|    1|    0.35|      NA| mu1      | outcome      | TRUE   |    1|     0.313|  0.003|  0.057|   0.201|   0.425|
|            2|           NA|    1|    0.40|      NA| mu1      | outcome      | TRUE   |    1|     0.344|  0.003|  0.058|   0.231|   0.457|
|            1|            1|   NA|    0.35|    0.35| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|   NA|    0.40|    0.35| CE       | contrast     | TRUE   |    1|    -0.012|  0.000|  0.009|  -0.030|   0.007|
|            1|            2|   NA|    0.35|    0.40| CE       | contrast     | TRUE   |    1|     0.012|  0.000|  0.009|  -0.007|   0.030|
|            2|            2|   NA|    0.40|    0.40| CE       | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            1|            1|    0|    0.35|    0.35| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|    0|    0.40|    0.35| QE0      | contrast     | TRUE   |    1|    -0.020|  0.000|  0.010|  -0.040|   0.001|
|            1|            2|    0|    0.35|    0.40| QE0      | contrast     | TRUE   |    1|     0.020|  0.000|  0.010|  -0.001|   0.040|
|            2|            2|    0|    0.40|    0.40| QE0      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            1|            1|    1|    0.35|    0.35| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|
|            2|            1|    1|    0.40|    0.35| QE1      | contrast     | TRUE   |    1|     0.031|  0.000|  0.005|   0.021|   0.040|
|            1|            2|    1|    0.35|    0.40| QE1      | contrast     | TRUE   |    1|    -0.031|  0.000|  0.005|  -0.040|  -0.021|
|            2|            2|    1|    0.40|    0.40| QE1      | contrast     | TRUE   |    1|     0.000|  0.000|  0.000|   0.000|   0.000|

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
