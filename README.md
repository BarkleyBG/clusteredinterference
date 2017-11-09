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

You'll need to use the `devtools` package to install :

``` r
devtools::install_github("BarkleyBG/clusteredinterference")
```

### Architectures supported

This package has been checked on several operating systems with package.

| --- Platform --- | -- Latest succesful check -- | |macOS 10.11 El Capitan, R-release (experimental)| v0.1 | |Windows Server 2008 R2 SP1, R-release, 32/64 bit | v0.1|

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
#>   Outcome Treatment Cluster_ID      Age Distance
#> 1       1         0          1 37.62826 6.485258
#> 2       0         1          1 36.61508 6.928957
#> 3       1         1          1 31.74776 6.659470
#> 4       1         0          2 34.79259 8.138802
#> 5       0         0          2 48.05607 7.736209
#> 6       0         0          2 42.21215 8.023865
```

Estimation
----------

Estimation is carried out with one function:

``` r
causal_fx <- policyFX(
  data = toy_data,
  formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
  alphas = c(.15, .25), 
  k_samps = 1
)
```

``` r
knitr::kable(causal_fx$estimates, digits = 3)
```

| estimand       |  estimate|    var|     se|     LCI|     UCI|  alpha1|  alpha2|  trt| estimand\_type | effect\_type |  k\_samps|
|:---------------|---------:|------:|------:|-------:|-------:|-------:|-------:|----:|:---------------|:-------------|---------:|
| mu(0.15)       |     0.699|  0.008|  0.091|   0.521|   0.877|    0.15|      NA|   NA| mu             | outcome      |         1|
| mu(0.25)       |     0.667|  0.005|  0.072|   0.526|   0.808|    0.25|      NA|   NA| mu             | outcome      |         1|
| mu0(0.15)      |     0.716|  0.009|  0.093|   0.533|   0.899|    0.15|      NA|    0| mu0            | outcome      |         1|
| mu0(0.25)      |     0.688|  0.006|  0.079|   0.532|   0.843|    0.25|      NA|    0| mu0            | outcome      |         1|
| mu1(0.15)      |     0.162|  0.002|  0.043|   0.078|   0.246|    0.15|      NA|    1| mu1            | outcome      |         1|
| mu1(0.25)      |     0.244|  0.003|  0.054|   0.139|   0.349|    0.25|      NA|    1| mu1            | outcome      |         1|
| CE(0.15,0.15)  |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|   NA| CE             | contrast     |         1|
| CE(0.25,0.15)  |    -0.032|  0.001|  0.028|  -0.086|   0.022|    0.25|    0.15|   NA| CE             | contrast     |         1|
| CE(0.15,0.25)  |     0.032|  0.001|  0.028|  -0.022|   0.086|    0.15|    0.25|   NA| CE             | contrast     |         1|
| CE(0.25,0.25)  |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|   NA| CE             | contrast     |         1|
| QE0(0.15,0.15) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|    0| QE0            | contrast     |         1|
| QE0(0.25,0.15) |    -0.029|  0.001|  0.028|  -0.084|   0.026|    0.25|    0.15|    0| QE0            | contrast     |         1|
| QE0(0.15,0.25) |     0.029|  0.001|  0.028|  -0.026|   0.084|    0.15|    0.25|    0| QE0            | contrast     |         1|
| QE0(0.25,0.25) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|    0| QE0            | contrast     |         1|
| QE1(0.15,0.15) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|    1| QE1            | contrast     |         1|
| QE1(0.25,0.15) |     0.082|  0.000|  0.013|   0.057|   0.107|    0.25|    0.15|    1| QE1            | contrast     |         1|
| QE1(0.15,0.25) |    -0.082|  0.000|  0.013|  -0.107|  -0.057|    0.15|    0.25|    1| QE1            | contrast     |         1|
| QE1(0.25,0.25) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|    1| QE1            | contrast     |         1|

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
