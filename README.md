About the `clusteredinterference` package
=========================================

This package implements the estimators proposed in [Barkley et al. (2017), *Causal Inference from Observational Studies with Clustered Interference*.](https://arxiv.org/abs/1711.04834)

What is clustered interference?
-------------------------------

### What is interference?

In causal inference, when one individual's treatment may affect another individual's outcome, it's often called **interference**. In most applications, it is assumed that there is no interference whatsoever. In some applications this must be relaxed; i.e., as in infectious disease research.

### Clustered interference is partial interference

A relaxation of the assumption of "no interference" is to assume that individuals may be partitioned into distinct clusters of individuals (e.g., households, or classrooms, etc.) such that there may be interference within the clusters, but not between the clusters. Historically, this assumption has been referred to as **partial interference** after Sobel (2006).

We introduce the terminology **clustered interference** to refer to this same assumption. To us, it is a descriptive moniker for the assumption, perhaps clarifying the restriction of interference to clusters.

About the method
----------------

[Barkley et al. (2017)](https://arxiv.org/abs/1711.04834) propose new causal estimands for defining treatment effects in the context of observational studies when there may be interference or spillover effects between units in the same cluster. The manuscript also introduces IPTW estimators for thos estimands, which are implemented in `clusteredinterference`.

The manuscript
--------------

A version of this manuscript is available [on arXiv at 1711.04834](https://arxiv.org/abs/1711.04834).

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

| Platform                                         | Latest succesful check |
|--------------------------------------------------|------------------------|
| macOS 10.11 El Capitan, R-release (experimental) | v0.3.0                 |
| Windows Server 2008 R2 SP1, R-release, 32/64 bit | v0.3.0                 |

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
set.seed(1113)
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
| mu(0.15)       |     0.699|  0.008|  0.089|   0.523|   0.874|    0.15|      NA|   NA| mu             | outcome      |         1|
| mu(0.25)       |     0.666|  0.005|  0.070|   0.529|   0.804|    0.25|      NA|   NA| mu             | outcome      |         1|
| mu0(0.15)      |     0.716|  0.008|  0.092|   0.536|   0.895|    0.15|      NA|    0| mu0            | outcome      |         1|
| mu0(0.25)      |     0.687|  0.006|  0.078|   0.535|   0.839|    0.25|      NA|    0| mu0            | outcome      |         1|
| mu1(0.15)      |     0.162|  0.002|  0.043|   0.078|   0.246|    0.15|      NA|    1| mu1            | outcome      |         1|
| mu1(0.25)      |     0.244|  0.003|  0.054|   0.139|   0.349|    0.25|      NA|    1| mu1            | outcome      |         1|
| OE(0.15,0.15)  |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|   NA| OE             | contrast     |         1|
| OE(0.25,0.15)  |    -0.032|  0.001|  0.028|  -0.086|   0.022|    0.25|    0.15|   NA| OE             | contrast     |         1|
| OE(0.15,0.25)  |     0.032|  0.001|  0.028|  -0.022|   0.086|    0.15|    0.25|   NA| OE             | contrast     |         1|
| OE(0.25,0.25)  |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|   NA| OE             | contrast     |         1|
| SE0(0.15,0.15) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|    0| SE0            | contrast     |         1|
| SE0(0.25,0.15) |    -0.029|  0.001|  0.028|  -0.084|   0.027|    0.25|    0.15|    0| SE0            | contrast     |         1|
| SE0(0.15,0.25) |     0.029|  0.001|  0.028|  -0.027|   0.084|    0.15|    0.25|    0| SE0            | contrast     |         1|
| SE0(0.25,0.25) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|    0| SE0            | contrast     |         1|
| SE1(0.15,0.15) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.15|    0.15|    1| SE1            | contrast     |         1|
| SE1(0.25,0.15) |     0.082|  0.000|  0.013|   0.057|   0.107|    0.25|    0.15|    1| SE1            | contrast     |         1|
| SE1(0.15,0.25) |    -0.082|  0.000|  0.013|  -0.107|  -0.057|    0.15|    0.25|    1| SE1            | contrast     |         1|
| SE1(0.25,0.25) |     0.000|  0.000|  0.000|   0.000|   0.000|    0.25|    0.25|    1| SE1            | contrast     |         1|

Vignette
--------

The vignette provides more information on the formal arguments:

``` r
browseVignettes("clusteredinterference")
```

News and version history
------------------------

A changelog is found in the `NEWS.md` file. Version history is also tracked by the [release tags](https://github.com/BarkleyBG/clusteredinterference/releases) for this GitHub repo.

References and acknowledgments
------------------------------

-   The terminology of **partial interference** is attributed to Sobel (2006):
    -   Sobel, M. E. (2006). *What do randomized studies of housing mobility demonstrate? Causal inference in the face of interference.* Journal of the American Statistical Association, 101(476), 1398-1407. [doi: 10.1198/016214506000000636](http://dx.doi.org/10.1198/016214506000000636)
-   Please see the [`inferference`](https://cran.r-project.org/package=inferference) package for related estimators from the following articles:
    -   Perez‐Heydrich, C., Hudgens, M. G., Halloran, M. E., Clemens, J. D., Ali, M., & Emch, M. E. (2014). Assessing effects of cholera vaccination in the presence of interference. Biometrics, 70(3), 731-741. [doi: 10.1111/biom.12184](doi.wiley.com/10.1111/biom.12184)
    -   Tchetgen, E. J. T., & VanderWeele, T. J. (2012). On causal inference in the presence of interference. Statistical methods in medical research, 21(1), 55-75. [doi: 10.1177/0962280210386779](https://doi.org/10.1177/0962280210386779)
-   An earlier version of the methods implemented in `clusteredinterference` was implemented using the [`geex`](https://github.com/bsaul/geex) package for estimating equations.
-   Thanks to [Bradley Saul](https://github.com/bsaul) for [`inferference`](https://cran.r-project.org/package=inferference), [`geex`](https://github.com/bsaul/geex), and for comments and suggestions that were helpful in the creation of `clusteredinterference`.
