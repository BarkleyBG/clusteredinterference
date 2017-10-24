## ------------------------------------------------------------------------
library(clusteredinterference)

## ------------------------------------------------------------------------
library(inferference)
data("vaccinesim", package = "inferference") 
head(vaccinesim)


## ------------------------------------------------------------------------
type_C <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  alphas = c(.3, .5), 
  k_samps = 1,
  verbose = FALSE
)


## ---- eval = FALSE-------------------------------------------------------
#  Y | A ~ X1 + X2 + (1|group) | group
#  outcome | treatment ~ predictors and random intercept | clustering specification

## ---- eval = FALSE-------------------------------------------------------
#  A ~ X1 + X2 + (1|group)

## ---- eval = FALSE-------------------------------------------------------
#  root_options = list(atol = 1e-7)

## ------------------------------------------------------------------------
my_grid <- makeTargetGrid(alphas = c(0.3, 0.5), small_grid = TRUE)
my_grid <- my_grid[is.na(my_grid$trt), ]
my_grid

## ------------------------------------------------------------------------
type_C <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  alphas = c(.3, .5), 
  k_samps = 1,
  verbose = FALSE
)


## ---- eval = FALSE-------------------------------------------------------
#  Y | A ~ X1 + X2 + (1|group) | group
#  outcome | treatment ~ predictors and random intercept | clustering specification

## ---- eval = FALSE-------------------------------------------------------
#  A ~ X1 + X2 + (1|group)

## ---- eval = FALSE-------------------------------------------------------
#  root_options = list(atol = 1e-7)

## ------------------------------------------------------------------------
my_grid <- makeTargetGrid(alphas = c(0.3, 0.5), small_grid = TRUE)
my_grid <- my_grid[is.na(my_grid$trt), ]
my_grid

## ------------------------------------------------------------------------
type_C <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  # alphas = c(.3, .5), 
  target_grid = my_grid,
  k_samps = 1,
  verbose = TRUE
)


## ------------------------------------------------------------------------
type_C <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  # alphas = c(.3, .5), 
  target_grid = my_grid,
  k_samps = 1,
  verbose = TRUE
)


