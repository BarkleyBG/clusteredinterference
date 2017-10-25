## ------------------------------------------------------------------------
library(clusteredinterference)

## ------------------------------------------------------------------------
library(inferference)
data("vaccinesim", package = "inferference") 
head(vaccinesim)


## ------------------------------------------------------------------------
causal_fx <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  alphas = c(.3, .5), 
  k_samps = 1,
  verbose = FALSE
)

causal_fx$estimates

## ---- eval = FALSE-------------------------------------------------------
#  Y | A ~ X1 + X2 + (1|group) | group
#  outcome | treatment ~ predictors and random intercept | clustering specification

## ---- eval = FALSE-------------------------------------------------------
#  A ~ X1 + X2 + (1|group)

## ---- eval = FALSE-------------------------------------------------------
#  root_options = list(atol = 1e-7)

## ------------------------------------------------------------------------
my_grid <- makeTargetGrid(alphas = (4:10)/20, small_grid = TRUE) 
my_grid

## ------------------------------------------------------------------------
my_grid$estVar <- FALSE

## ------------------------------------------------------------------------
causal_fx2 <- policyFX(
  data = vaccinesim[1:300,],
  formula = Y | A ~ X1 + X2 + (1|group) | group,
  # alphas = c(.3, .5), 
  target_grid = my_grid,
  k_samps = 5,
  verbose = FALSE
)

knitr::kable(causal_fx2$estimates, digits = 3)

## ------------------------------------------------------------------------
library(ggplot2)
library(magrittr)
library(dplyr)

causal_fx2$estimates %>% 
  filter(is.na(alpha2)) %>% ##Only the mu's
  ggplot(aes(
    x = alpha1, 
    y = estimate, 
    group = estimand,
    color = estimand,
    linetype = estimand
    )) + 
  geom_point() +
  geom_line() + 
  geom_hline(yintercept=0)+
  theme_bw() + 
  labs(
    title = "Estimated Population Means"
  ) + 
  theme(legend.position = "bottom") + 
  coord_cartesian(ylim = c(0,0.3))
  
  

