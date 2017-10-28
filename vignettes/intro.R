## ------------------------------------------------------------------------
library(clusteredinterference)

## ------------------------------------------------------------------------
data("toy_data") 

## ------------------------------------------------------------------------
causal_fx <- policyFX(
  data = toy_data,
  formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
  alphas = c(.3, .5), 
  k_samps = 1
)

causal_fx$estimates

## ---- eval = FALSE-------------------------------------------------------
#  outcome | treatment ~ predictors and random intercept | clustering specification

## ---- eval = FALSE-------------------------------------------------------
#  Treatment ~ Age + Distance + (1 | Cluster_ID)

## ---- eval = FALSE-------------------------------------------------------
#  root_options = list(atol = 1e-7)

## ------------------------------------------------------------------------
my_grid <- makeTargetGrid(alphas = (4:10)/20, small_grid = TRUE) 
my_grid

## ------------------------------------------------------------------------
my_grid$estVar <- FALSE

## ------------------------------------------------------------------------
causal_fx2 <- policyFX(
  data = toy_data,
  formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
  # alphas = c(.3, .5), 
  target_grid = my_grid,
  k_samps = 5,
  verbose = FALSE,
  root_options = list(atol=1e-4)
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
  
  

