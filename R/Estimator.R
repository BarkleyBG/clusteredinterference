
#' Estimate Causal Effects of Population Treatment Policies Assuming Clustered
#' Interference
#'
#' This function implements the estimators from Barkley et al. (2017) for
#' estimating causal effects \emph{("FX")} of treatment policies from an
#' observational study when \strong{clustered interference} is assumed.
#' Clustered interference is also often known as "partial" interference. For the
#' manuscript introducing the methods in \pkg{clusteredinterference}, see: URL
#' \url{https://arxiv.org/abs/1711.04834}
#'
#'
#' @param data A \code{data.frame} (not a \code{tibble}). Columns of
#'   \code{factor} types are not recommended and will sometimes throw
#'   (defensive) errors.
#' @param formula The \code{formula} defines the different components of the
#'   method. The components are specified by \code{outcome | treatment ~
#'   f(covariates) + (1|cluster_id) | cluster_id}. The middle component is
#'   passed to \code{\link[lme4]{glmer}}, so \code{treatment ~ f(covariates) +
#'   (1|cluster_id)} specifies the model form for the propensity score (i.e.,
#'   treatment) model. See Details.
#' @param alphas A numeric vector for the probabilities corresponding to the
#'   policies of interest. Each entry must be between 0 and 1.
#' @param k_samps The maximum number of vectors to evaluate to estimate the
#'   counterfactual probabilities (i.e., \eqn{\omega(A,N,\alpha)}). Setting to 0
#'   avoids approximation at the cost of increased computation time. Recommended
#'   to set <= 5.
#' @param ... The dots argument. The user may supply their own
#'   \code{target_grid} through the dots argument. The \code{target_grid} can be
#'   made through exported function \code{\link{makeTargetGrid}}
#' @param nAGQ This is the number of Adaptive Gaussian Quadrature points used in
#'   the \code{\link[lme4]{glmer}} model fitting computation. Defaults to 2. It
#'   is recommended to use more than 1.
#' @param root_options These are passed to \code{\link[rootSolve]{multiroot}}
#'   function.
#' @param return_matrices A Boolean on whether to return the "bread" and "meat"
#'   matrices in the sandwich variance. Defaults to \code{FALSE}.
#' @param verbose A Boolean on whether to print output to \code{stderr}.
#'   Defaults to FALSE.
#'
#' @details These estimators are based on inverse probability-weighting by the
#'   propensity score for treatment (IPW) to estimate causal effects of
#'   counterfactual policies of interest (i.e., \emph{"policy effects"}) when
#'   clustered interference is assumed. The policies of interest correspond to
#'   counterfactual scenarios in which treatment may be correlated within
#'   clusters.
#'
#'   This method estimates causal contrasts of these policies by estimating the
#'   counterfactual treatment probabilities; taking the correlation structures
#'   into account requires heavy computational resources, so the user should be
#'   patient.
#'
#'   The modeling formula for the propensity score (i.e., treatment) model is
#'   specified via the \code{formula} formal argument. An example of a model
#'   logit-linear fixed effects would be \code{Y | A ~ X1 + X2 + (1 |
#'   cluster_ID) | cluster_ID}. A similar model that also includes an
#'   interaction term is \code{Y | A ~ X1 + X2 + X1:X2 + (1 | cluster_ID) |
#'   cluster_ID}.
#'
#' @author Brian G. Barkley, \email{BarkleyBG@@unc.edu}.
#'
#'   Some of the plumbing functions for estimating the sandwich variance matrix
#'   were adapted from Bradley Saul's
#'   \href{https://cran.r-project.org/package=geex}{\pkg{geex}} package.
#'
#'   Some of the plumbing functions for the logistic mixed model likelihood were
#'   adapted from Bradley Saul's
#'   \href{https://cran.r-project.org/package=inferference}{\pkg{inferference}}
#'   package.
#'
#' @return A \code{list} object including: \enumerate{ \item \code{estimates}: A
#'   tidy \code{data.frame} with columns \code{estimand}, \code{estimate},
#'   \code{var}, \code{se}, \code{LCI} and \code{UCI} for 95\% CI's, and more
#'   information. \item \code{parameters}: An untidy \code{list} of the point
#'   estimates of all (target and nuisance) parameters. \item
#'   \code{variance_matrices}: When \code{return_matrices} is \code{TRUE} this
#'   is a \code{list} object for the "bread" and "meat" matrices in the sandwich
#'   variance calculations for each estimand. Otherwise, it is a \code{list}
#'   object with length 0. \item \code{propensity_scores}: The estimated
#'   propensity scores for each cluster. \item \code{model}: The treatment model
#'   object. \item \code{formula}: The full formula argument provided, after
#'   coercion to a \code{\link[Formula]{Formula}} object}
#'
#' @seealso Please see the main package vignette at
#'   \code{vignette("estimate-policyFX")}. It describes the necessary arguments,
#'   as well as some extra functionality.
#'
#' @references Barkley, B. G., Hudgens, M. G., Clemens, J. D., Ali, M., and
#'   Emch, M. E. (2017). Causal Inference from Observational Studies with
#'   Clustered Interference. \emph{arXiv preprint arXiv:1711.04834}. (URL:
#'   \url{https://arxiv.org/abs/1711.04834}.)
#'
#'   Bradley C Saul and Michael G Hudgens (2017). A Recipe for
#'   \code{inferference}: Start with Causal Inference. Add Interference. Mix
#'   Well with R. \emph{Journal of Statistical Software} \strong{82}(2), pp.
#'   1-21. doi: <10.18637/jss.v082.i02> (URL:
#'   \url{http://doi.org/10.18637/jss.v082.i02}).
#'   \url{https://cran.r-project.org/package=inferference}.
#'   \url{https://github.com/bsaul/inferference}.
#'
#'   Bradley Saul (2017). \code{geex}: An API for M-Estimation.
#'   \url{https://cran.r-project.org/package=geex}.
#'   \url{https://github.com/bsaul/geex}.
#'
#'
#' @author Brian G. Barkley

#' @examples
#' \dontrun{
#' toy_data <- clusteredinterference::toy_data
#' causal_fx <- policyFX(
#'   data = toy_data,
#'   formula = Outcome | Treatment ~ Age + Distance + (1 | Cluster_ID) | Cluster_ID,
#'   alphas = c(.3, .5),
#'   k_samps = 1,
#'   verbose = FALSE
#' )}
#' @export
policyFX <- function(
  data,
  formula,
  alphas,
  k_samps=NULL,
  ...,
  verbose = FALSE,
  root_options = NULL,
  nAGQ=2,
  return_matrices = FALSE
){

  if (is.null(k_samps)) {
    stop("Please specify a numeric for `k_samps` argument. 0 is for full sample, but <=5 is recommended")
  }
  max_k_evals <- k_samps

  dots <- list(...)
  if ("target_grid" %in% names(dots)) {
    target_grid <- dots$target_grid
    alphas_in_grid <- sort(unique(target_grid$alpha1))
    if (missing("alphas")){
      alphas <- alphas_in_grid
    }
    if (!all(alphas == alphas_in_grid)){
      warning("target_grid has different alphas than those specified; the alphas from the target_grid argument will be used")
      alphas <- alphas_in_grid
    }
  } else {
    target_grid <- makeTargetGrid(alphas = alphas)
  }

  target_grid$k <- max_k_evals

  formule <- Formula::Formula(formula)

  arguments <- argTests(
    data = data,
    formule = formule,
    alphas = alphas,
    max_k_evals = max_k_evals
  )

  data <- arguments$data ##sorted
  alphas <- arguments$alphas ##sorted
  modeling_formula <- arguments$modeling_formula
  grouping_var_name<- arguments$grouping_var_name

  ## verbose
  if (verbose){print("I. Estimate GLMER model")}
  ## verbose

  glmer_fit <- lme4::glmer(
    data = data,
    formula = modeling_formula,
    family = stats::binomial,
    nAGQ = nAGQ
  )

  model_parms <- unlist(lme4::getME(glmer_fit, c("beta","theta")))


  split_data_list <- splitData(
    data = data,
    formule = formule,
    glmer_fit = glmer_fit,
    grouping_var_name = grouping_var_name
  )

  ## verbose
  if (verbose){print("IB. get MCFP vecs")}
  ## verbose

  vec_dfm_no_alphas <- sampleVecs4Trt(
    split_data_list = split_data_list,
    max_k_evals = max_k_evals
  )


  index_references <- makeIndexReferences(
    model_parms = model_parms,vec_dfm_no_alphas = vec_dfm_no_alphas)

  typical_args <- list(
    split_data_list = split_data_list,
    model_parms = model_parms,
    verbose = verbose
  )
  CFBI_args <- append(typical_args, list(
    alphas = alphas,
    glmer_fit = glmer_fit,
    root_options = root_options,
    vec_dfm_no_alphas = vec_dfm_no_alphas
  ))
  CFBI_MCFP_list <- do.call(what=estimateNuisance,args=CFBI_args)

  CFBI_estimates <-  CFBI_MCFP_list$CFBI_estimates
  MCFP_tidy_estimates <- CFBI_MCFP_list$MCFP_tidy_estimates


  target_args <- append(typical_args, list(
    CFBI_estimates = CFBI_estimates,
    MCFP_tidy_estimates = MCFP_tidy_estimates,
    alphas = alphas,
    target_grid = target_grid
  ))
  ## to output
  targets_and_cps <- do.call(what=estimateTargets,args=target_args)


  target_estimates <- targets_and_cps$target_grid

  var_args <- target_args
  var_args$target_grid <- NULL
  var_args$target_grid_estimates <- target_estimates


  ## Diagnostics
  MCFP_tidier_estimates <- tidyOutMCFP(MCFP_tidy_estimates)

  ## to output
  parameters <-  list(
    model_parms = model_parms,
    CFBI_estimates = CFBI_estimates,
    MCFP_estimates = MCFP_tidier_estimates,
    target_estimates = target_estimates
  )
  if (verbose){print('Beginning variance calculations. This may take some time.')}

  var_args$index_references <- index_references
  var_args$return_matrices <- return_matrices



  var_ests_list <- do.call(estimateVariance, args = var_args)

  ## to ouptut
  tidy_estimates <- tidyEstimates(var_ests_list$tidy_grid)
  ## to output
  if (return_matrices== TRUE){
    variance_matrices <-  var_ests_list$var_matrices
  } else{
    variance_matrices <- list()
  }



  out_list <- list(
    estimates = tidy_estimates,
    parameters = parameters,
    variance_matrices = variance_matrices,
    propensity_scores = targets_and_cps$prop_scores,
    model = glmer_fit,
    formula = formule
  )
  class(out_list) <- "policyFX"
  out_list

}
