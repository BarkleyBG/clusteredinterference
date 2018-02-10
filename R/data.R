#' A "toy" dataset for illustrating the estimator
#'
#' A dataset containing simulated data for outcome, treatment status, cluster
#' membership, and two pretreatment covariates. There are 22 clusters that have
#' 3 individuals each, and 8 clusters that have 4 individuals each.
#'
#' @format A data frame with 98 rows and 5 variables: \describe{
#'   \item{Outcome}{Individual outcome status observed at the end of follow-up.
#'   E.g., infection status} \item{Treatment}{Individual treatment status
#'   observed at the end of follow-up. E.g., vaccination status}
#'   \item{Cluster_ID}{Unique identifier for the cluster of which the individual
#'   is a member} \item{Age}{Individual's age, in years}
#'   \item{Distance}{Individual's distance to river, in kilometers} }
"toy_data"
