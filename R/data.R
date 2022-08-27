#' Meta-analysis about meat consumption
#'
#' Meta-analysis of the effectiveness of educational behavior interventions that
#' attempt to reduce meat consumption by appealing to animal welfare.
#'
#' @format A data frame with 100 rows and 4 columns:
#' \describe{
#'   \item{yi}{Point estimate} # TODO MM: on what scale
#'   \item{vi}{Variance of point estimate}
#'   \item{cluster}{Paper from which estimate comes}
#'   \item{randomized}{Logical indicating whether study was randomized}
#' }
#'
#' @references
#' \insertRef{mathur2021}{mbma}
"meta_meat"
