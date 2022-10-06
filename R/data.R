#' Meta-analysis about meat consumption
#'
#' Meta-analysis of the effectiveness of educational behavior interventions that
#' attempt to reduce meat consumption by appealing to animal welfare.
#'
#' @format A data frame with 100 rows and 4 columns:
#' \describe{
#'   \item{yi}{Point estimate on log-risk ratio scale}
#'   \item{vi}{Variance of point estimate}
#'   \item{cluster}{Paper that contributed the point estimate}
#'   \item{randomized}{Logical indicating whether study was randomized}
#' }
#'
#' @references
#' \insertRef{mathur2021}{multibiasmeta}
"meta_meat"
