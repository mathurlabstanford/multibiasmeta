#' E-value for meta-analysis with multiple biases
#'
#' @inheritParams multibias_meta
#' @param q The attenuated value to which to shift the point estimate or CI.
#'   Should be specified on the same scale as \code{dat$yi} (e.g., if
#'   \code{dat$yi} is on the log-RR scale, then \code{q} should be as well).
#' @param bias_max The largest value of \code{bias}, on the additive scale,
#'   that should be included in the grid search. The bias has the same units as
#'   \code{yi}.
#' @param assumed_bias_type List of biases to consider for computing evalues
#'   (objects of \code{bias} as returned by \code{EValue::confounding()},
#'   \code{EValue::selection()}, \code{EValue::misclassification()}) (defaults
#'   to NULL, i.e. agnostic as to the nature of the internal bias). If not NULL,
#'   the \code{yi} argument must be on the log-RR scale (if \code{yi} is not
#'   already on that scale, use \code{EValue::convert_measures()} to make it
#'   so).
#'
#' @return
#'
#' @details For more on the functions passed as \code{assumed_bias_type}, see the
#'   EValue package vignette on multiple biases:
#'   <https://cran.r-project.org/web/packages/EValue/vignettes/multiple-bias.html>.
#'
#' @references
#' \insertRef{mathur2022}{multibiasmeta}
#'
#' \insertRef{ding2016}{multibiasmeta}
#'
#' \insertRef{smith2019}{multibiasmeta}
#'
#' \insertRef{vanderweele2019}{multibiasmeta}
#'
#' @export
#' @examples
#' multibias_evalue(yi = meta_meat$yi,
#'                  vi = meta_meat$vi,
#'                  selection_ratio = 4,
#'                  biased = !meta_meat$randomized)
#'
#' # specify confounding as internal bias
#' multibias_evalue(yi = meta_meat$yi,
#'                  vi = meta_meat$vi,
#'                  selection_ratio = 4,
#'                  biased = !meta_meat$randomized,
#'                  assumed_bias_type = list(EValue::confounding()))
multibias_evalue <- function(yi, # data
                             vi,
                             sei,
                             cluster = 1:length(yi),
                             biased = TRUE,

                             selection_ratio, # params
                             q = 0,

                             favor_positive = TRUE, # opts
                             alpha_select = 0.05,
                             ci_level = 0.95,
                             small = TRUE,
                             bias_max = 20,
                             assumed_bias_type = NULL) {

  # warn if naive estimate is in opposite direction than favor_positive
  naive_pos <- metafor::rma(yi, vi, method = "FE")$beta > 0
  if (naive_pos != favor_positive)
    warning("Favored direction is opposite of the pooled estimate.")

  # set up multibias_meta with either vi or sei passed
  # if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  # if (missing(sei)) corrected_fun <- partial(multibias_meta, vi = vi)
  # else corrected_fun <- partial(multibias_meta, sei = sei)

  # resolve vi and sei
  if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(vi)) vi <- sei ^ 2
  if (missing(sei)) sei <- sqrt(vi)

  # compute expectation of bias for parameter (mu_hat or ci_lower)
  compute_eb <- function(param, tolerance = 0.0001) {

    # function to optimize: multibias_meta for a given bias
    bias_factor <- function(bias_shared) {
      # all arguments passed through except biases being set to bias_shared
      corrected <- suppressWarnings(
        multibias_meta(yi = yi,
                       vi = vi,
                       bias_affirmative = bias_shared,
                       bias_nonaffirmative = bias_shared,
                       selection_ratio = selection_ratio,
                       cluster = cluster,
                       biased = biased,
                       favor_positive = favor_positive,
                       alpha_select = alpha_select,
                       ci_level = ci_level,
                       small = small))
      # difference between parameter estimate and q
      return(abs(corrected$stats[[param]] - q))
    }

    # optimize to find smallest difference
    opt <- optimize(f = bias_factor, interval = c(0, bias_max))
    # check that search stayed within bounds
    if (abs(opt$minimum - bias_max) < tolerance & opt$objective > tolerance)
      return(paste(">", bias_max))
    return(opt$minimum)
  }

  # find biases for mu_hat and ci_lower
  bias_est <- compute_eb("mu_hat")
  bias_ci <- compute_eb("ci_lower")

  if (is.null(assumed_bias_type)) {
    evalue_est <- evalue_ci <- NULL
  } else {
    # turn list of internal biases into EValue::multi_bias object
    biases <- exec(EValue::multi_bias, !!!assumed_bias_type)
    evalue_est <- summary(EValue::multi_evalue(est = EValue::RR(exp(bias_est)),
                                               biases = biases))
    # convert biases to evalues
    evalue_ci <- summary(EValue::multi_evalue(est = EValue::RR(exp(bias_ci)),
                                              biases = biases))
  }

  dat <- tibble(yi, vi, sei, cluster, biased)
  stats <- list(bias_est = bias_est, bias_ci = bias_ci,
                evalue_est = evalue_est, evalue_ci = evalue_ci)

  vals <- list(selection_ratio = selection_ratio,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small,
               q = q,
               bias_max = bias_max)

  results <- list(data = dat, values = vals, stats = stats, fit = list())
  class(results) <- "metabias"
  return(results)

}
