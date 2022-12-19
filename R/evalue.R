#' E-value for meta-analysis with multiple biases
#'
#' @inheritParams metabias::params
#' @inheritParams multibias_meta
#' @param bias_max The largest value of `bias`, on the additive scale, that
#'   should be included in the grid search. The bias has the same units as `yi`.
#' @param assumed_bias_type List of biases to consider for computing evalues
#'   (objects of `bias` as returned by [EValue::confounding()],
#'   [EValue::selection()], [EValue::misclassification()]) (defaults to NULL,
#'   i.e. agnostic as to the nature of the internal bias). If not NULL, the `yi`
#'   argument must be on the log-RR scale (if `yi` is not already on that scale,
#'   use [EValue::convert_measures()] to make it so).
#'
#' @return An object of class [metabias::metabias()].
#'
#' @details For more on the functions passed as `assumed_bias_type`, see the
#'   `EValue` package multiple-bias vignette:
#'   `vignette("multiple-bias", package = "EValue")`
#'
#' @references
#' \insertRef{mathur2022multibias}{metabias}
#'
#' \insertRef{ding2016}{metabias}
#'
#' \insertRef{smith2019}{metabias}
#'
#' \insertRef{vanderweele2019}{metabias}
#'
#' @export
#' @example inst/examples/multibias_evalue.R
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

  # flip yi if not favor_positive
  if (!favor_positive) yi <- -yi

  # resolve vi and sei
  if (missing(vi) && missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
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
    if (abs(opt$minimum - bias_max) < tolerance && opt$objective > tolerance)
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
  stats <- tibble(bias_est = bias_est, bias_ci = bias_ci,
                  evalue_est = evalue_est, evalue_ci = evalue_ci)

  vals <- list(selection_ratio = selection_ratio,
               q = q,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small,
               bias_max = bias_max)

  metabias::metabias(data = dat, values = vals, stats = stats, fits = list())

}
