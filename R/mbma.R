#' Correction for meta-analysis with multiple biases
#'
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances (i.e., squared standard errors) for
#'   the point estimates.
#' @param sei A vector of estimated standard errors for the point estimates.
#'   (Only one of \code{vi} or \code{sei} needs to be specified).
#' @param selection_ratio Ratio by which publication bias favors affirmative
#'   studies.
#' @param bias_affirmative Mean internal bias, on the additive scale, among published affirmative
#'   studies. The bias has the same units as \code{yi}.
#' @param bias_nonaffirmative Mean internal bias, on the additive scale, among published nonaffirmative
#'   studies. The bias has the same units as \code{yi}.
#' @param cluster Vector of the same length as the number of rows in the data,
#'   indicating which cluster each study should be considered part of (defaults
#'   to each study being in its own cluster).
#' @param biased Boolean indicating whether each study is considered internally
#'   biased, either single value used for all studies or a vector the same
#'   length as the number of rows in the data (defaults to all studies).
#' @param favor_positive \code{TRUE} if publication bias are
#'   assumed to favor significant positive estimates; \code{FALSE} if assumed to
#'   favor significant negative estimates.
#' @param alpha_select Alpha level at which an estimate's probability of being
#'   favored by p-hacking and/or by publication bias is assumed to change (i.e.,
#'   the threshold at which study investigators, journal editors, etc., consider
#'   an estimate to be significant).
#' @param ci_level Confidence interval level (as proportion) for the corrected
#'   point estimate. (The alpha level for inference on the corrected point
#'   estimate will be calculated from \code{ci_level}.)
#' @param small Should inference allow for a small meta-analysis? We recommend
#'   always using \code{TRUE}.
#'
#' @return
#' @export
#'
#' @examples
#' # no internal bias
#' corrected_meta_mbma(yi = meta_meat$yi,
#'                     vi = meta_meat$vi,
#'                     bias_affirmative = 0,
#'                     bias_nonaffirmative = 0,
#'                     selection_ratio = 4,
#'                     biased = !meta_meat$randomized)
#'
#' # internal bias
#' corrected_meta_mbma(yi = meta_meat$yi,
#'                     vi = meta_meat$vi,
#'                     bias_affirmative = log(1.5),
#'                     bias_nonaffirmative = log(1.1),
#'                     selection_ratio = 4,
#'                     biased = !meta_meat$randomized)
#'
#' # treat all studies as biased, not just non-randomized ones
#' corrected_meta_mbma(yi = meta_meat$yi,
#'                     vi = meta_meat$vi,
#'                     bias_affirmative = log(1.5),
#'                     bias_nonaffirmative = log(1.1),
#'                     selection_ratio = 4,
#'                     biased = TRUE)
corrected_meta_mbma <- function(yi,
                                vi,
                                sei,
                                selection_ratio,
                                bias_affirmative,
                                bias_nonaffirmative,
                                cluster = 1:length(yi),
                                biased = TRUE,
                                favor_positive = TRUE,
                                alpha_select = 0.05,
                                ci_level = 0.95,
                                small = TRUE) {

  # stopifnot(length(cluster) == nrow(dat))
  # stopifnot(length(biased) == 1 || length(biased) == nrow(dat))
  # stopifnot("yi" %in% names(dat))
  # stopifnot("vi" %in% names(dat))
  # stopifnot("pval" %in% names(dat))
  stopifnot(length(bias_affirmative) == 1)
  stopifnot(length(bias_nonaffirmative) == 1)

  if (!favor_positive) yi <- -yi

  if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(vi)) vi <- sei ^ 2
  if (missing(sei)) sei <- sqrt(vi)
  if (length(sei) != length(yi)) stop(
    "Length of 'vi' or 'sei' must match that of 'yi'."
  )
  if (any(sei < 0)) stop("vi or sei should never be negative.")

  tcrit <- qnorm(1 - alpha_select / 2)

  dat <- tibble(yi = yi, vi = vi, sei = sei, biased = biased, cluster = cluster) |>
    mutate(affirmative = (.data$yi / .data$sei) > tcrit)

  p_affirm_pub <- dat |> filter(.data$biased) |> pull(.data$affirmative) |> mean()
  p_nonaffirm_pub <- 1 - p_affirm_pub
  denominator <- p_affirm_pub + selection_ratio * p_nonaffirm_pub
  mhat_b <- (p_nonaffirm_pub * selection_ratio * bias_nonaffirmative +
               p_affirm_pub * bias_affirmative ) / denominator

  dat <- dat |> mutate(yi_adj = .data$yi - .data$biased * mhat_b)

  tau2 <- metafor::rma.uni(yi = dat$yi_adj, vi = dat$vi)$tau2

  dat <- dat |>
    mutate(weight = if_else(.data$affirmative, 1, selection_ratio),
           userweight = .data$weight / (.data$vi + tau2))

  meta_mbma <- robumeta::robu(yi_adj ~ 1,
                              data = dat,
                              studynum = cluster,
                              var.eff.size = dat$vi,
                              userweights = dat$userweight,
                              small = small)

  alpha <- 1 - ci_level
  stats <- meta_mbma$reg_table |>
    mutate(ci_width = qt(1 - alpha / 2, .data$dfs) * .data$SE,
           ci_lower = .data$b.r - .data$ci_width,
           ci_upper = .data$b.r + .data$ci_width) |>
    select(mu_hat = .data$b.r, .data$ci_lower, .data$ci_upper, p_value = .data$prob)

  vals <- list(selection_ratio = selection_ratio,
               bias_affirmative = bias_affirmative,
               bias_nonaffirmative = bias_nonaffirmative,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small)

  results <- list(data = dat, values = vals, stats = stats, fit = meta_mbma)
  class(results) <- "metabias"
  return(results)

}


#' E-value for meta-analysis with multiple biases
#'
#' @inheritParams corrected_meta_mbma
#' @param q The attenuated value to which to shift the point estimate or CI.
#'   Should be specified on the same scale as \code{dat$yi} (e.g., if
#'   \code{dat$yi} is on the log-RR scale, then \code{q} should be as well).
#' @param bias_grid_hi The largest value of \code{bias}, on the additive scale, that should be included
#'   in the grid search. The bias has the same units as \code{yi}.
#' @param evalue_transformation TODO MM: Please see multi_evalue_example.R
#'
#' @return
#' @export
#'
#' @examples
#' evalue_mbma(yi = meta_meat$yi,
#'             vi = meta_meat$vi,
#'             selection_ratio = 4,
#'             biased = !meta_meat$randomized)
evalue_mbma <- function(yi,
                        vi,
                        sei,
                        selection_ratio,
                        q = 0,
                        cluster = 1:length(yi),
                        biased = TRUE,
                        favor_positive = TRUE,
                        alpha_select = 0.05,
                        ci_level = 0.95,
                        small = TRUE,
                        bias_grid_hi = 20,
                        # TODO MM: Please see multi_evalue_example.R
                        evalue_transformation = function(b) b + sqrt(b ^ 2 - b)) {


  if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(sei)) corrected_fun <- purrr::partial(corrected_meta_mbma, vi = vi)
  else corrected_fun <- purrr::partial(corrected_meta_mbma, sei = sei)

  compute_eb <- function(opt_col) {
    bias_factor <- function(bias_shared) {
      # all arguments passed through except biases being set to bias_shared
      corrected <- corrected_fun(yi = yi,
                                 bias_affirmative = bias_shared,
                                 bias_nonaffirmative = bias_shared,
                                 selection_ratio = selection_ratio,
                                 cluster = cluster,
                                 biased = biased,
                                 favor_positive = favor_positive,
                                 alpha_select = alpha_select,
                                 ci_level = ci_level,
                                 small = small)
      return(abs(corrected$stats[[opt_col]] - q))
    }

    opt <- optimize(f = bias_factor, interval = c(0, bias_grid_hi))
    if (abs(opt$minimum - bias_grid_hi) < 0.0001 & opt$objective > 0.0001)
      return(paste(">", bias_grid_hi))
    return(opt$minimum)
  }

  bias_est <- compute_eb("mu_hat")
  # MM: Below, we'd call multi_evalue instead of evalue_transformation
  evalue_est <- evalue_transformation(exp(bias_est))

  bias_ci <- compute_eb("ci_lower")
  evalue_ci <- evalue_transformation(exp(bias_ci))

  stats <- list(bias_est = bias_est, bias_ci = bias_ci,
                evalue_est = evalue_est, evalue_ci = evalue_ci)

  vals <- list(selection_ratio = selection_ratio,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small,
               q = q,
               bias_grid_hi = bias_grid_hi)

  results <- list(data = NULL, values = vals, stats = stats, fit = list())
  class(results) <- "metabias"
  return(results)

}
