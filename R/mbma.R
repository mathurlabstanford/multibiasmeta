#' Add CI with specified ci_level to reg_table output from robumeta::robu()
#'
#' @keywords internal
robu_ci <- function(reg_table, ci_level) {
  alpha <- 1 - ci_level
  reg_table |>
    mutate(ci_width = qt(1 - alpha / 2, .data$dfs) * .data$SE,
           ci_lower = .data$b.r - .data$ci_width,
           ci_upper = .data$b.r + .data$ci_width) |>
    select(mu_hat = .data$b.r, .data$ci_lower, .data$ci_upper, p_value = .data$prob)
}

#' @keywords internal
get_multibias_meta <- function(yi,
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
                               small = TRUE,
                               model_label,
                               ...) {

  vals <- list(selection_ratio = selection_ratio,
               bias_affirmative = bias_affirmative,
               bias_nonaffirmative = bias_nonaffirmative,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small)

  # validate inputs
  stopifnot(selection_ratio >= 1)
  stopifnot(length(bias_affirmative) == 1)
  stopifnot(length(bias_nonaffirmative) == 1)
  stopifnot(length(cluster) == length(yi))
  stopifnot(length(biased) == 1 || length(biased) == length(yi))

  # resolve vi and sei
  if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(vi)) vi <- sei ^ 2
  if (missing(sei)) sei <- sqrt(vi)
  if (length(sei) != length(yi)) stop(
    "Length of 'vi' or 'sei' must match that of 'yi'."
  )
  if (any(sei < 0)) stop("vi or sei should never be negative.")

  # flip yi if not favor_positive
  if (!favor_positive) yi <- -yi

  # combine data vectors and get each study's affirmative status
  tcrit <- qnorm(1 - alpha_select / 2)
  dat <- tibble(yi, vi, sei, biased, cluster) |>
    mutate(affirmative = (.data$yi / .data$sei) > tcrit)

  # adjust yi values based on biases and selection ratio
  p_affirm_pub <- dat |> filter(.data$biased) |> pull(.data$affirmative) |> mean()
  p_nonaffirm_pub <- 1 - p_affirm_pub
  denominator <- p_affirm_pub + selection_ratio * p_nonaffirm_pub
  mhat_b <- (p_nonaffirm_pub * selection_ratio * bias_nonaffirmative +
               p_affirm_pub * bias_affirmative ) / denominator
  dat <- dat |> mutate(yi_adj = .data$yi - .data$biased * mhat_b)

  # compute weights
  tau2 <- metafor::rma.uni(yi = dat$yi_adj, vi = dat$vi)$tau2
  dat <- dat |>
    mutate(weight = if_else(.data$affirmative, 1, selection_ratio),
           userweight = .data$weight / (.data$vi + tau2))

  # fit meta analysis using adjusted yi values and weights
  meta_mbma <- robumeta::robu(yi_adj ~ 1,
                              data = dat,
                              studynum = cluster,
                              var.eff.size = dat$vi,
                              userweights = dat$userweight,
                              small = small)

  # extract stats from meta analysis
  stats <- robu_ci(meta_mbma$reg_table, ci_level) |>
    mutate(model = model_label, .before = everything())

  fit <- list()
  fit[[model_label]] <- meta_mbma
  results <- list(data = dat, values = vals, stats = stats, fit = fit)
  class(results) <- "metabias"
  return(results)

}


#' Correction for meta-analysis with multiple biases
#'
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances (i.e., squared standard errors) for
#'   the point estimates.
#' @param sei A vector of estimated standard errors for the point estimates.
#'   (Only one of \code{vi} or \code{sei} needs to be specified).
#' @param cluster Vector of the same length as the number of rows in the data,
#'   indicating which cluster each study should be considered part of (defaults
#'   to treating studies as independent; i.e., each study is in its own cluster).
#' @param biased Boolean indicating whether each study is considered internally
#'   biased; either single value used for all studies or a vector the same
#'   length as the number of rows in the data (defaults to all studies).
#' @param selection_ratio Ratio by which publication bias favors affirmative
#'   studies (i.e., studies with p-values less than \code{alpha_select} and
#'   estimates in the direction indicated by \code{favor_positive}).
#' @param bias_affirmative Mean internal bias, on the additive scale, among
#'   published affirmative studies. The bias has the same units as \code{yi}.
#' @param bias_nonaffirmative Mean internal bias, on the additive scale, among
#'   published nonaffirmative studies. The bias has the same units as \code{yi}.
#' @param favor_positive \code{TRUE} if publication bias are
#'   assumed to favor significant positive estimates; \code{FALSE} if assumed to
#'   favor significant negative estimates.
#' @param alpha_select Alpha level at which an estimate's probability of being
#'   favored by publication bias is assumed to change (i.e.,
#'   the threshold at which study investigators, journal editors, etc., consider
#'   an estimate to be significant).
#' @param ci_level Confidence interval level (as proportion) for the corrected
#'   point estimate. (The alpha level for inference on the corrected point
#'   estimate will be calculated from \code{ci_level}.)
#' @param small Should inference allow for a small meta-analysis? We recommend
#'   always using \code{TRUE}.
#' @param return_worst_meta Should the worst-case meta-analysis of only the
#'   nonaffirmative studies be returned?
#' @param return_pubbias_meta Should a meta-analysis correcting for publication
#'   but not for confounding be returned?
#'
#' @return
#'
#' @references
#' \insertRef{mathur2022}{mbma}
#'
#' @export
#' @examples
#' # publication bias without internal bias
#' multibias_corrected_meta(yi = meta_meat$yi,
#'                          vi = meta_meat$vi,
#'                          biased = !meta_meat$randomized,
#'                          selection_ratio = 4,
#'                          bias_affirmative = 0,
#'                          bias_nonaffirmative = 0)
#'
#' # publication bias and internal bias in the non-randomized studies
#' multibias_corrected_meta(yi = meta_meat$yi,
#'                          vi = meta_meat$vi,
#'                          biased = !meta_meat$randomized,
#'                          selection_ratio = 4,
#'                          bias_affirmative = log(1.5),
#'                          bias_nonaffirmative = log(1.1))
#'
#' # treat all studies as biased, not just non-randomized ones
#' multibias_corrected_meta(yi = meta_meat$yi,
#'                          vi = meta_meat$vi,
#'                          biased = TRUE,
#'                          selection_ratio = 4,
#'                          bias_affirmative = log(1.5),
#'                          bias_nonaffirmative = log(1.1))
multibias_corrected_meta <- function(yi, # data
                                     vi,
                                     sei,
                                     cluster = 1:length(yi),
                                     biased = TRUE,

                                     selection_ratio, # params
                                     bias_affirmative,
                                     bias_nonaffirmative,

                                     favor_positive = TRUE, # opts
                                     alpha_select = 0.05,
                                     ci_level = 0.95,
                                     small = TRUE,
                                     return_worst_meta = FALSE,
                                     return_pubbias_meta = FALSE) {

  args <- as.list(environment())
  args <- args |> discard(\(a) class(a) == "name")

  # get meta corrected for multibias for internal bias and publication bias
  meta_multibias <- exec(get_multibias_meta, !!!args, model_label = "multibias")

  # get meta corrected for only publication bias (no internal bias)
  if (return_pubbias_meta) {
    pubbias_args <- args
    pubbias_args[c("bias_affirmative", "bias_nonaffirmative")] <- 0
    meta_pubbias <- exec(get_multibias_meta, !!!pubbias_args, model_label = "pubbias")

    meta_multibias$stats <- bind_rows(meta_multibias$stats, meta_pubbias$stats)
    meta_multibias$fit <- append(meta_multibias$fit, meta_pubbias$fit)
  }

  # get worst case meta
  if (return_worst_meta) {
    worst_label <- "worst_case"
    dat <- meta_multibias$data
    meta_worst <- robumeta::robu(yi ~ 1,
                                 data = dat |> filter(!.data$affirmative),
                                 studynum = cluster,
                                 var.eff.size = vi,
                                 small = small)

    worst_stats <- robu_ci(meta_worst$reg_table, ci_level) |>
      mutate(model = worst_label, .before = everything())

    meta_multibias$stats <- bind_rows(meta_multibias$stats, worst_stats)
    meta_multibias$fit[[worst_label]] <- meta_worst
  }

  return(meta_multibias)
}

#' E-value for meta-analysis with multiple biases
#'
#' @inheritParams multibias_corrected_meta
#' @param q The attenuated value to which to shift the point estimate or CI.
#'   Should be specified on the same scale as \code{dat$yi} (e.g., if
#'   \code{dat$yi} is on the log-RR scale, then \code{q} should be as well).
#' @param bias_max The largest value of \code{bias}, on the additive scale,
#'   that should be included in the grid search. The bias has the same units as
#'   \code{yi}.
#' @param internal_biases List of biases to consider for computing evalues
#'   (objects of \code{bias} as returned by \code{EValue::confounding()},
#'   \code{EValue::selection()}, \code{EValue::misclassification()}) (defaults
#'   to NULL, i.e. agnostic as to the nature of the internal bias). If not NULL,
#'   the \code{yi} argument must be on the log-RR scale (if \code{yi} is not
#'   already on that scale, use \code{EValue::convert_measures()} to make it
#'   so).
#'
#' @return
#'
#' @details For more on the functions passed as \code{internal_biases}, see the
#'   EValue package vignette on multiple biases:
#'   <https://cran.r-project.org/web/packages/EValue/vignettes/multiple-bias.html>.
#'
#' @references
#' \insertRef{mathur2022}{mbma}
#'
#' \insertRef{ding2016}{mbma}
#'
#' \insertRef{smith2019}{mbma}
#'
#' \insertRef{vanderweele2019}{mbma}
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
#'                  internal_bias = list(EValue::))
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
                             internal_biases = NULL) { # TODO: clearer default

  # set up multibias_corrected_meta with either vi or sei passed
  if (missing(vi) & missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(sei)) corrected_fun <- partial(multibias_corrected_meta, vi = vi)
  else corrected_fun <- partial(multibias_corrected_meta, sei = sei)

  # compute expectation of bias for parameter (mu_hat or ci_lower)
  compute_eb <- function(param, tolerance = 0.0001) {

    # function to optimize: multibias_corrected_meta for a given bias
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

  if (is.null(internal_biases)) {
    evalue_est <- evalue_ci <- NULL
  } else {
    # turn list of internal biases into EValue::multi_bias object
    biases <- exec(EValue::multi_bias, !!!internal_biases)
    evalue_est <- summary(EValue::multi_evalue(est = EValue::RR(exp(bias_est)),
                                               biases = biases))
    # convert biases to evalues
    evalue_ci <- summary(EValue::multi_evalue(est = EValue::RR(exp(bias_ci)),
                                              biases = biases))
  }

  stats <- list(bias_est = bias_est, bias_ci = bias_ci,
                evalue_est = evalue_est, evalue_ci = evalue_ci)

  vals <- list(selection_ratio = selection_ratio,
               favor_positive = favor_positive,
               alpha_select = alpha_select,
               ci_level = ci_level,
               small = small,
               q = q,
               bias_max = bias_max)

  results <- list(data = NULL, values = vals, stats = stats, fit = list())
  class(results) <- "metabias"
  return(results)

}
