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
  if (missing(vi) && missing(sei)) stop("Must specify 'vi' or 'sei' argument.")
  if (missing(vi)) vi <- sei ^ 2
  if (missing(sei)) sei <- sqrt(vi)
  if (length(sei) != length(yi)) stop(
    "Length of 'vi' or 'sei' must match that of 'yi'."
  )
  if (any(sei < 0)) stop("vi or sei should never be negative.")
  stopifnot(length(vi) == length(yi))
  stopifnot(length(sei) == length(yi))

  # warn if naive estimate is in opposite direction than favor_positive
  naive_pos <- metafor::rma(yi, vi, method = "FE")$beta > 0
  if (naive_pos != favor_positive)
    warning("Favored direction is opposite of the pooled estimate.")

  # flip yi if not favor_positive
  if (!favor_positive) yi <- -yi

  # combine data vectors and get each study's affirmative status
  tcrit <- qnorm(1 - alpha_select / 2)
  dat <- tibble(yi, vi, sei, biased, cluster) |>
    mutate(affirmative = (.data$yi / .data$sei) > tcrit)

  # adjust yi values based on biases and selection ratio
  p_affirm_pub <- dat |>
    filter(.data$biased) |>
    pull(.data$affirmative) |> mean()
  p_nonaffirm_pub <- 1 - p_affirm_pub
  denominator <- p_affirm_pub + selection_ratio * p_nonaffirm_pub
  mhat_b <- (p_nonaffirm_pub * selection_ratio * bias_nonaffirmative +
               p_affirm_pub * bias_affirmative) / denominator
  dat <- dat |> mutate(yi_adj = .data$yi - .data$biased * mhat_b)

  # compute weights
  tau2 <- metafor::rma.uni(yi = dat$yi_adj, vi = dat$vi)$tau2
  dat <- dat |>
    mutate(weight = if_else(.data$affirmative, 1, selection_ratio),
           userweight = .data$weight / (.data$vi + tau2))

  # fit meta analysis using adjusted yi values and weights
  meta_multibias <- robumeta::robu(yi_adj ~ 1,
                                   data = dat,
                                   studynum = cluster,
                                   var.eff.size = dat$vi,
                                   userweights = dat$userweight,
                                   small = small)

  # extract stats from meta analysis
  stats <- metabias::robu_ci(meta_multibias, ci_level) |>
    mutate(model = model_label, .before = everything())

  fit <- list()
  fit[[model_label]] <- meta_multibias

  metabias::metabias(data = dat, values = vals, stats = stats, fits = fit)

}


#' Correction for meta-analysis with multiple biases
#'
#' @inheritParams metabias::params
#' @param biased Boolean indicating whether each study is considered internally
#'   biased; either single value used for all studies or a vector the same
#'   length as the number of rows in the data (defaults to all studies).
#' @param bias_affirmative Mean internal bias, on the additive scale, among
#'   published affirmative studies. The bias has the same units as `yi`.
#' @param bias_nonaffirmative Mean internal bias, on the additive scale, among
#'   published nonaffirmative studies. The bias has the same units as `yi`.
#' @param return_worst_meta Boolean indicating whether the worst-case
#'   meta-analysis of only the nonaffirmative studies be returned.
#' @param return_pubbias_meta Boolean indicating whether a meta-analysis
#'   correcting for publication but not for confounding be returned.
#'
#' @return An object of class [metabias::metabias()].
#'
#' @references \insertRef{mathur2022multibias}{metabias}
#'
#' @export
#' @example inst/examples/multibias_meta.R
multibias_meta <- function(yi, # data
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
    meta_pubbias <- exec(get_multibias_meta, !!!pubbias_args,
                         model_label = "pubbias")

    meta_multibias$stats <- bind_rows(meta_multibias$stats, meta_pubbias$stats)
    meta_multibias$fits <- append(meta_multibias$fit, meta_pubbias$fit)
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

    worst_stats <- metabias::robu_ci(meta_worst, ci_level) |>
      mutate(model = worst_label, .before = everything())

    meta_multibias$stats <- bind_rows(meta_multibias$stats, worst_stats)
    meta_multibias$fits[[worst_label]] <- meta_worst
  }

  return(meta_multibias)
}
