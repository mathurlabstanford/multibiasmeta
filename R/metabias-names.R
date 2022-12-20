#' @keywords internal
meta_names <- function(component) {
  names_list <- list(
    data = c("yi", "vi", "sei", "biased", "cluster", "affirmative", "yi_adj",
             "weight", "userweight"),
    values = c("selection_ratio", "bias_affirmative", "bias_nonaffirmative",
               "favor_positive", "alpha_select", "ci_level", "small"),
    stats = c("model", "estimate", "se", "ci_lower", "ci_upper", "p_value"),
    fits = c("multibias", "pubbias", "worst_case"))
  names_list[[component]]
}

#' @keywords internal
meta_names_str <- function(component) {
  cnames <- meta_names(component)
  paste(paste0("`", cnames, "`"), collapse = ", ")
}

#' @keywords internal
evalue_names <- function(component) {
  names_list <- list(
    data = c("yi", "vi", "sei", "cluster", "biased"),
    values = c("selection_ratio", "q", "favor_positive", "alpha_select",
               "ci_level", "small", "bias_max"),
    stats = c("bias_est", "bias_ci", "evalue_est", "evalue_ci"))
  names_list[[component]]
}

#' @keywords internal
evalue_names_str <- function(component) {
  cnames <- evalue_names(component)
  paste(paste0("`", cnames, "`"), collapse = ", ")
}
