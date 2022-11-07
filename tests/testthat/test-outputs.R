test_that("correct output structure for multibias_meta", {

  # change values of all params and opts
  params <- list(selection_ratio = 1, bias_affirmative = log(1.5),
                 bias_nonaffirmative = log(1.1), favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE)

  mm <- rlang::exec(multibias_meta, yi = -meta_meat$yi, vi = meta_meat$vi,
                    cluster = meta_meat$cluster,
                    biased = !meta_meat$randomized, !!!params,
                    return_worst_meta = TRUE, return_pubbias_meta = TRUE)

  expect_s3_class(mm, "metabias")
  expect_length(mm, 4)
  expect_named(mm, c("data", "values", "stats", "fit"))

  expect_s3_class(mm$data, "data.frame")
  expect_equal(nrow(mm$data), nrow(meta_meat))
  expect_equal(mm$data$yi, meta_meat$yi)
  expect_equal(mm$data$vi, meta_meat$vi)
  expect_equal(mm$data$cluster, meta_meat$cluster)
  expect_equal(mm$data$biased, !meta_meat$randomized)

  expect_type(mm$values, "list")
  expect_equal(mm$values, params)

  expect_s3_class(mm$stats, "data.frame")
  expect_equal(nrow(mm$stats), 3)
  expect_named(mm$stats, c("model", "mu_hat", "se", "ci_lower", "ci_upper", "p_value"))
  expect_true(all(mm$stats$se < 1))
  expect_true(all(mm$stats$ci_lower < mm$stats$mu_hat))
  expect_true(all(mm$stats$ci_upper > mm$stats$mu_hat))
  expect_true(all(mm$stats$p_value < 1))

  expect_length(mm$fit, 3)
  expect_named(mm$fit, c("multibias", "pubbias", "worst_case"))
  expect_s3_class(mm$fit$multibias, "robu")
  expect_s3_class(mm$fit$pubbias, "robu")
  expect_s3_class(mm$fit$worst_case, "robu")

})

test_that("correct output structure for evalue", {

  # change values of all params and opts
  params <- list(selection_ratio = 4, q = 0.1, favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE,
                 bias_max = 10)

  ev <- rlang::exec(multibias_evalue, yi = -meta_meat$yi, vi = meta_meat$vi,
                    cluster = meta_meat$cluster,
                    biased = !meta_meat$randomized, !!!params,
                    assumed_bias_type = list(EValue::confounding()))

  expect_s3_class(ev, "metabias")
  expect_length(ev, 4)
  expect_named(ev, c("data", "values", "stats", "fit"))

  expect_s3_class(ev$data, "data.frame")
  expect_equal(nrow(ev$data), nrow(meta_meat))
  expect_equal(ev$data$yi, meta_meat$yi)
  expect_equal(ev$data$vi, meta_meat$vi)
  expect_equal(ev$data$cluster, meta_meat$cluster)
  expect_equal(ev$data$biased, !meta_meat$randomized)

  expect_type(ev$values, "list")
  expect_equal(ev$values, params)

  expect_type(ev$stats, "list")
  expect_named(ev$stats, c("bias_est", "bias_ci", "evalue_est", "evalue_ci"))
  expect_gte(ev$stats$bias_est, ev$stats$bias_ci)
  expect_gte(ev$stats$evalue_est, ev$stats$evalue_ci)

  expect_length(ev$fit, 0)

})
