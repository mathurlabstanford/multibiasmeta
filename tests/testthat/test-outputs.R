test_that("correct output structure for multibias_meta", {

  # change values of all params and opts
  params <- list(selection_ratio = 1, bias_affirmative = log(1.5),
                 bias_nonaffirmative = log(1.1), favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE)

  meta <- rlang::exec(multibias_meta, yi = -meta_meat$yi, vi = meta_meat$vi,
                      cluster = meta_meat$cluster,
                      biased = !meta_meat$randomized, !!!params,
                      return_worst_meta = TRUE, return_pubbias_meta = TRUE)

  expect_s3_class(meta, "metabias")
  expect_named(meta, c("data", "values", "stats", "fits"))

  expect_s3_class(meta$data, "data.frame")
  expect_equal(nrow(meta$data), nrow(meta_meat))
  expect_named(meta$data, meta_names("data"))
  expect_equal(meta$data$yi, meta_meat$yi)
  expect_equal(meta$data$vi, meta_meat$vi)
  expect_equal(meta$data$cluster, meta_meat$cluster)
  expect_equal(meta$data$biased, !meta_meat$randomized)

  expect_type(meta$values, "list")
  expect_named(meta$values, meta_names("values"))
  purrr::walk(names(params), \(p) expect_equal(params[[p]], meta$values[[p]]))

  expect_s3_class(meta$stats, "data.frame")
  expect_equal(nrow(meta$stats), 3)
  expect_named(meta$stats, meta_names("stats"))
  expect_true(all(meta$stats$se < 1))
  expect_true(all(meta$stats$ci_lower < meta$stats$estimate))
  expect_true(all(meta$stats$ci_upper > meta$stats$estimate))
  expect_true(all(meta$stats$p_value < 1))

  expect_length(meta$fit, 3)
  expect_named(meta$fit, meta_names("fits"))
  expect_s3_class(meta$fit$multibias, "robu")
  expect_s3_class(meta$fit$pubbias, "robu")
  expect_s3_class(meta$fit$worst_case, "robu")

})

test_that("correct output structure for multibias_evalue", {

  # change values of all params and opts
  params <- list(selection_ratio = 4, q = 0.1, favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE,
                 bias_max = 10)

  evalue <- rlang::exec(multibias_evalue, yi = -meta_meat$yi, vi = meta_meat$vi,
                        cluster = meta_meat$cluster,
                        biased = !meta_meat$randomized, !!!params,
                        assumed_bias_type = list(EValue::confounding()))

  expect_s3_class(evalue, "metabias")
  expect_named(evalue, c("data", "values", "stats", "fits"))

  expect_s3_class(evalue$data, "data.frame")
  expect_equal(nrow(evalue$data), nrow(meta_meat))
  expect_named(evalue$data, evalue_names("data"))
  expect_equal(evalue$data$yi, meta_meat$yi)
  expect_equal(evalue$data$vi, meta_meat$vi)
  expect_equal(evalue$data$cluster, meta_meat$cluster)
  expect_equal(evalue$data$biased, !meta_meat$randomized)

  expect_type(evalue$values, "list")
  expect_named(evalue$values, evalue_names("values"))
  expect_equal(evalue$values, params)

  expect_s3_class(evalue$stats, "data.frame")
  expect_named(evalue$stats, evalue_names("stats"))
  expect_gte(evalue$stats$bias_est, evalue$stats$bias_ci)
  expect_gte(evalue$stats$evalue_est, evalue$stats$evalue_ci)

  expect_length(evalue$fit, 0)

})
