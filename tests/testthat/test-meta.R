test_that("compare with direct calculation", {
  # compare to direct calculation of user-weighted estimator
  # line 75 of helper_applied_mbma.R
  meta <- multibias_meta(yi = meta_meat$yi,
                         vi = meta_meat$vi,
                         biased = !meta_meat$randomized,
                         selection_ratio = 4,
                         bias_affirmative = log(1.5),
                         bias_nonaffirmative = log(1.1))

  meta_estimate <- meta$stats$mu_hat
  meta_data <- meta$data
  t2hat <- metafor::rma.uni(yi = meta_data$yi_adj, vi = meta_data$vi)$tau2

  calculated_estimate <- sum(meta_data$yi_adj *
                               (meta_data$weight / (meta_data$vi + t2hat))) /
    sum(meta_data$weight / (meta_data$vi + t2hat))

  expect_equal(meta_estimate, calculated_estimate)
})


test_that("compare to to pubbias", {
  # multibiasmeta should agree with SAPB when there's no confounding
  # line 205 of analyze_applied_mbma.R
  meta <- multibias_meta(yi = meta_meat$yi,
                         vi = meta_meat$vi,
                         cluster = meta_meat$cluster,
                         biased = !meta_meat$randomized,
                         selection_ratio = 4,
                         bias_affirmative = 0,
                         bias_nonaffirmative = 0)

  meta_pubbias <- PublicationBias::pubbias_meta(yi = meta_meat$yi,
                                                vi = meta_meat$vi,
                                                cluster = meta_meat$cluster,
                                                selection_ratio = 4)

  expect_equal(meta$stats$mu_hat, meta_pubbias$stats$estimate)
  expect_equal(meta$stats$se, meta_pubbias$stats$se)
  expect_equal(meta$stats$ci_lower, meta_pubbias$stats$ci_lower)
  expect_equal(meta$stats$ci_upper, meta_pubbias$stats$ci_upper)
  expect_equal(meta$stats$p_value, meta_pubbias$stats$p_value)

})

test_that("compare estimate to evalue for all confounded", {
  # if all studies are confounded, E-value for estimate should just be equal to
  #   Mhat from SAPB since lambda = 1
  # line 268 of analyze_applied_mbma.R

  meta <- multibias_meta(yi = meta_meat$yi,
                         vi = meta_meat$vi,
                         biased = TRUE,
                         selection_ratio = 4,
                         bias_affirmative = 0,
                         bias_nonaffirmative = 0)

  evalue <- multibias_evalue(yi = meta_meat$yi,
                             vi = meta_meat$vi,
                             selection_ratio = 4,
                             biased = TRUE)

  expect_equal(meta$stats$mu_hat, evalue$stats$bias_est, tolerance = 0.001)
  expect_equal(meta$stats$ci_lower, evalue$stats$bias_ci, tolerance = 0.001)

})
