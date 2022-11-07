test_that("catch bad inputs", {

  # selection ratio must be at least 1
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi = meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 0,
                   bias_affirmative = 0, bias_nonaffirmative = 0))

  # biases must be one value
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi = meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 1:4, bias_nonaffirmative = 0))
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi = meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 1:4))

  # need at least one of vi and sei
  expect_error(
    multibias_meta(yi = meta_meat$yi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0))

  # yi, vi, sei, and cluster should all be the same length
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi = rep(1, 10),
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0))
  expect_error(
    multibias_meta(yi = meta_meat$yi, sei = rep(1, 10),
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0))
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi =  meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0,
                   cluster = rep(1, 10)))

  # biased indicator should be either one value or length yi
  expect_error(
    multibias_meta(yi = meta_meat$yi, vi =  meta_meat$vi,
                   biased = rep(1, 10), selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0)
  )

  # sei values must be positive
  expect_error(
    multibias_meta(yi = meta_meat$yi, sei = -meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0))

  # warning for favored direction opposite of naive estimate
  expect_warning(
    multibias_meta(yi = -meta_meat$yi, vi = meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0))
  expect_warning(
    multibias_meta(yi = meta_meat$yi, sei = meta_meat$vi,
                   biased = !meta_meat$randomized, selection_ratio = 1,
                   bias_affirmative = 0, bias_nonaffirmative = 0,
                   favor_positive = FALSE))

})
