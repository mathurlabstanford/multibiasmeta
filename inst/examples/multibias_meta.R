# publication bias without internal bias
meta_0 <- multibias_meta(yi = meta_meat$yi,
                         vi = meta_meat$vi,
                         selection_ratio = 4,
                         bias_affirmative = 0,
                         bias_nonaffirmative = 0)
meta_0$stats

# publication bias and internal bias in the non-randomized studies
meta_4 <- multibias_meta(yi = meta_meat$yi,
                         vi = meta_meat$vi,
                         biased = !meta_meat$randomized,
                         selection_ratio = 4,
                         bias_affirmative = log(1.5),
                         bias_nonaffirmative = log(1.1))
meta_4$stats

# treat all studies as biased, not just non-randomized ones
meta_all <- multibias_meta(yi = meta_meat$yi,
                           vi = meta_meat$vi,
                           biased = TRUE,
                           selection_ratio = 4,
                           bias_affirmative = log(1.5),
                           bias_nonaffirmative = log(1.1))
meta_all$stats
