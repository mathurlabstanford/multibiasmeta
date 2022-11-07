# publication bias without internal bias
multibias_meta(yi = meta_meat$yi,
               vi = meta_meat$vi,
               biased = !meta_meat$randomized,
               selection_ratio = 4,
               bias_affirmative = 0,
               bias_nonaffirmative = 0)

# publication bias and internal bias in the non-randomized studies
multibias_meta(yi = meta_meat$yi,
               vi = meta_meat$vi,
               biased = !meta_meat$randomized,
               selection_ratio = 4,
               bias_affirmative = log(1.5),
               bias_nonaffirmative = log(1.1))

# treat all studies as biased, not just non-randomized ones
multibias_meta(yi = meta_meat$yi,
               vi = meta_meat$vi,
               biased = TRUE,
               selection_ratio = 4,
               bias_affirmative = log(1.5),
               bias_nonaffirmative = log(1.1))
