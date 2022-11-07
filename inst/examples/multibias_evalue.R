multibias_evalue(yi = meta_meat$yi,
                 vi = meta_meat$vi,
                 selection_ratio = 4,
                 biased = !meta_meat$randomized)

# specify confounding as internal bias
multibias_evalue(yi = meta_meat$yi,
                 vi = meta_meat$vi,
                 selection_ratio = 4,
                 biased = !meta_meat$randomized,
                 assumed_bias_type = list(EValue::confounding()))
