multibias_evalue(yi = meta_meat$yi,
                 vi = meta_meat$vi,
                 biased = !meta_meat$randomized,
                 selection_ratio = 4)

# specify confounding as internal bias
multibias_evalue(yi = meta_meat$yi,
                 vi = meta_meat$vi,
                 biased = !meta_meat$randomized,
                 selection_ratio = 4,
                 assumed_bias_type = list(EValue::confounding()))
