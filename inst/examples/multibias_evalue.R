# specify confounding as internal bias
evalues <- multibias_evalue(yi = meta_meat$yi,
                            vi = meta_meat$vi,
                            biased = !meta_meat$randomized,
                            selection_ratio = 4)
evalues$stats

# specify confounding as internal bias
evalues_confounding <- multibias_evalue(yi = meta_meat$yi,
                                        vi = meta_meat$vi,
                                        biased = !meta_meat$randomized,
                                        selection_ratio = 4,
                                        assumed_bias_type = list(EValue::confounding()))
evalues_confounding$stats
