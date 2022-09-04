
library(EValue)

# ~ Example 1: confounding ------------------------

# instead of the argument evalue_transformation, we could have the user (optionally) pass a bias object as follows, per the EValue package
# confounding() requires no additional arguments, but some other biases do (see Example 2 below)
# relevant vignette: https://cran.r-project.org/web/packages/EValue/vignettes/multiple-bias.html
my_biases = confounding()

# equivalent to evalue_transformation in current version of corrected_meta_mbma
# RR(4) would be replaced with exp(stats$mu_hat); see also nuance #2 in later section here
summary( multi_evalue( est = RR(4), biases = my_biases ) )

# would be good to also return to the user the interpretation of the E-value
#  sensitivity parameter, as follows:
biases = multi_bias( my_biases )
summary(biases)$output

# the resulting parameters are described here (let's cite in Details):
#  https://cran.r-project.org/web/packages/EValue/vignettes/multiple-bias.html

# and more fully in the following papers, which would be good to cite in Details:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4820664/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6553568/
# https://academic.oup.com/aje/article/188/10/1823/5506602?casa_token=5ZyiVJp9_5UAAAAA:40rpOH1mRz0IDeRJ35atRRk9x6MJgIHMNOxLCcsnfouzN3qWXrght0XVWNIHQcRwWP1Bhgl8vY9B


# ~ Example 2: outcome misclassification ------------------------

# for outcome misclassification, the transformation is just the identity,
#  so we'll get back the original RR
my_biases = misclassification("outcome")

multi_evalue( est = RR(4), biases = my_biases )

biases = multi_bias( my_biases )
summary(biases)
summary(biases)$output

# ~ Nuances ------------------------

# 1.) The bounds that multi_evalue uses internally work with effect sizes on the RR
#   scale, but EValue::convert_measures has fns to convert other effect sizes
#   (e.g., odds ratio, hazard ratio, etc.) to RRs.
#
#  corrected_meta_mbma can do all of its work *prior* to converting to the E-value
#  scale without making any transformations to the effect sizes. For the E-value,
#  we need to know what scale the effect sizes are on.
#
# I think the simplest solution would be to simply explain in the docs that *if* an E-value transformation is being used, the yi argument must be on the log-RR scale. If the yi's are not already on that scale, the user can use EValue::convert_measures themselves. I worry that if we try to build the conversions into corrected_meta_mbma itself, it will be too many arguments and too confusing.


# 2.) multi_evalue is mainly designed for considering multiple internal biases within a single study,
#  but as a special case also handles considering a single bias (as in examples above).
#  For corrected_meta_mbma, I think it would usually be overkill and too hard to interpret to
#  consider multiple internal biases on top of pub bias. So I think we don't have to *prevent* the user
#  from passing multiple biases, but the docs and examples should focus on cases with just 1 internal bias.










































