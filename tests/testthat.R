# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(mbma)

test_check("mbma")



# Test 1
# line 75 of helper_applied_MBMA.R
# compare to direct calculation of user-weighted estimator
expect_equal( as.numeric( meta.mbma$b.r ),
              sum( dat$yi.adj.est * ( weights / (dat$vi + t2hat.naive) ) ) / sum( weights / (dat$vi + t2hat.naive) ) )




# Test 2
# MBMA should agree with SAPB when there's no confounding
# line 205 of analyze_applied_MBMA.R
meta.SAPB.check = PublicationBias::corrected_meta(yi = d$yi,
                                                  vi = d$vi,
                                                  cluster = d$cluster,
                                                  eta = eta,
                                                  model = "robust",
                                                  favor.positive = TRUE)

expect_equal( meta.SAPB.row$Mhat, transf(meta.SAPB.check$est) )
expect_equal( meta.SAPB.row$MLo, transf(meta.SAPB.check$lo) )
expect_equal( meta.SAPB.row$MHi, transf(meta.SAPB.check$hi) )
expect_equal( meta.SAPB.row$MPval, meta.SAPB.check$pval )




# Test 3
# if all studies are confounded, E-value for estimate should just be equal to Mhat from
#  SAPB since lambda = 1
# line 268 of analyze_applied_MBMA.R
expect_equal( exp( res2$EB_est[ res2$eta_assumed == eta ] ),
              res$Mhat[ res$Analysis == "mbma-unadj" ],
              tol = 0.001 )

expect_equal( exp( res2$EB_ci[ res2$eta_assumed == eta ] ),
              res$MLo[ res$Analysis == "mbma-unadj" ],
              tol = 0.001 )



