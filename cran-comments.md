## Test environments

The package was tested with `R CMD check --as-cran` on the following platforms:

* Ubuntu Linux GCC (oldrel, release, devel),
* Windows x86_64  (oldrel, release, devel),
* MacOS (release).

## R CMD check results


There were no ERRORs or WARNINGs. The following NOTE was displayed:

```
Maintainer: ‘Peter Solymos <peter@analythium.io>’

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2023-06-10 as issues were not corrected
    in time.
```

The package has new maintainer and is being submitted after archival.

The package was previously removed from CRAN because one of its dependencies, phacking, was archived.

```
Error: processing vignette ‘tutorial.Rmd’ failed with diagnostics:
    there is no package called 'phacking'
```

The phacking package is back on CRAN, thus the issues are resolved.
