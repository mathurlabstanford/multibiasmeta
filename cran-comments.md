## Test environments

The package was tested with `R CMD check --as-cran` on the following platforms:

* Ubuntu Linux GCC (oldrel, release, devel),
* Windows x86_64  (oldrel, release, devel),
* MacOS (release).

## R CMD check results


There were no ERRORs or WARNINGs. The following NOTE was displayed:

```
Maintainer: ‘Peter Solymos <peter@analythium.io>’

Days since last update: 1
```

I am submitting an update because of the following follow up email from CRAN:

```
Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_multibiasmeta.html>.

Please correct before 2023-09-02 to safely retain your package on CRAN.

Packages in Suggests should be used conditionally: see 'Writing R Extensions'.
This needs to be corrected even if the missing package(s) become available.
It can be tested by checking with _R_CHECK_DEPENDS_ONLY_=true.

The CRAN Team
```

The check results indicated ERROR on Fedora:

```
Version: 0.2.0
Check: package dependencies
Result: NOTE
    Package suggested but not available for checking: ‘phacking’
Flavors: r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc
```

The phacking package is not yet available on Fedora, therefore we made the tests and tutorial that depend on phacking conditional on the package being present. We tested the new package with the `_R_CHECK_DEPENDS_ONLY_=true` setting as suggested.
