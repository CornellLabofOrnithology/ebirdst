# ebirdst 1.2020.1

- CRAN checks found files created and left behind in ~/Desktop, relocated test files to tempdir() and deleting after test completion with withr::defer()

## Test environments

- local OS X install, R 4.2
- Windows (github actions), R 4.2
- ubuntu 20.04 (github actions), R 4.2
- win-builder (devel and release)
- R-hub (Solaris), R 4.2

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (1.2020.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2020, so we've included that year in our version number to indicate that this package only works with the 2020 version of the data.
- NOTE: win-builder highlights the following two links as "Service Unavailable"; however, I've repeatedly checked them and they're both working.
  - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019EA000658
  - https://doi.org/10.1002/eap.2056
