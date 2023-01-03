# ebirdst 1.2021.1

- fix bug causing species with same base code to be downloaded together, e.g. leafly also downloads leafly2, GitHub issue #43

# ebirdst 1.2021.1

- fix bug with extent in `load_fac_map_parameters()`, GitHub issue #40
- use dynamic PAT cutoff in PPM calculations
- update species list to account for second release of eBird data this year

## Test environments

- local OS X install, R 4.2
- Windows (github actions), R 4.2
- ubuntu 20.04 (github actions), R 4.2
- win-builder (devel and release)
- R-hub (Solaris), R 4.2

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (1.2021.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2021, so we've included that year in our version number to indicate that this package only works with the 2021 version of the data.
- NOTE: win-builder highlights the following two links as "Service Unavailable"; however, I've repeatedly checked them and they're both working.
  - https://doi.org/10.1029/2019EA000658
  - https://doi.org/10.1002/eap.2056
