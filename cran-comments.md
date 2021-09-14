# ebirdst 0.3.2

- again try to prevent tests and examples from leaving files behind to pass CRAN checks

## Test environments

- local OS X install, R 4.1
- OS X (github actions), R 4.1
- Windows (github actions), R 4.1
- ubuntu 14.04 (github actions), R 4.1
- win-builder (devel and release)
- R-hub (Solaris), R 4.1

## R CMD check results

0 errors | 0 warnings | 2 notes

notes:
- eBird was incorrectly identified as a misspelled word
- this release should fix the CRAN policy violation that caused the package to be archived
