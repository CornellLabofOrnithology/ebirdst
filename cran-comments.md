# Resubmission

Fixing link redirects. Removing unused dependencies.

# ebirdst 0.2.2

- add support for variable ensemble support in `compute_ppms()`
- wrap all examples in dontrun to avoid writing files locally when compiling package

## Test environments

- local OS X install, R 4.0.3
- ubuntu 16.04 (travis-ci), R 4.0.3
- Windows (appveyor), R 4.0.3
- Rhub
  - Ubuntu Linux 16.04 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

- NOTE Namespaces in Imports field not imported from: ‘rappdirs’
rappdirs is used as a function argument for ebirdst_download()
