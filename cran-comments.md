# RESUBMISSION: ebirdst 0.2.0

- `compute_ppms()` example had long run time, wrapped in \donttest{}

# ebirdst 0.2.0

- change maintainer to Matthew Strimas-Mackey
- update to access 2019 status and trends data
- partial dependence data no longer available, all references to PDs removed
- bug fix: `load_raster()` gave incorrect names to seasonal rasters
- bug fix: didn't properly implement quantile binning
- `date_to_st_week()` gets the status and trends week for a give vector of dates

## Test environments

- local OS X install, R 3.6.2
- OS X (travis-ci), R 3.6.2
- ubuntu 16.04 (travis-ci), R 3.6.2
- Windows (appveyor), R 3.6.2
- Rhub
  - Ubuntu Linux 16.04 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note
- NOTE: changed maintainer from Tom Auer <mta45@cornell.edu> to Matthew Strimas-Mackey <mes335@cornell.edu>
