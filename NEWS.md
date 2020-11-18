# ebirdst 1.0.0

- add support for new data structures used for 2020 eBird Status and Trends
- overhaul of 

# ebirdst 0.2.2

- add support for variable ensemble support in `compute_ppms()`

# ebirdst 0.2.1

- bug fix: corrected date types in seasonal definitions
- bug fix: fixed possibility that ebirdst_extent could produce invalid date (day 366 of 2015)
- added import of pipe operator
- `velox` was archived, removed dependency from Suggests
- `fasterize` was archived, removed dependency from Imports

# ebirdst 0.2.0

- change maintainer to Matthew Strimas-Mackey
- update to access 2019 status and trends data
- partial dependence data no longer available, all references to PDs removed
- bug fix: `load_raster()` gave incorrect names to seasonal rasters
- bug fix: didn't properly implement quantile binning
- `date_to_st_week()` gets the status and trends week for a give vector of dates

# ebirdst 0.1.0

- first CRAN release
