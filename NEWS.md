# ebirdst 1.2021.3

- patch to fix a bug introduced in last release causing missing config files from data downloads [issue #44]

# ebirdst 1.2021.2

- fix bug causing species with same base code to be downloaded together, e.g. leafly also downloads leafly2 [issue #43]

# ebirdst 1.2021.1

- fix bug with extent in `load_fac_map_parameters()`, GitHub issue #40
- use dynamic PAT cutoff in PPM calculations
- update species list to account for second release of eBird data this year

# ebirdst 1.2021.0

- update for the v2021 eBird Status and Trends data

# ebirdst 1.2020.1

- CRAN checks found files created and left behind in ~/Desktop, relocated test files to tempdir() and deleting after test completion with withr::defer()  

# ebirdst 1.2020.0

- major update to align with the new eBird Status and Trends API
- update to align with the 2020 eBird Status Data Products
- transition from rappdirs to tools::R_user_dir() for handling download directories
- all new vignettes

# ebirdst 0.3.5

- bug fix: API update is causing all data downloads to fail

# ebirdst 0.3.4

- rename master branch to main on GitHub requires different download path for example data

# ebirdst 0.3.3

- move example data to GitHub

# ebirdst 0.3.2

- again try to prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.1

- prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.1

- prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.0

- add support for new data structures used for 2020 eBird Status and Trends
- functionality to handle partial dependence data added
- overhaul of package API to be more intuitive and streamlined
- all documentation and vignettes updated

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
