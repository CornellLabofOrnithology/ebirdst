## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "\\textwidth", 
                      fig.height = 5, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes local and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----libraries----------------------------------------------------------------
library(ebirdst)
library(raster)
library(velox)
library(sf)
library(smoothr)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# resolve namespace conflicts
select <- dplyr::select

## ----st-download--------------------------------------------------------------
# download to a temp directory for the vigette
# in practice, change to permanent directory the status and trends downloads
sp_path <- ebirdst_download(species = "example_data")
# load the abundance data
# this automaticaaly labels layers with their dates
abd <- load_raster("abundance", path = sp_path)

