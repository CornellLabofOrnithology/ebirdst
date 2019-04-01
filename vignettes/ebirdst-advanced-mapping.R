## ----setup, include=FALSE------------------------------------------------
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

## ----libraries-----------------------------------------------------------
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

