% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebirdst-loading.R
\name{date_to_st_week}
\alias{date_to_st_week}
\title{Get the status and trends week that a date falls into}
\usage{
date_to_st_week(dates)
}
\arguments{
\item{dates}{a vector of dates.}
}
\value{
An integer vector of weeks numbers from 1-52.
}
\description{
Get the status and trends week that a date falls into
}
\examples{
d <- as.Date(c("2016-04-08", "2018-12-31", "2014-01-01"))
date_to_st_week(d)
}