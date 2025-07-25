% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetPanels.R
\name{GetPanels}
\alias{GetPanels}
\title{Extract channels and markers from FCS files}
\usage{
GetPanels(fnames)
}
\arguments{
\item{fnames}{string vector. Full paths to FCS files}
}
\value{
named \code{list} per file of \code{data.frame}s with columns \emph{'Index'},
\emph{'Channel'}, and \emph{'Marker'}
}
\description{
Retrieves channel and marker indices and names per FCS sample without parsing
the files themselves.
By only reading the FCS file headers, this offers a speed-up and allows to
check for any inconsistencies between sample panels to prevent errors or
confusion downstream.
}
