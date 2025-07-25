% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParallelPreprocess.R
\name{ParallelPreprocess}
\alias{ParallelPreprocess}
\title{Parallel FCS file pre-processing}
\usage{
ParallelPreprocess(
  fnames,
  fpath_out,
  compensate = TRUE,
  spillover = NULL,
  transform = TRUE,
  tf_list = NULL,
  tf_match = FALSE,
  cores = parallel::detectCores() - 1,
  verbose = TRUE
)
}
\arguments{
\item{fnames}{string vector. Full paths to FCS files}

\item{fpath_out}{string. Path to directory where results should be saved. If
it does not exist, it will be created via \code{dir.create}}

\item{compensate}{logical. Whether to apply compensation of signal spillover
to the FCS files. Defaults to \code{TRUE}}

\item{spillover}{numeric matrix with column and row names. Spillover matrix
to use for compensation, unless one specified within each FCS file should be
used instead. Defaults to \code{NULL}}

\item{transform}{logical. Whether to apply signal transformation to the FCS
files. If set to \code{TRUE}, \code{tf_list} must be specified as well. Defaults to
\code{TRUE}}

\item{tf_list}{\code{flowCore::transformList}. A list of transformation
instructions per channel to apply to the FCS files. Defaults to \code{NULL}}

\item{tf_match}{logical. Whether to use string matching to align the channel
names specified in \code{tf_list} to those actually found in the FCS files.
Defaults to \code{FALSE}}

\item{cores}{integer. Number of CPU cores to use for multi-threading (at
least 2). Defaults to number of detectable cores minus 1}

\item{verbose}{logical. Whether to indicate progress. Defaults to \code{TRUE}}
}
\value{
nothing is returned
}
\description{
Applies pre-processing (compensation and/or transformation) to FCS files.
The pre-processed files are saved into a specified output directory.
}
\details{
This function uses parallelisation via a SNOW cluster for speed-up.

All FCS file names (even without the entire path) must be unique, so that
they can be saved in \code{fpath_out} under their original names.

For compensation, the spillover matrix can either be specified directly (as
\code{spillover}).
The spillover matrix (from which the compensation matrix is computed via
matrix inversion) must be square and have (identical) row and column names
(channel names).
If \code{spillover} is not specified, we attempt to extract it from each
corresponding FCS file.
In the latter case, the \verb{$SPILL}, \code{spillover}, and \verb{$SPILLOVER} slots of the
list retrieved via \code{flowCore::spillover} are searched (in that order) for a
valid spillover matrix to use.

For transformation, a \code{flowCore::transformList} object must be provided,
specifying transformation per channel.
If the channel names in that list are not exactly identical to those in the
FCS files (but clearly matchable to them), set \code{tf_match} to \code{TRUE} to
automatically resolve this using the \code{MatchChannels} function from this
package.

For flow cytometry data, we recommend using \code{flowCore::estimateLogicle} to
create initial specs for the transformation and then tuning the linearisation
width (\code{w}) parameter for each channel's transformation until satisfactory.
For mass cytometry data, we recommend using \code{flowCore::arcsinhTransform} with
a scale factor (\code{b}) of about 5.
}
