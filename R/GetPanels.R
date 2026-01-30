
#' Extract channels and markers from FCS files
#'
#' Retrieves channel and marker indices and names per FCS sample without parsing
#' the files themselves.
#' By only reading the FCS file headers, this offers a speed-up and allows to
#' check for any inconsistencies between sample panels to prevent errors or
#' confusion downstream.
#'
#' @param fnames string vector. Full paths to FCS files
#' @param ... optional additional named parameters for 
#' `flowCore::read.FCSheader`
#'
#' @return named `list` per file of `data.frame`s with columns *'Index'*,
#' *'Channel'*, and *'Marker'*
#'
#' @export
GetPanels <- function(
  fnames,
  ...
) {

  ## Validate inputs

  stopifnot('`fnames` must be a non-empty non-list string vector' =
              is.vector(fnames) && is.character(fnames) && !is.list(fnames) &&
              length(fnames)>0)
  stopifnot('Some files in `fnames` are missing' = all(file.exists(fnames)))
  
  if (Sys.getenv('DUPLICATE_EXCEPTION')!='TRUE') {
    
    stopifnot('Some files in `fnames` are duplicates' = all(!duplicated(fnames)))
    stopifnot('Some files in `fnames` are not FCS files' =
                tolower(substr(fnames, nchar(fnames)-3, nchar(fnames))) == '.fcs')
  }

  headers <- flowCore::read.FCSheader(fnames, ...)
  res <- lapply(headers, function(h) {
    channels <- h[grepl('\\$P[0-9]+N$', names(h))]
    markers  <- h[grepl('\\$P[0-9]+S$', names(h))]

    names(channels) <-
      substr(names(channels), start = 3, stop = nchar(names(channels))-1)
    names(markers) <-
      substr(names(markers), start = 3, stop = nchar(names(markers))-1)

    d <- data.frame(
      'Index' = as.numeric(names(channels)),
      'Channel' = as.vector(channels),
      'Marker' = NA
    )
    d$Marker[match(names(markers), d$Index)] <- as.vector(markers)

    d[order(d$Index), ]
  })
  names(res) <- fnames

  res
}
