
#' Match channel names
#'
#' Uses basic string matching to align two sets of channel names that should in
#' principle clearly map onto each other (same panel).
#' For example, if a *'Comp-'* prefix has been added to channel names stored in
#' an FCS file, the channel names specified in a `flowCore::transformList`
#' without the prefix are still clearly matchable to them.
#'
#' @param ref string vector. Tentative channel names to compare against the
#' ones in the FCS files (*e.g.*, from a `flowCore::transformList`)
#' @param real string vector. Channel names from the FCS files (retrievable
#' with `flowCore::colnames`)
#'
#' @details
#' Multiple matching schemes are used, to try and find a 1-to-1 matching.
#' The channels are matched first by their numeric part.
#' If this is insufficient, they are matched alphanumerically (by a letter
#' followed by numbers, as commonly used in naming channels).
#' If this is insufficient, an exact match is sought.
#' If all three of these approaches fail to yield a 1-to-1 matching, an error
#' message is generated.
#'
#' @return `data.frame` with columns *'Reference'*, *'Real'*, and *'Index'*
#' (mapping of `ref` onto `real`)
#'
#' @export
MatchChannels <- function(
    ref,
    real
) {
  
  matched_idcs <- rep(NA, length(ref))
  
  for (idx in seq_along(ref)) {
    
    one_ref <- ref[idx]
    
    ## Match by numeric part
    idcs <- which(
      stringr::str_extract(
        one_ref, '[[:digit:]]+'
      ) == stringr::str_extract(
        real, '[[:digit:]]+'
      )
    )
    
    if (length(idcs) != 1) {
      
      ## Match by alphanumeric part
      idcs <- which(
        stringr::str_extract(
          one_ref, '[A-Za-z][[0-9]]{3}'
        ) == stringr::str_extract(
          real,  '[A-Za-z][[0-9]]{3}'
        )
      )
    }
    
    if (length(idcs) != 1) {
      
      ## Match exactly
      idcs <- which(real == one_ref)
    }
    
    if (length(idcs) == 0) {
      
      idcs <- NA
    } else if (length(idcs) > 1) {
      
      stop(
        'transformList channel "', one_ref,
        '" was matched to multiple flowFrame channels:\n\t',
        paste(real[idcs], collapse = ',\n\r')
      )
    }
    
    matched_idcs[idx] <- idcs
  }
  
  matched <- real[matched_idcs]
  
  ## Aligned channel names
  
  data.frame(
    'Reference' = ref,
    'Real'      = matched,
    'Index'     = matched_idcs
  )
}

#' Parallel FCS file pre-processing
#'
#' Applies pre-processing (compensation, and/or transformation) to FCS files.
#' The pre-processed files are saved into a specified output directory.
#'
#' This function uses parallelisation via a SNOW cluster for speed-up.
#'
#' @param fnames string vector. Full paths to FCS files
#' @param fpath_out string. Path to directory where results should be saved. If
#' it does not exist, it will be created via `dir.create`
#' @param compensate logical. Whether to apply compensation of signal spillover
#' to the FCS files. Defaults to `TRUE`
#' @param spillover numeric matrix with column and row names. Spillover matrix
#' to use for compensation, unless one specified within each FCS file should be
#' used instead. Defaults to `NULL`
#' @param transform logical. Whether to apply signal transformation to the FCS
#' files. If set to `TRUE`, `tf_list` must be specified as well. Defaults to
#' `TRUE`
#' @param tf_list `flowCore::transformList`. A list of transformation
#' instructions per channel to apply to the FCS files. Defaults to `NULL`
#' @param tf_match logical. Whether to use string matching to align the channel
#' names specified in `tf_list` to those actually found in the FCS files.
#' Defaults to `FALSE`
#' @param cores integer. Number of CPU cores to use for multi-threading (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for `flowCore::read.FCS`
#'
#' @details
#' 
#' All FCS file names (even without the entire path) must be unique, so that
#' they can be saved in `fpath_out` under their original names.
#' 
#' For compensation, the spillover matrix can either be specified directly (as
#' `spillover`).
#' The spillover matrix (from which the compensation matrix is computed via
#' matrix inversion) must be square and have (identical) row and column names
#' (channel names).
#' If `spillover` is not specified, we attempt to extract it from each
#' corresponding FCS file.
#' In the latter case, the `$SPILL`, `spillover`, and `$SPILLOVER` slots of the
#' list retrieved via `flowCore::spillover` are searched (in that order) for a
#' valid spillover matrix to use.
#' 
#' For transformation, a `flowCore::transformList` object must be provided,
#' specifying transformation per channel.
#' If the channel names in that list are not exactly identical to those in the
#' FCS files (but clearly matchable to them), set `tf_match` to `TRUE` to
#' automatically resolve this using the `MatchChannels` function from this
#' package.
#' 
#' For flow cytometry data, we recommend using `flowCore::estimateLogicle` to
#' create initial specs for the transformation and then tuning the linearisation
#' width (`w`) parameter for each channel's transformation until satisfactory.
#' For mass cytometry data, we recommend using `flowCore::arcsinhTransform` with
#' a scale factor (`b`) of about 5.
#'
#' @return nothing is returned
#'
#' @export
ParallelPreprocess <- function(
    fnames,
    fpath_out,
    compensate = TRUE,
    spillover  = NULL,
    transform  = TRUE,
    tf_list    = NULL,
    tf_match   = FALSE,
    cores      = parallel::detectCores()-1,
    verbose    = TRUE,
    ...
) {

  ## Validate inputs

  stopifnot('`fnames` must be a non-empty non-list string vector' =
              is.vector(fnames) && is.character(fnames) && !is.list(fnames) &&
              length(fnames)>0)
  stopifnot('Some files in `fnames` are missing' = all(file.exists(fnames)))
  
  if (Sys.getenv('DUPLICATE_EXCEPTION')!='TRUE') {
    
    stopifnot('Some file paths in `fnames` are duplicates' =
                all(!duplicated(fnames)))
    stopifnot('Some file names (without paths) in `fnames` are duplicates' =
                all(!duplicated(basename(fnames))))
  }
  
  stopifnot('Some files in `fnames` are not FCS files' =
              tolower(substr(fnames, nchar(fnames)-3, nchar(fnames))) == '.fcs')
  stopifnot('`fpath_out` must be a single string' =
              is.character(fpath_out) && length(fpath_out)==1)
  if (!file.exists(fpath_out)) {
    res_fpath_out <- tryCatch(
      expr    = { dir.create(fpath_out) },
      error   = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!res_fpath_out) {
      stop('Cannot create the results directory "', fpath_out, '"')
    }
  }
  stopifnot('`compensate` must be a single logical' = 
              is.logical(compensate) && length(compensate)==1)
  if (!is.null(spillover)) {
    stopifnot(
        'If specified, `spillover` must be a square numeric matrix with (identical) row and column names' = 
        is.matrix(spillover) &&
        nrow(spillover)==ncol(spillover) &&
        !is.null(colnames(spillover)) &&
        !is.null(rownames(spillover)) &&
        all(colnames(spillover)==rownames(spillover))
    )
  }
  stopifnot('`transform` must be a single logical' = 
              is.logical(transform) && length(transform)==1)
  
  panel <- GetPanels(fnames[1])[[1]]
  
  if (!is.null(tf_list)) {
    stopifnot(
      'If specified, `tf_list` must be a flowCore::transformList' = 
       class(tf_list)=='transformList' &&
       attr(class(tf_list), 'package')=='flowCore'
    )
    stopifnot(
      'Some channels specified in `tf_list` were not found in the panel' =
       all(names(tf_list@transforms) %in% panel$Channel)
    )
  } else {
    if (isTRUE(transform)) {
      stop('If `transform` is set to TRUE, `tf_list` needs to be specified')
    }
  }
  stopifnot('`tf_match` must be a single logical' = 
              is.logical(tf_match) && length(tf_match)==1)
  stopifnot('`cores` must be a single integer' =
              is.numeric(cores) && cores%%cores==0 && length(cores)==1)
  stopifnot('`cores` must be at least 2' = cores>=2)
  stopifnot('`verbose` must be a single logical' =
              is.logical(verbose) && length(verbose)==1)
  
  ## Init
  
  n_fcs <- length(fnames)
  
  if (verbose) {
    message(paste0(
      'Applying ',
      paste(
        if (compensate) 'compensation' else NULL,
        if (transform) 'transformation' else NULL,
        sep = ' and '
      ),
      ' to ', n_fcs, ' FCS ',
      if (n_fcs>1) 'files' else 'file',
      '\nUsing ', cores, ' CPU cores'
    ))
  }
  
  ## Extract channels and markers from first file
  
  panel <- GetPanels(fnames[1])[[1]]
  
  ## Import multi-threading functions
  
  makeCluster    <- snow::makeCluster
  registerDoSNOW <- doSNOW::registerDoSNOW
  foreach        <- foreach::foreach
  `%dopar%`      <- foreach::`%dopar%`
  stopCluster    <- snow::stopCluster
  
  ## Set up SNOW cluster & progress bar
  
  clu <- makeCluster(cores)
  registerDoSNOW(clu)
  
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n_fcs, style = 3)
    opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
  } else {
    opts <- list()
  }

  ## Pre-process files in parallel

  packages <- c('flowCore')
  export   <- c('ValidateFCS')

  tmp <- foreach( # iterate over FCS files asynchronously
    i             = seq_along(fnames),
    .combine      = c,
    .export       = export,
    .inorder      = FALSE,
    .final        = NULL,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {

    ## Import FCS file
    
    fname <- fnames[i]
    ff <- flowCore::read.FCS(fname, ...)
    
    ## Resolve spillover matrix
    
    if (compensate) {
      if (is.null(spillover)) {
        
        spillover <- flowCore::spillover(ff)$`$SPILL`
        if (is.null(spillover)) {
          spillover <- flowCore::spillover(ff)$`spillover`
        }
        if (is.null(spillover)) {
          spillover <- flowCore::spillover(ff)$`$SPILLOVER`
        }
        if (is.null(spillover)) {
          stop(
            'Error at "', fname,
            '"\nSpillover matrix was not provided and could not be extracted ',
            'from the FCS file'
          )
        }
        colnames(spillover) <- rownames(spillover)
      }
      
      ch      <- colnames(spillover)
      missing <- !ch%in%flowCore::colnames(ff)
      if (any(missing)) {
        stop(
          'Error at "', fname,
          '"\nSample lacks channels specified in spillover matrix:\n\t',
          paste(ch[missing], collapse = '\n\t')
        )
      }
    }
    
    
    ## Resolve transformation instructions
    
    if (!is.null(tf_list)) {
      
      ch_ff <- flowCore::colnames(ff)
      ch_tf <- names(tf_list@transforms)
      
      if (tf_match) {
        
        m <- MatchChannels(ref = ch_tf, real = ch_ff)
        flowCore::colnames(ff)[m$Index] <- m$Reference
      }
      
      missing <- !ch_tf%in%ch_ff
      if (any(missing)) {
        stop(
          'Error at "', fname,
          '"\nSample lacks channels specified in transformation list:\n\t',
          paste(ch_tf[missing], collapse = '\n\t')
        )
      }
    }
    
    ## Compensate
    
    if (compensate) {
      
      ff <- flowCore::compensate(ff, spillover)
    }
    
    ## Transform
    
    if (transform) {
      
      ff <- flowCore::transform(ff, tf_list)
      if (tf_match) {
        
        flowCore::colnames(ff) <- ch_ff
      }
    }
    
    ## Validate and write FCS data
    
    flowCore::write.FCS(
      x        = ValidateFCS(ff@exprs, ff),
      filename = file.path(fpath_out, basename(fnames[i]))
    )
    
    NULL
  }

  if (verbose) {
    close(pb)
  }
  stopCluster(clu)
  invisible(gc())
  rm(clu)
  invisible(gc())

  return(invisible(NULL))
}
