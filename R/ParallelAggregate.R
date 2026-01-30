
#' Parallel flowFrame aggregation
#'
#' Creates a `flowCore::flowFrame` or a matrix with named columns that
#' aggregates data from specified FCS files.
#' The size of this aggregate may be capped at an approximately fixed size.
#'
#' This function uses parallelisation via a SNOW cluster for speed-up.
#'
#' @param fnames string vector. Full paths to FCS files
#' @param N integer or `Inf` or `NULL`. Target aggregate size (>2), or `Inf` or
#' `NULL` to take all events. Defaults to `NULL`
#' @param batches integer vector or `NULL`. Batch numbers per FCS file. Defaults
#' to `NULL`
#' @param cols integer or string vector or `NULL`. Indices or names of
#' channels to use, or `NULL` to use all. Defaults to `NULL`
#' @param keep_file_order logical. Whether to keep the given order of `FCS`
#' files in the aggregate. Defaults to `FALSE`
#' @param keep_cell_order logical. Whether to keep the order of cells within
#' each sample in the aggregate. If `N` is `Inf` or `NULL`, implicitly `TRUE`
#' (due to no subsampling). Defaults to `FALSE`
#' @param seed integer or `NULL`. Seed for random number generation for
#' subsampling, or `NULL` to not set a seed. Defaults to `NULL`
#' @param cores integer. Number of CPU cores to use for multi-threading (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param as_flowFrame logical. Whether to return a `flowCore::flowFrame` rather
#' than an expression matrix with channel names for columns. Defaults to `TRUE`
#' @param guid string. Value of the `GUID` descriptor if
#' `flowCore::flowFrame` is to be returned. Defaults to *'Aggregate.fcs'*
#' @param fsom FlowSOM object or `NULL`. FlowSOM model with metaclustering with
#' channels matching to those in all FCS files, to only extract events that map
#' onto a specific metacluster, or `NULL` to not filter events that way.
#' Defaults to `NULL`
#' @param metacluster integer. Metacluster number to extract events, if FlowSOM
#' model is specified. Defaults to 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for `flowCore::read.FCS`
#'
#' @details
#' For a target aggregate size `N` and `n` FCS files, each `flowFrame`
#' contributes about `N/n`, unless one or more files have fewer than `N/n`
#' cells, in which case all their events are taken and the sample sizes for
#' other files is increased to compensate for this.
#'
#' If a corresponding FlowSOM model with metaclustering is given, we can choose
#' to extract only events that get mapped onto a specific metacluster.
#' In this case, up to `N/n` events per sample are picked greedily.
#'
#' The aggregate gets extra columns: *'File'* with index of FCS file of origin
#' per event and *'OriginalIndex'* with index of each event within its FCS file
#' of origin.
#'
#' If batch per file is specified, the *'Batch'* column is added with the batch
#' number of origin per event.
#'
#' The output gets attributes `fnames`, `seed`, and `metacluster` (which gets
#' `NULL` if no FlowSOM model was given) with the corresponding argument values.
#' These can be accessed as `attributes(output_object)$fnames`, etc.
#'
#' @return `flowCore::flowFrame` or matrix with named columns
#'
#' @export
ParallelAggregate <- function(
    fnames,
    N,
    batches         = NULL,
    cols            = NULL,
    keep_file_order = FALSE,
    keep_cell_order = FALSE,
    seed            = NULL,
    cores           = parallel::detectCores()-1,
    as_flowFrame    = TRUE,
    guid            = 'Aggregate.fcs',
    fsom            = NULL,
    metacluster     = 1,
    verbose         = TRUE,
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
  stopifnot('`N` must be NULL, Inf or integer greater than 2' =
              is.null(N) ||
              (is.infinite(N) && length(N)==1) ||
              (is.numeric(N) && N%%N==0 && length(N)==1 && N>2))
  stopifnot('If specified, `N` must be greater than length of `fnames`' =
              is.null(N) || is.infinite(N) || N>length(fnames))
  stopifnot('If specified, `batches` must be a non-list vector of integers of
             the same length as `fnames`' =
              is.null(batches) ||
              (is.vector(batches) && !is.list(batches)) &&
              (sum(batches%%batches)==0) && length(batches)==length(fnames))
  if (!is.null(cols)) {
    stopifnot('If specified, `cols` must be a vector of column names or indices' =
                (is.vector(cols) && length(cols)>0 &&
                   (is.numeric(cols)||is.character(cols))))
  }
  stopifnot('`keep_file_order` must be a single logical' =
              is.logical(keep_file_order) && length(keep_cell_order)==1)
  stopifnot('`keep_cell_order` must be a single logical' =
              is.logical(keep_cell_order) && length(keep_cell_order)==1)
  stopifnot('`cores` must be a single integer' =
              is.numeric(cores) && cores%%cores==0 && length(cores)==1)
  stopifnot('`cores` must be at least 2' = cores>=2)
  stopifnot('`as_flowFrame` must be a single logical' =
              is.logical(as_flowFrame) && length(as_flowFrame)==1)
  stopifnot('`guid` must be a single string' =
              is.character(guid) && length(guid)==1)
  stopifnot('`fsom` must be a single FlowSOM object or `NULL`' =
              is.null(fsom) || class(fsom)=='FlowSOM')
  stopifnot('If specified, `fsom` must have metaclustering' =
              is.null(fsom) || !is.null(fsom$metaclustering))
  stopifnot('`metacluster` must be a single positive integer' =
              is.numeric(metacluster) && metacluster%%metacluster==0 &&
              length(metacluster)==1 && metacluster>0)
  stopifnot('`metacluster` is out of bounds of `fsom` metaclustering' =
              is.null(fsom) || metacluster%in%levels(fsom$metaclustering))
  stopifnot('`verbose` must be a single logical' =
              is.logical(verbose) && length(verbose)==1)

  ## Init

  note_batch <- !is.null(batches)
  select_mc  <- !is.null(fsom)
  take_all   <- is.null(N) || is.infinite(N)
  n_files    <- length(fnames)

  if (verbose) {
    if (select_mc) {
      msg <- paste0('Aggregating metacluster ', metacluster,
                    ' expression data from ', n_files, ' FCS files')
    } else {
      msg <- paste0('Aggregating expression data from ', n_files, ' FCS files')
    }
    msg <- paste0(msg, '\nUsing ', cores, ' CPU cores')
    if (select_mc) {
      msg <- paste0(
        msg, '\nUsing FlowSOM model from ',
        strptime(fsom$info$date, format = '%Y-%m-%d %H:%M')
      )
    }
    message(msg)
  }

  ## Import multi-threading functions

  makeCluster    <- snow::makeCluster
  registerDoSNOW <- doSNOW::registerDoSNOW
  foreach        <- foreach::foreach
  `%dopar%`      <- foreach::`%dopar%`
  stopCluster    <- snow::stopCluster

  ## Determine sample size per file

  if (is.null(fsom)) {
    s_full    <- sapply( # FCS file sizes
      fnames, function(f) {
        as.integer(flowCore::read.FCSheader(f, keyword = '$TOT')[[1]][['$TOT']])
      })
    if (take_all) { # N is NULL or Inf? take all data
      s_perfile <- s_full
      N <- sum(s_perfile)
    } else { # N is finite? take N/n_files per file with adjustment
      s_generic <- ceiling(N/n_files)
      s_perfile <- rep(NA, times = n_files)
      mask_small <- s_full<s_generic
      mask_large <- !mask_small
      n_small    <- sum(mask_small)
      n_large    <- sum(mask_large)
      s_perfile[mask_small] <- s_full[mask_small]
      s_adjusted            <- ceiling((N-sum(s_perfile, na.rm = TRUE))/n_large)
      s_perfile[mask_large] <- s_adjusted
    }
  } else { # metacluster selected? take up to N/n_files greedily
    s_perfile <- rep(ceiling(N/n_files), times = n_files)
  }

  ## Set up SNOW cluster & progress bar

  cl <- makeCluster(cores)
  registerDoSNOW(cl)

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
    opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
  } else {
    opts <- list()
  }

  ## Aggregate expression data

  if (keep_file_order) { # maybe sort by file at the end
    f_sort <- function(res) res[order(res[, 'File']), , drop = FALSE]
  } else {
    f_sort <- NULL
  }
  if (take_all) { # maybe subsample
    f_cell_idcs <- function(n, s) seq_len(n)
  } else {
    if (keep_cell_order) { # maybe sort by cell index at each pass
      f_cell_idcs <- function(n, s) sort(sample.int(n, s))
    } else {
      f_cell_idcs <- function(n, s) sample.int(n, s)
    }
  }
  if (note_batch) { # maybe note batch of origin per even
    f_annotate <- function(d, i, idcs)
      cbind(d, 'File' = i, 'Batch' = batches[i], 'OriginalIndex' = idcs)
  } else {
    f_annotate <- function(d, i, idcs)
      cbind(d, 'File' = i, 'OriginalIndex' = idcs)
  }

  packages <- c('flowCore')
  if (select_mc) { # metacluster selected? use FlowSOM & model parameters
    select_mc <- TRUE
    packages <- c(packages, 'FlowSOM')
    codes <- fsom$map$codes
    meta <- fsom$metaclustering
  }

  res <- foreach( # iterate over FCS files asynchronously
    i             = seq_along(fnames),
    .combine      = rbind,
    .inorder      = FALSE,
    .final        = f_sort,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {

    if (!is.null(cols)&&length(cols)>0) {
      d <- flowCore::read.FCS(fnames[i], ...)[, cols, drop = FALSE]@exprs
    } else {
      d <- flowCore::read.FCS(fnames[i], ...)@exprs
    }
    if (select_mc) {
      mc_idcs <- which(
        meta[FlowSOM:::MapDataToCodes(codes, d)[, 1]]==metacluster
      )
      d <- d[mc_idcs, , drop = FALSE]
    }
    nd <- nrow(d)
    if (nd<2) {
      return(f_annotate(d, i, idcs))
    }
    s <- min(s_perfile[i], nd)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    idcs <- f_cell_idcs(nd, s)
    d <- d[idcs, , drop = FALSE]
    return(f_annotate(d, i, idcs))
  }

  if (verbose) {
    close(pb)
  }
  stopCluster(cl)
  invisible(gc())
  rm(cl)
  invisible(gc())

  ## Convert to flowFrame if requested

  if (as_flowFrame) {

    ff <- flowCore::read.FCS(fnames[1], which.lines = 1, ...)
    if (!is.null(cols) && length(cols)>0) {
      ff <- ff[, cols, drop = FALSE]
      ff@exprs <- res[, cols, drop = FALSE]
    } else {
      ff@exprs <- res[, flowCore::colnames(ff), drop = FALSE]
    }
    if (note_batch) {
      ff <- flowCore::fr_append_cols(ff, res[, c('File', 'Batch', 'OriginalIndex')])
    } else {
      ff <- flowCore::fr_append_cols(ff, res[, c('File', 'OriginalIndex')])
    }
    ff@description$GUID <- 'Aggregate.fcs'
    res <- ff
  }

  ## Add attributes for reproducibility

  attributes(res)[['fnames']]      <- fnames
  attributes(res)[['seed']]        <- seed
  attributes(res)[['metacluster']] <- if (select_mc) metacluster else NULL

  res
}
