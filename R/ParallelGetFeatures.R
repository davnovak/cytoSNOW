
#' Parallel FlowSOM feature extraction
#'
#' Computes matrices of abundance and state features per FCS sample.
#' This is done by mapping the expression matrices of selected FCS samples onto
#' an existing FlowSOM object.
#'
#' This function uses parallelization via a SNOW cluster for speed-up.
#'
#' @param fsom FlowSOM object. Trained FlowSOM model with channels matching 
#' those in all FCS files (same panel of markers)
#' @param fnames  vector. Full paths to FCS files
#' @param level  vector. Resolution level(s) of interest for feature
#' extraction. One or both of `c('clusters', 'metaclusters')`. Defaults to both
#' @param type  vector. Type(s) of features to extract. One or more of
#' `c('counts', 'proportions', 'medians', 'phenopositivities')`
#' State feature types (medians and phenopositivities) require `state_markers`
#' to be specified. Additionally, phenopositivities require `thresholds` to be
#' specified
#' @param state_markers  vector or NULL. Channels or markers to
#' consider for the extraction of state features. Required for computing medians
#' and/or phenopositivities. Defaults to NULL
#' @param thresholds named numeric vector or NULL. Expression values per channel
#' or marker that are the upper bounds for strictly negative phenotype. Required
#' for computing phenopositivities. Defaults to NULL
#' @param cores integer. Number of CPU cores to use for parallelization (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param sample_names string vector or NULL. Names of samples in `fnames` to
#' use as the rownames of features matrices. Otherwise, `fnames` themselves are
#' used. Defaults to NULL
#' @param verbose logical. Whether to indicate progress. Defaults to TRUE
#' @param ... optional additional named parameters for [flowCore::read.FCS()]
#'
#' @details
#' Two groups of features can be extracted: **abundance** features and **state**
#' features.
#' **Abundance** features can be used to describe compositional changes between
#' samples in terms of cell type representation.
#' Those can be tested using **differential abundance** modelling in *diffcyt*.
#' **State** features can be used to describe changes changes in cell state
#' between samples across different cell types.
#' Those can be tested using **differential state** modelling in *diffcyt* and
#' *iidx*.
#'
#' Possible **abundance** feature types are *'counts'* and *'proportions'*.
#' They correspond, respectively, to the absolute count and relative proportion
#' of cells per (meta)cluster in each sample.
#'
#' Possible **state** feature types are *'medians'* and *'phenopositivities'*.
#' Here, 'medians' correspond to the median expression values of a defined
#' subset of markers per (meta)cluster in each sample.
#' In contrast, *'phenopositivities'* correspond to the proportions of cells
#' that exceed a threshold for phenotypically positive for a defined subset of
#' markers, per (meta)cluster in each sample.
#' The phenopositivity thresholds then need to specified manually.
#'
#' The feature type names, in the order described here, correspond to
#' *'counts'*, *'percentages'*, *'MFIs'*, and *'percentages_positive'* in the
#' FlowSOM package.
#'
#' @return list with named elements corresponding to matrices per requested
#' feature `type`, each with prefix *'clusters_'* or *'metaclusters_'*
#'
#' @seealso [FlowSOM::GetFeatures()] for the original FlowSOM feature-extraction 
#' function
#'
#' @export
ParallelGetFeatures <- function(
    fsom,
    fnames,
    level           = c('clusters', 'metaclusters'),
    type            = c('counts', 'proportions', 'medians', 'phenopositivities'),
    state_markers   = NULL,
    thresholds      = NULL,
    cores           = parallel::detectCores()-1,
    sample_names    = NULL,
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
  
  stopifnot('`level` must be non-empty non-list string vector' =
              is.vector(level) && is.character(level) && !is.list(level) &&
              length(level)>0)
  level <- unique(as.vector(
    sapply(level, function(x) match.arg(x, c('clusters', 'metaclusters')))
  ))
  stopifnot('Invalid elements of `level` were found' = sum(is.na(level))==0)
  stopifnot('`fsom` must be a single FlowSOM object' = class(fsom)=='FlowSOM')
  stopifnot(
    'To get metacluster-level features, `fsom` must have metaclustering' =
      !('metaclusters'%in%level && is.null(fsom$metaclustering)))

  if (!is.null(state_markers) && length(state_markers)>0) {
    stopifnot(
      'If specified, `state_markers` must be a non-list string vector' =
        is.vector(state_markers) && is.character(state_markers) &&
        !is.list(state_markers) && length(state_markers)>0)
  }
  if (!is.null(thresholds) && length(thresholds)>0) {
    stopifnot(
      'If specified, `thresholds` must be a non-list named numeric vector' =
        is.numeric(thresholds) && !is.list(thresholds) &&
        length(thresholds)>0 && !is.null(names(thresholds)))
    if (is.null(state_markers) || length(state_markers)==0) {
      state_markers <- names(thresholds)
    }
  }
  if (!is.null(sample_names) && length(sample_names)>0) {
    stopifnot(
      'If specified, `sample_names` must be a non-list string vector' =
        is.vector(sample_names) && is.character(sample_names) &&
        !is.list(sample_names) && length(sample_names)>0)
    stopifnot(
      'If specified, `sample_names` must be of the same length as `fnames' =
        length(sample_names)==length(fnames)
    )
    
    if (Sys.getenv('DUPLICATE_EXCEPTION')!='TRUE') {
      
      stopifnot(
        'If specified, `sample_names` must not contain duplicates' =
          length(unique(sample_names))==length(sample_names)
      )
    }
  } else {
    sample_names <- fnames
  }
  stopifnot('`type` must be non-empty non-list string vector' =
              is.vector(type) && is.character(type) && !is.list(type) &&
              length(type)>0)
  type <- unique(as.vector(
    sapply(type, function(x) match.arg(x,
    c('counts', 'proportions', 'medians', 'phenopositivities')))
  ))
  stopifnot('Invalid elements of `type` were found' = sum(is.na(type))==0)
  stopifnot('`cores` must be a single integer' =
              is.numeric(cores) && cores%%cores==0 && length(cores)==1)
  stopifnot('`cores` must be at least 2' = cores>=2)
  stopifnot('`verbose` must be a single logical' =
              is.logical(verbose) && length(verbose)==1)

  fsom_channels <- as.vector(colnames(fsom$map$codes))
  fsom_markers  <- as.vector(FlowSOM::GetMarkers(fsom, fsom_channels))
  
  ff <- flowCore::read.FCS(fnames[1], ...)
  channels <- flowCore::colnames(ff)
  markers  <- flowCore::markernames(ff)[channels]

  if (!is.null(state_markers)) {
    state_map <- pmatch(state_markers, markers)
    nas <- is.na(state_map)
    if (any(nas)) {
      state_map[nas] <- pmatch(state_markers[nas], channels)
    }
    stopifnot('Some `state_markers` could not be mapped to channels/markers' =
                sum(is.na(state_map))==0)
    state_markers  <- markers[state_map]
    state_channels <- channels[state_map]
  }

  if (!is.null(thresholds)) {
    thr_map <- pmatch(names(thresholds), state_markers)
    nas <- is.na(thr_map)
    if (any(nas)) {
      thr_map[nas] <- pmatch(names(thresholds)[nas], state_channels)
    }
    stopifnot('Some names of `thresholds` could not be mapped to state channels/markers' =
                sum(is.na(thr_map))==0)
    
    thresholds  <- thresholds[thr_map]
    thr_markers <- `names<-`(thresholds, state_markers[thr_map])
    thr_channels <- `names<-`(thresholds, state_channels[thr_map])
  }

  ## Resolve resolution (cluster/metacluster)

  cl      <- 'clusters'%in%level
  mcl     <- 'metaclusters'%in%level
  n_files <- length(fnames)
  n_cl    <- FlowSOM::NClusters(fsom)
  if (mcl) {
    n_mcl   <- FlowSOM::NMetaclusters(fsom)
  }
  codes   <- fsom$map$codes
  meta <- NULL
  if (mcl) {
    meta <- fsom$metaclustering
  }

  ## Resolve feature types and feature matrix column names per type

  abund <- any(type %in% c('counts', 'proportions'))
  state <- any(type %in% c('medians', 'phenopositivities'))

  abund_cl_comps  <- if (abund && cl) paste0('C', seq_len(n_cl)) else NULL
  abund_mcl_comps <- if (abund && mcl) paste0('MC', seq_len(n_mcl)) else NULL

  meds_cl_comps  <- NULL
  meds_mcl_comps <- NULL
  if ('medians'%in%type) {
    meds_cl_comps <-
      if (cl) {
        apply(
          X = expand.grid(state_markers, seq_len(n_cl))[, c(2, 1)],
          MARGIN = 1, FUN = function(x) paste0('C', as.numeric(x[1]), ' ', x[2])
        )
      } else NULL
    meds_mcl_comps <-
      if (mcl) {
        apply(
          X = expand.grid(state_markers, seq_len(n_mcl))[, c(2, 1)],
          MARGIN = 1, FUN = function(x) paste0('MC', as.numeric(x[1]), ' ', x[2])
        )
      } else NULL
  }
  if ('phenopositivities'%in%type) {
    pheno_cl_comps <-
      if (cl) {
        apply(
          X = expand.grid(names(thr_markers), seq_len(n_cl))[, c(2, 1)],
          MARGIN = 1, FUN = function(x) paste0('C', as.numeric(x[1]), ' ', x[2])
        )
      } else NULL
    pheno_mcl_comps <-
      if (mcl) {
        apply(
          X = expand.grid(names(thr_markers), seq_len(n_mcl))[, c(2, 1)],
          MARGIN = 1, FUN = function(x) paste0('MC', as.numeric(x[1]), ' ', x[2])
        )
      } else NULL
  }

  if (verbose) {
    msg <- paste0(
      'Extracting FlowSOM features for ', n_files, ' FCS files using a FlowSOM',
      ' model from ', strptime(fsom$info$date, format = '%Y-%m-%d %H:%M'), '\n',
      'Resolution level: ', paste(level, collapse = ', '), '\n',
      'Feature types: ', paste(type, collapse = ', '), '\n',
      'Using ', cores, ' CPU cores'
    )
    message(msg)
  }

  ## Define function to combine results on the go

  f_combine <- function(res1, res2) {

    `names<-`(
      lapply(names(res1), function(x) rbind(res1[[x]], res2[[x]])),
      names(res1)
    )
  }

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
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
    opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
  } else {
    opts <- list()
  }

  ## Extract features in parallel

  packages <- c('flowCore', 'FlowSOM')

  res <- foreach( # iterate over FCS files asynchronously
    i             = seq_along(fnames),
    .combine      = f_combine,
    .inorder      = FALSE,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {

    fname <- fnames[i]
    expr <- flowCore::read.FCS(fname, ...)[, channels, drop = FALSE]@exprs

    ## Get (meta)cluster mappings

    cl_map <- factor(FlowSOM:::MapDataToCodes(
      codes = codes, newdata = expr[, fsom_channels, drop = FALSE], distf = 2
    )[, 1], levels = seq_len(n_cl))

    if (mcl) {
      mcl_map <- meta[cl_map]
    } else {
      mcl_map <- NULL
    }
    
    ## Initialise results
    
    cl_counts  <- NULL
    cl_props   <- NULL
    cl_meds    <- NULL
    cl_pheno   <- NULL
    mcl_counts <- NULL
    mcl_props  <- NULL
    mcl_meds   <- NULL
    mcl_phen   <- NULL
    
    ## Compute abundances

    if (abund) {

      cl_counts  <-
        if (cl) `rownames<-`(matrix(table(cl_map), nrow = 1), fname) else NULL
      mcl_counts <-
        if (mcl) `rownames<-`(matrix(table(mcl_map), nrow = 1), fname) else NULL

      if ('proportions'%in%type) {

        cl_props  <- if (cl) prop.table(cl_counts) else NULL
        mcl_props <- if (mcl) prop.table(mcl_counts) else NULL
      }
      
      if (!'counts'%in%type) {
        
        cl_counts  <- NULL
        mcl_counts <- NULL
      }
    }

    ## Compute states

    if (state) {

        cl_meds <-
          if ('medians'%in%type && cl) {
            `rownames<-`(
              matrix(unlist(
                lapply(
                  seq_len(n_cl),
                  function(i) {
                    d <- expr[cl_map==i, state_channels, drop = FALSE]
                    if (nrow(d)==0) {
                      return(rep(NA, length(state_channels)))
                    }
                    apply(d, 2, stats::median)
                  }
                )
              ), nrow = 1),
              fname
            )
          } else NULL

        mcl_meds <-
          if ('medians'%in%type && mcl) {
            `rownames<-`(
              matrix(unlist(
                lapply(
                  seq_len(n_mcl),
                  function(i) {
                    d <- expr[mcl_map==i, state_channels, drop = FALSE]
                    if (nrow(d)==0) {
                      return(rep(NA, length(state_channels)))
                    }
                    apply(d, 2, median)
                  }
                )
              ), nrow = 1),
              fname
            )
          } else NULL

        cl_pheno <- if ('phenopositivities'%in%type && cl) {
          `rownames<-`(
            matrix(unlist(
              lapply(
                seq_len(n_cl),
                function(i) {
                  d <- expr[cl_map==i, state_channels, drop = FALSE]
                  if (nrow(d)==0) {
                    return(rep(NA, length(state_channels)))
                  }
                  apply(
                    X = rbind(
                      expr[cl_map==i, names(thr_channels), drop = FALSE],
                      thr_channels
                    ),
                    MARGIN = 2,
                    FUN = function(x)
                      sum(utils::head(x, -1)>utils::tail(x, 1))/(length(x)-1)
                  )
                }
              )
            ), nrow = 1),
            fname
          )
        } else NULL

        mcl_pheno <- if ('phenopositivities'%in%type && mcl) {
          `rownames<-`(
            matrix(unlist(
              lapply(
                seq_len(n_mcl),
                function(i) {
                  d <- expr[mcl_map==i, state_channels, drop = FALSE]
                  if (nrow(d)==0) {
                    return(rep(NA, length(state_channels)))
                  }
                  apply(
                    X = rbind(
                      expr[mcl_map==i, names(thr_channels), drop = FALSE],
                      thr_channels
                    ),
                    MARGIN = 2,
                    FUN = function(x)
                      sum(utils::head(x, -1)>utils::tail(x, 1))/(length(x)-1)
                  )
                }
              )
            ), nrow = 1),
            fname
          )
        } else NULL
    }

    return(
      list(
        'clusters_counts'                = cl_counts,
        'clusters_proportion'            = cl_props,
        'clusters_medians'               = cl_meds,
        'clusters_phenopositivities'     = cl_pheno,
        'metaclusters_counts'            = mcl_counts,
        'metaclusters_proportions'       = mcl_props,
        'metaclusters_medians'           = mcl_meds,
        'metaclusters_phenopositivities' = mcl_pheno
      )
    )
  }

  if (verbose) {
    close(pb)
  }
  stopCluster(clu)
  invisible(gc())
  rm(clu)
  invisible(gc())

  ## Resolve row order and names

  res <- `names<-`(
    lapply(
      names(res),
      function(x) {
        if (is.null(res[[x]])) { return(NULL) }
        `rownames<-`(res[[x]][match(rownames(res[[x]]), fnames), ], sample_names)
      }
    ),
    names(res)
  )

  ## Resolve column names

  if (!is.null(res$clusters_counts)) {
    colnames(res$clusters_counts) <- abund_cl_comps
  }
  if (!is.null(res$metaclusters_counts)) {
    colnames(res$metaclusters_counts) <- abund_mcl_comps
  }
  if (!is.null(res$clusters_proportions)) {
    colnames(res$clusters_proportions) <- abund_mcl_comps
  }
  if (!is.null(res$metaclusters_proportions)) {
    colnames(res$metaclusters_proportions) <- abund_mcl_comps
  }
  if (!is.null(res$clusters_medians)) {
    colnames(res$clusters_medians) <- meds_cl_comps
  }
  if (!is.null(res$metaclusters_medians)) {
    colnames(res$metaclusters_medians) <- meds_mcl_comps
  }
  if (!is.null(res$clusters_phenopositivities)) {
    colnames(res$clusters_phenopositivities) <- pheno_cl_comps
  }
  if (!is.null(res$metaclusters_phenopositivities)) {
    colnames(res$metaclusters_phenopositivities) <- pheno_mcl_comps
  }

  res
}



