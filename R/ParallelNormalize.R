
#' Split cytometry samples by metaclusters
#'
#' Splits FCS files into metaclusters defined by a trained FlowSOM model. 
#'
#' This function uses parallelization via a SNOW cluster for speed-up.
#' 
#' This is an internal cytoSNOW function that does not need to be called by the 
#' user. Input validation is not carried out, unlike in the (safer) user-level 
#' functions [cytoSNOW::ParallelNormalize.Train()] and 
#' [cytoSNOW::ParallelNormalize.apply()].
#'
#' @param fnames vector. Full paths to FCS files that should be split
#' @param batches vector. Batch labels per file in `fnames`
#' @param fsom FlowSOM object. Trained FlowSOM model with channels matching 
#' those in all FCS files (same panel of markers) to use for distinguishing 
#' cell compartments (metaclusters) to normalize separately
#' @param fpath_out string. Path to directory where the split FCS files should 
#' be saved
#' @param cols integer or string vector or `NULL`. Indices or names of channels 
#' which should be normalized, or `NULL` to use all. Defaults to `NULL`
#' @param remove_na logical. Whether empty sample-metacluster combinations 
#' should be left out of the returned data frame. Defaults to `TRUE`
#' @param cores integer. Number of CPU cores to use for parallelization (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for [flowCore::read.FCS()]
#'
#' @details
#' The original FCS file names are used, along with the postfix 
#' *'_metaclusterXX'* where *XX* stands for the metacluster number.
#'
#' @return data frame with columns *'Metacluster'*, *'FileName'*, *'Sample'* 
#' (*i.e.*, the original file name), and *'Batch'* that identifies the split FCS 
#' files written to disk
#'
#' @seealso [cytoSNOW::ParallelNormalize.Train()] for parallel training of 
#' CytoNorm normalization models and [cytoSNOW::ParallelNormalize.apply()] for 
#' the parallel application of CytoNorm normalization models
#' 
#' @export
.SplitSamplesByMetaclusters <- function(
  fnames,
  batches,
  fsom,
  fpath_out,
  cols      = NULL,
  remove_na = TRUE,
  cores     = parallel::detectCores()-1,
  verbose   = TRUE,
  ...
) {
  
  ## Import multi-threading functions
  
  makeCluster    <- snow::makeCluster
  registerDoSNOW <- doSNOW::registerDoSNOW
  foreach        <- foreach::foreach
  `%dopar%`      <- foreach::`%dopar%`
  stopCluster    <- snow::stopCluster
  
  ## Split samples by metaclusters in parallel
  
  codes <- fsom$map$codes
  if (!is.null(cols) && length(cols)>1) {
    codes <- codes[, cols, drop = FALSE]
  }
  meta  <- fsom$metaclustering
  umcl  <- unique(meta)
  nmcl  <- length(umcl)
  
  packages <- c('FlowSOM')
  
  clu <- makeCluster(cores)
  registerDoSNOW(clu)
  
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = length(fnames), style = 3)
    opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
  } else {
    opts <- list()
  }
  
  res <- foreach( # iterate over FCS files asynchronously
    i             = seq_along(fnames),
    .combine      = rbind,
    .inorder      = FALSE,
    .final        = NULL,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {
    
    fname <- fnames[i]
    batch <- batches[i]
    
    ff <- flowCore::read.FCS(fname, ...)
    d <- ff@exprs
    
    if (!is.null(cols) && length(cols)>0) {
      d <- d[, cols, drop = FALSE]
    } else {
      d <- d[, intersect(colnames(d), colnames(codes)), drop = FALSE]
    }
    
    mcls <- # metacluster label per cell
      meta[FlowSOM:::MapDataToCodes(
        codes   = codes[, intersect(colnames(codes), colnames(d))],
        newdata = d
      )[, 1]]
    
    fnames_split <- c()
    
    for (mcl in umcl) {
      
      mask <- mcls==mcl
      if (sum(mask)>0) {
        
        fname_split <- tempfile(
          pattern = paste0(
            gsub('[.]fcs$', '', gsub('[:/]', '_', fname)),
            '_MC', mcl, '_'
          ),
          tmpdir = fpath_out,
          fileext = '.fcs'
        )
        fnames_split <- c(fnames_split, fname_split)
        suppressWarnings(
          flowCore::write.FCS(
            ff[mask, ],
            file = fname_split
          )
        )
      } else {
        
        fnames_split <- c(fnames_split, NA)
      }
    }
    
    meta_df <- data.frame(
      'Metacluster' = umcl,
      'FileName'    = fnames_split,
      'Sample'      = fname,
      'Batch'       = batch
    )
    
    if (remove_na) {
      meta_df <- meta_df[!is.na(meta_df$FileName), ]
    }
    
    return(
      meta_df
    )
  }
  
  if (verbose) {
    close(pb)
  }
  stopCluster(clu)
  invisible(gc())
  rm(clu)
  invisible(gc())
  
  res
}

#' Split cytometry samples by gates
#'
#' Splits FCS files into gates defined by a gating function. 
#'
#' This function uses parallelization via a SNOW cluster for speed-up.
#' 
#' This is an internal cytoSNOW function that does not need to be called by the 
#' user. Input validation is not carried out, unlike in the (safer) user-level 
#' functions [cytoSNOW::ParallelNormalize.Train()] and 
#' [cytoSNOW::ParallelNormalize.apply()].
#'
#' @param fnames vector. Full paths to FCS files that should be split
#' @param batches vector. Batch labels per file in `fnames`
#' @param gating_fun function. Takes a [flowCore::flowFrame] object as input and
#' returns a gate label per event
#' @param gate_names vector. Names of all gates from `gating_fun`
#' @param fpath_out string. Path to directory where the split FCS files should 
#' be saved
#' @param cols integer or string vector or `NULL`. Indices or names of channels 
#' which should be normalized, or `NULL` to use all. Defaults to `NULL`
#' @param remove_na logical. Whether empty sample-gate combinations should be 
#' left out of the returned data frame. Defaults to `TRUE`
#' @param cores integer. Number of CPU cores to use for parallelization (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for [flowCore::read.FCS()]
#'
#' @details
#' The original FCS file names are used, along with the postfix 
#' *'_GateXX'* where *XX* stands for an internally generated gate number.
#'
#' @return data frame with columns *'Gate'*, *'FileName'*, *'Sample'* 
#' (*i.e.*, the original file name), and *'Batch'* that identifies the split FCS 
#' files written to disk
#'
#' @seealso [cytoSNOW::ParallelNormalize.Train()] for parallel training of 
#' CytoNorm normalization models and [cytoSNOW::ParallelNormalize.apply()] for 
#' the parallel application of CytoNorm normalization models
#' 
#' @export
.SplitSamplesByGates <- function(
  fnames,
  batches,
  gating_fun,
  gate_names,
  fpath_out,
  cols      = NULL,
  remove_na = TRUE,
  cores     = parallel::detectCores()-1,
  verbose   = TRUE,
  ...
) {
  
  ## Import multi-threading functions

  makeCluster    <- snow::makeCluster
  registerDoSNOW <- doSNOW::registerDoSNOW
  foreach        <- foreach::foreach
  `%dopar%`      <- foreach::`%dopar%`
  stopCluster    <- snow::stopCluster

  ## Create gate codes
  
  gate_codes <- `names<-`(sprintf('%03d', seq_along(gate_names)), gate_names)
  
  ## Split samples by gates in parallel
  
  clu <- makeCluster(cores)
  registerDoSNOW(clu)
  
  packages <- c()
  
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = length(fnames), style = 3)
    opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
  } else {
    opts <- list()
  }

  res <- foreach( # iterate over FCS files asynchronously
    i             = seq_along(fnames),
    .combine      = rbind,
    .inorder      = FALSE,
    .final        = NULL,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {

    fname <- fnames[i]
    batch <- batches[i]

    ff <- flowCore::read.FCS(fname, ...)
    gates <- gating_fun(ff)
    
    invalid_gates <- unique(gates[!gates%in%names(gate_codes)])
    
    if (length(invalid_gates)>0) {
      stop(
        'The following gates from `gating_fun` ',
        'were not found in `gate_names`: ',
        paste0(invalid_gates)
      )
    }

    fnames_split <- c()

    for (i in seq_along(gate_codes)) {

      gate <- names(gate_codes)[i]
      code <- gate_codes[i]
      mask <- gates==gate
      
      if (sum(mask)>0) {

        fname_split <- tempfile(
          pattern = paste0(
            gsub('[.]fcs$', '', gsub('[:/]', '_', fname)),
            '_Gate', code, '_'
          ),
          tmpdir = fpath_out,
          fileext = '.fcs'
        )
        fnames_split <- c(fnames_split, fname_split)
        
        ff_out <- ff[mask, ]
        if (!is.null(cols) && length(cols)>0) {
          ff_out <- ff_out[, cols]
        }
        
        suppressWarnings(
          flowCore::write.FCS(
            ff_out,
            file = fname_split
          )
        )
      } else {

        fnames_split <- c(fnames_split, NA)
      }
    }

    meta_df <- data.frame(
      'Gate'        = names(gate_codes),
      'FileName'    = fnames_split,
      'Sample'      = fname,
      'Batch'       = batch
    )
    if (remove_na) {
      meta_df <- meta_df[!is.na(meta_df$FileName), ]
    }

    return(
      meta_df
    )
  }

  if (verbose) {
    close(pb)
  }
  stopCluster(clu)
  invisible(gc())
  rm(clu)
  invisible(gc())

  res
}

#' Parallel training of normalization model
#'
#' Trains a batch-effect correction model using CytoNorm. This model can later 
#' be applied using [cytoSNOW::ParallelNormalize.Apply()]. Instead of using 
#' FlowSOM clustering for splitting files into compartments, a custom gating 
#' function can be used (see the `gating_fun` and `gate_names` parameters and 
#' the *Details* section).
#'
#' @param fnames vector. Full paths to training FCS files
#' @param batches vector. Batch labels per file in `fnames`
#' @param cols integer or string vector or `NULL`. Indices or names of channels 
#' which should be normalized, or `NULL` to use all. Defaults to `NULL`
#' @param fsom FlowSOM object. Trained FlowSOM model with channels matching 
#' those in all FCS files (same panel of markers) to use for distinguishing 
#' cell compartments (metaclusters) to normalize separately. Or `NULL` to train 
#' a new one using `fsom_params` as the training parameters, or to use a gating 
#' model instead of FlowSOM. Defaults to `NULL`
#' @param fsom_params named list. Parameters to pass to 
#' [CytoNorm::prepareFlowSOM()] for the FlowSOM model training. Ignored if 
#' `fsom` is specified or using a gating model instead. If specified, 
#' `fsom_params` must give at least the *'nCells'* specified. See *Details* 
#' for defaults
#' @param gating_fun function. Takes a [flowCore::flowFrame] object as input and
#' returns a gate label per event. Or `NULL` to use FlowSOM for splitting files
#' instead. Defaults to `NULL`
#' @param gate_names vector. Names of all gates from `gating_fun`. Or `NULL` to 
#' use FlowSOM for splitting files. Must be specified if `gating_fun` is 
#' specified. Defaults to `NULL`
#' @param normMethod.train training function of the signal normalization method 
#' to use. Defaults to `CytoNorm::QuantileNorm.train`
#' @param norm_params named list. Parameters to pass to the normalization method 
#' function specified in `norm_method`. Defaults to `list('nQ' = 99)`
#' @param seed integer or `NULL`. Seed for reproducibility
#' @param cores integer. Number of CPU cores to use for parallelization (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for [flowCore::read.FCS()]
#'
#' @details
#' This function uses parallelization via a SNOW cluster for speed-up.
#' 
#' If no precomputed FlowSOM model is given (`fsom` is not specified), the 
#' default training parameters are `nCells = 1e6`, `xdim = 15`, `ydim = 15`,
#' `nClus = 30`, and `scale = FALSE`. If `colsToUse` (channels to use for 
#' clustering) are not specified here, the value of `cols` (all channels that 
#' should be normalized) is used. If `cols` is also not specified, all channels 
#' in the panel are used. It is strongly recommended to specify both the `cols` 
#' parameter, and the `colsToUse` slot of `fsom_params`.
#' 
#' A gating model can be used to split FCS files instead of a FlowSOM model. To 
#' do this, two input parameters need to be specified. First, a gating function 
#' that takes a [flowCore::flowFrame] as input and returns a vector of gate 
#' labels as output (`gating_fun`). Second, a vector with the names of all gates
#' that `gating_fun` can output (`gate_names`). If both are specified, the 
#' FlowSOM model and/or parameters are ignored. Instead of splitting training 
#' FCS files by metacluster to train a separate normalization function for each 
#' of them, they are split by gate, as defined by the user's function. This 
#' gives more control and ensures that the model doesn't skew genuine 
#' differential abundances of known cell types across batches.
#' 
#' @note
#' The original [CytoNorm::CytoNorm.train()] function includes a `transformList` 
#' parameter that lets the user apply transformation in conjunction with the 
#' normalization model training/application. In contrast, this function does not
#' include that parameter. This is to encourage doing preprocessing separately 
#' first, checking whether it's working correctly, and only then proceeding with
#' normalization.
#'
#' Note that this function generates temporary files in the [tempdir()] 
#' directory, which it then automatically removes when done.
#'
#' If a FlowSOM model (not gating) is used for the splitting of files, the 
#' trained normalization model is compatible with the original 
#' [CytoNorm::CytoNorm.normalize()] function for applying the normalization. 
#' However, we recommend the use of [cytoSNOW::ParallelNormalize.Apply()].
#'
#' @return list with named elements for the CytoNorm normalization model 
#' (*'clusterRes'*) and the splitting model: either FlowSOM (*'fsom'*) or the 
#' custom gating model (*'gating'*).
#'
#' @seealso [CytoNorm::CytoNorm.train()] for the original CytoNorm training
#' function
#' @seealso [cytoSNOW::ParallelNormalize.Apply()] for applying the trained 
#' normalization model
#' 
#' @export
ParallelNormalize.Train <- ParallelNormalise.Train <- function(
  fnames,
  batches,
  cols             = NULL,
  fsom             = NULL,
  fsom_params      = list(
    'nCells' = 1e6,
    'xdim'   = 15,
    'ydim'   = 15,
    'nClus'  = 30,
    'scale'  = FALSE
  ),
  gating_fun       = NULL,
  gate_names       = NULL,
  normMethod.train = CytoNorm::QuantileNorm.train,
  norm_params      = list(
    'nQ' = 99
  ),
  seed             = NULL,
  cores            = parallel::detectCores()-1,
  verbose          = TRUE,
  ...
) {
  
  ## Validate inputs
  
  stopifnot('`fnames` must be a non-empty non-list string vector' =
              is.vector(fnames) && is.character(fnames) && !is.list(fnames) &&
              length(fnames)>0)
  stopifnot('Some files in `fnames` are missing' = all(file.exists(fnames)))
  
  stopifnot('`batches` must be a non-empty non-list string or numeric vector' =
              is.vector(batches) &&
              (is.character(batches) || is.numeric(batches)) &&
              !is.list(batches) && length(batches)>0)
  stopifnot(
    '`batches` must be of the same length as `fnames' =
      length(batches)==length(fnames)
  )
  if (Sys.getenv('DUPLICATE_EXCEPTION')!='TRUE') {
    
    stopifnot(
      '`fnames` must not contain duplicates' =
        length(unique(fnames))==length(fnames)
    )
  }
  if (!is.null(cols)) {
    stopifnot('If specified, `cols` must be a vector of column names or indices' =
                (is.vector(cols) && length(cols)>0 &&
                   (is.numeric(cols)||is.character(cols))))
  }
  stopifnot('If specified, `fsom` must be a single FlowSOM object' =
              is.null(fsom) || class(fsom)=='FlowSOM')
  stopifnot(
    'If specified, `fsom_params` must be a non-duplicitly named or empty list' =
      is.list(fsom_params) && (length(fsom_params)==0 || (
        !is.null(names(fsom_params)) &&
          length(unique(names(fsom_params)))==length(fsom_params))
    ))
  
  if (!is.null(gating_fun)) {
    stopifnot(
      'If `gating_fun` is specified, `gate_names` must be specified as well' =
        !is.null(gate_names)
    )
    stopifnot(
      'If specified, `gating_fun` must be a function' = is.function(gating_fun)
    )
  }
  if (!is.null(gate_names)) {
    stopifnot(
      'If `gate_names` is specified, `gating_fun` must be specified as well' =
        !is.null(gating_fun)
    )
    stopifnot(
      'If specified, `gate_names` must be a non-empty non-list vector of strings numbers or factors' =
      is.vector(gate_names) &&
        (is.character(gate_names) ||
           is.factor(gate_names) ||
           is.numeric(gate_names)) &&
        length(gate_names)>0 &&
        !is.list(gate_names)
    )
  }
  
  if (length(fsom_params)>0) {
    stopifnot('If specified, `fsom_params` must give `nCells` as a parameter' =
                'nCells'%in%names(fsom_params))
  }
  stopifnot('`normMethod.train` must be a function`' = 
              is.function(normMethod.train))
  stopifnot(
    'If specified, `norm_params` must be a non-duplicitly named or empty list' =
      is.list(norm_params) && (length(norm_params)==0 || (
        !is.null(names(norm_params)) &&
          length(unique(names(norm_params)))==length(norm_params))
      ))
  
  if ('nQ' %in% names(norm_params)) {
    stopifnot(
      'Quantile count in goal distribution and CytoNorm call is different' =
        nrow(norm_params[['goal']][[1]])==norm_params[['nQ']]
    )
  }
  
  ## Determine a location for temporary files
  
  tmp <- tempdir(check = TRUE)
  
  ## Init
  
  n_fcs      <- length(fnames)
  n_batch    <- length(unique(batches))
  panel      <- GetPanels(fnames[1], ...)[[1]]
  use_gating <- !is.null(gating_fun) && !is.null(gate_names)
  if (use_gating) {
    gate_names <- sort(as.character(gate_names))
  }
  
  n_cols  <-
    if (!is.null(cols)) {
      length(cols)
    } else {
      nrow(panel)
    }
  
  if (verbose) {
    message(paste0(
      'Training CytoNorm normalization model for batch-effect correction',
      '\n• Number of normalizable channels: ', n_cols,
      '\n• Input data: ', n_fcs, ' FCS ', if (n_fcs>1) 'files' else 'file',
      ' across ', n_batch, ' experimental ',
      if (n_batch>1) 'batches' else 'batch',
      '\n• Separate normalization functions will be trained per ',
      if (!use_gating) { 'FlowSOM metacluster' } else { 'gate' },
      '\n• Using ', cores, ' CPU cores'
    ))
  }
  
  if (!use_gating) {
    
    ## Obtain a FlowSOM model
    
    if (is.null(fsom)) {
      
      if (verbose) {
        message('Training a FlowSOM normalization model')
      }
      
      fsom.nCells    <- fsom_params[['nCells']]
      fsom.colsToUse <- fsom_params[['colsToUse']]
      if (is.null(fsom.colsToUse) && !is.null(cols)) {
        fsom.colsToUse <- cols
      }
      pars           <- fsom_params[
        grep('nCells|channels|colsToUse', names(fsom_params), invert = TRUE)
      ]
      suppressMessages({
        fsom <- CytoNorm::prepareFlowSOM(
          files          = fnames,
          colsToUse      = fsom.colsToUse,
          nCells         = fsom.nCells,
          FlowSOM.params = pars,
          seed           = seed,
          verbose        = FALSE,
          ...
        )
      })
    } else {
      
      if (verbose) {
        message('Using a precomputed FlowSOM normalization model')
      }
    }
    
    ## Additional consistency check
    
    if ('goal'%in%names(norm_params) && is.list(norm_params[['goal']])) {
      stopifnot(
        'Metacluster count in goal distribution and FlowSOM model is different' =
          length(norm_params[['goal']])==FlowSOM::NMetaclusters(fsom)
      )
    }
    
    ## Split samples by metaclusters
    
    if (verbose) {
      message(
        'Splitting training samples into ', FlowSOM::NMetaclusters(fsom) ,
        ' FlowSOM metaclusters'
      )
    }
    
    meta_split <-
      .SplitSamplesByMetaclusters(
        fnames    = fnames,
        batches   = batches,
        fsom      = fsom,
        fpath_out = tmp,
        cols      = intersect(colnames(fsom$map$codes), cols),
        cores     = cores,
        verbose   = verbose,
        ...
      )
    
    ## Import multi-threading functions
    
    makeCluster    <- snow::makeCluster
    registerDoSNOW <- doSNOW::registerDoSNOW
    foreach        <- foreach::foreach
    `%dopar%`      <- foreach::`%dopar%`
    stopCluster    <- snow::stopCluster
    
    ## Model metacluster-channel distributions by batch
    
    if (verbose) {
      message('Modeling metacluster-channel distributions by batch')
    }
    
    if (is.null(cols)) {
      cols <- GetPanels(fnames[1])[[1]]$Channel
    }
    
    umcl <- levels(FlowSOM::GetMetaclusters(fsom))
    nmcl <- length(umcl)
    
    packages <- c()
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = nmcl, style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    mcl_res <- foreach( # iterate over metaclusters asynchronously
      i             = seq_along(umcl),
      .combine      = c,
      .inorder      = TRUE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      mcl <- umcl[i]
      
      mask   <- meta_split$Metacluster==mcl
      files  <- meta_split$FileName[mask]
      labels <- as.character(meta_split$Batch[mask])
      
      pars <- c(
        list(
          'files'         = files,
          'labels'        = labels,
          'channels'      = cols,
          'transformList' = NULL,
          'verbose'       = verbose,
          'plot'          = FALSE
        ),
        norm_params,
        list(...)
      )
      if (is.list(pars[['goal']])) {
        pars[['goal']] <- pars[['goal']][[mcl]]
      }
      
      suppressMessages({
        res <- list(do.call(
          normMethod.train, pars
        ))
      })
      names(res) <- mcl
      return(res)
    }
    
    if (verbose) {
      close(pb)
    }
    stopCluster(clu)
    invisible(gc())
    rm(clu)
    invisible(gc())
    
    invisible(suppressWarnings({
      file.remove(meta_split$FileName)  
    }))
    
    res <- list(
      'fsom'       = fsom,
      'clusterRes' = mcl_res
    )
  } else {
    
    ## Split samples by gates
    
    ng <- length(gate_names)
    
    if (verbose) {
      message(
        'Splitting training samples into ', ng, ' gates:\n\t',
        paste0(gate_names, collapse = ', ')
      )
    }
    
    gate_split <-
      .SplitSamplesByGates(
        fnames     = fnames,
        batches    = batches,
        gating_fun = gating_fun,
        gate_names = gate_names,
        fpath_out  = tmp,
        cols       = cols,
        cores      = cores,
        verbose    = verbose,
        ...
      )
    
    ## Import multi-threading functions
    
    makeCluster    <- snow::makeCluster
    registerDoSNOW <- doSNOW::registerDoSNOW
    foreach        <- foreach::foreach
    `%dopar%`      <- foreach::`%dopar%`
    stopCluster    <- snow::stopCluster
    
    ## Model gate-channel distributions by batch
    
    if (verbose) {
      message('Modeling gate-channel distributions by batch')
    }
    
    if (is.null(cols)) {
      cols <- GetPanels(fnames[1])[[1]]$Channel
    }
    
    packages <- c()
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = ng, style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    gate_res <- foreach( # iterate over gates asynchronously
      i             = seq_along(gate_names),
      .combine      = c,
      .inorder      = TRUE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      gate <- gate_names[i]
      
      mask   <- gate_split$Gate==gate
      files  <- gate_split$FileName[mask]
      labels <- as.character(gate_split$Batch[mask])
      
      pars <- c(
        list(
          'files'         = files,
          'labels'        = labels,
          'channels'      = cols,
          'transformList' = NULL,
          'verbose'       = verbose,
          'plot'          = FALSE
        ),
        norm_params,
        list(...)
      )
      if (is.list(pars[['goal']])) {
        pars[['goal']] <- pars[['goal']][[gate]]
      }
      
      suppressMessages({
        res <- list(do.call(
          normMethod.train, pars
        ))
      })
      names(res) <- gate
      return(res)
    }
    
    if (verbose) {
      close(pb)
    }
    stopCluster(clu)
    invisible(gc())
    rm(clu)
    invisible(gc())
    
    invisible(suppressWarnings({
      file.remove(gate_split$FileName)  
    }))
    
    attributes(gating_fun)$gate_names <- gate_names
    res <- list(
      'gating'     = gating_fun,
      'clusterRes' = gate_res
    )
  }
  
  
  attributes(res)$class           <- 'cytoSNOW_CytoNormModel'
  attributes(res)$splitting_model <- ifelse(use_gating, 'gating', 'fsom')
  attributes(res)$n_cols          <- n_cols
  attributes(res)$cols            <- cols
  attributes(res)$n_batch         <- n_batch
  attributes(res)$n_fcs           <- n_fcs
  
  return(res)
}

#' Parallel application of normalization model
#'
#' Applies a CytoNorm batch effect correction model, trained via 
#' [cytoSNOW::ParallelNormalize.Train()], to FCS files affected by batch effect.
#'
#' @param model normalization model generated with 
#' [cytoSNOW::ParallelNormalize.train()]
#' @param fnames vector. Full paths to FCS files which should be normalized
#' @param batches vector. Batch labels per file in `fnames`
#' @param fpath_out string. Path to directory to save normalized files
#' @param keep_full_paths logical. Whether to include the full path to each FCS 
#' file in the name of its normalized counterpart. Defaults to FALSE
#' @param prefix string. Prefix to glue to the beginning of the original file 
#' names in `fnames` to denote normalized. Defaults to *'Norm_'*
#' @param normMethod.normalize application function of the signal normalization 
#' method to use. Defaults to `CytoNorm::QuantileNorm.normalize`
#' @param cores integer. Number of CPU cores to use for parallelization (at
#' least 2). Defaults to number of detectable cores minus 1
#' @param verbose logical. Whether to indicate progress. Defaults to `TRUE`
#' @param ... optional additional named parameters for [flowCore::read.FCS()]
#'
#' @details
#' This function uses parallelization via a SNOW cluster for speed-up.
#'
#' @note
#' It is highly advisable to have unique file names for each input FCS sample 
#' and have all of them in the same directory (no nested directories). An 
#' alternative solution is to take the entire paths to each sample as unique 
#' identifiers (set `keep_full_paths` to TRUE), but this is discouraged.
#'
#' @return nothing is returned
#'
#' @seealso [CytoNorm::CytoNorm.normalize()] for the original CytoNorm 
#' normalization function
#' @seealso [cytoSNOW::ParallelNormalize.Train()] for training the normalization 
#' model that this function can then apply
#' 
#' @export
ParallelNormalize.Apply <- ParallelNormalise.Apply <- function(
  model,
  fnames,
  batches,
  fpath_out,
  keep_full_paths      = FALSE,
  prefix               = 'Norm_',
  normMethod.normalize = CytoNorm::QuantileNorm.normalize,
  cores                = parallel::detectCores()-1,
  verbose              = TRUE,
  ...
) {
  
  ## Validate inputs
  
  stopifnot(
    '`model` must be generated using `cytoSNOW::ParallelNormalize.Train`' =
      class(model)=='cytoSNOW_CytoNormModel'
  )
  stopifnot('`fnames` must be a non-empty non-list string vector' =
              is.vector(fnames) && is.character(fnames) && !is.list(fnames) &&
              length(fnames)>0)
  stopifnot('Some files in `fnames` are missing' = all(file.exists(fnames)))
  
  stopifnot('`batches` must be a non-empty non-list string or numeric vector' =
              is.vector(batches) &&
              (is.character(batches) || is.numeric(batches)) &&
              !is.list(batches) && length(batches)>0)
  stopifnot(
    '`batches` must be of the same length as `fnames' =
      length(batches)==length(fnames)
  )
  if (Sys.getenv('DUPLICATE_EXCEPTION')!='TRUE') {
    
    stopifnot(
      '`fnames` must not contain duplicates' =
        length(unique(fnames))==length(fnames)
    )
  }
  stopifnot('`fpath_out` must be a single string' =
              is.character(fpath_out) && length(fpath_out)==1)
  stopifnot('`keep_full_paths` must be a single logical' = 
              is.logical(keep_full_paths) && length(keep_full_paths)==1)
  stopifnot('`prefix` must be a single string' =
              is.character(prefix) && length(prefix)==1)
  stopifnot('`normMethod.normalize` must be a function' =
              is.function(normMethod.normalize))
  stopifnot('`verbose` must be a single logical' =
              is.logical(verbose) && length(verbose)==1)
  
  ## Determine a location for temporary files
  
  tmp <- tempdir(check = TRUE)
  
  ## Create results directory
  
  if (!dir.exists(fpath_out)) {
    res_fpath_out <- tryCatch(
      expr    = { dir.create(fpath_out) },
      error   = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!res_fpath_out) {
      stop('Cannot create the results directory "', fpath_out, '"')
    }
  }
  
  ## Init
  
  n_fcs   <- length(fnames)
  n_batch <- length(unique(batches))
  cols    <- attributes(model)$cols
  n_cols  <- attributes(model)$n_cols
  
  use_gating <- attributes(model)$splitting_model=='gating'
  
  if (verbose) {
    message(paste0(
      'Applying CytoNorm normalization model for batch-effect correction',
      '\n• Number of normalizable channels: ', n_cols,
      '\n• Input data: ', n_fcs, ' FCS ', if (n_fcs>1) 'files' else 'file',
      ' across ', n_batch, ' experimental ',
      if (n_batch>1) 'batches' else 'batch',
      '\n• Separate normalization functions will be applied per ',
      if (!use_gating) { 'FlowSOM metacluster' } else { 'gate' },
      '\n• Using ', cores, ' CPU cores'
    ))
  }
  
  if (!use_gating) {
    
    ## Split samples by metaclusters
    
    fsom <- model$fsom
    nmcl <- FlowSOM::NMetaclusters(fsom)
    
    if (verbose) {
      message('Splitting target samples into ', nmcl, ' FlowSOM metaclusters')
    }
    
    meta_split <-
      .SplitSamplesByMetaclusters(
        fnames    = fnames,
        batches   = batches,
        fsom      = fsom,
        fpath_out = tmp,
        cols      = intersect(colnames(fsom$map$codes), cols),
        cores     = cores,
        verbose   = verbose,
        ...
      )
    
    ## Import multi-threading functions
    
    makeCluster    <- snow::makeCluster
    registerDoSNOW <- doSNOW::registerDoSNOW
    foreach        <- foreach::foreach
    `%dopar%`      <- foreach::`%dopar%`
    stopCluster    <- snow::stopCluster
    
    ## Apply metacluster-channel normalizations
    
    if (verbose) {
      message('Applying metacluster-channel normalizations by batch')
    }
    
    if (is.null(cols)) {
      cols <- GetPanels(fnames[1])[[1]]$Channel
    }
    
    umcl <- levels(FlowSOM::GetMetaclusters(fsom))
    nmcl <- length(umcl)
    
    packages <- c()
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = nmcl, style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    meta_norm <- foreach( # iterate over metaclusters asynchronously
      i             = seq_along(umcl),
      .combine      = rbind,
      .inorder      = TRUE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      mcl <- umcl[i]
      mask <-
        meta_split$Metacluster==mcl & file.exists(meta_split$FileName)
      # ^ allow for non-existent split files if metacluster empty
      
      if (!any(mask)) {
        return(NULL)
      }
      
      fnames_mcl  <- meta_split$FileName[mask]
      batches_mcl <- meta_split$Batch[mask]
      samples_mcl <- meta_split$Sample[mask]
      
      normMethod.normalize(
        model                 = model$clusterRes[[mcl]],
        files                 = fnames_mcl,
        labels                = batches_mcl,
        outputDir             = fpath_out,
        prefix                = prefix,
        transformList         = NULL,
        transformList.reverse = NULL,
        removeOriginal        = FALSE,
        verbose               = FALSE,
        ...
      )
      
      return(
        data.frame(
          'Metacluster' = mcl,
          'FileName'    =
            file.path(fpath_out, paste0(prefix, basename(fnames_mcl))),
          'Sample'      = samples_mcl,
          'Batch'       = batches_mcl
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
    
    invisible(suppressWarnings({
      file.remove(meta_split$FileName)  
    }))
    
    ## Stitch normalized split data back together
    
    if (verbose) {
      message('Stitching normalized split data back together')
    }
    
    packages <- c('flowCore', 'FlowSOM')
    
    codes <- fsom$map$codes[
      , intersect(cols, colnames(fsom$map$codes)), drop = FALSE
    ]
    meta  <- fsom$metaclustering
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = length(fnames), style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    out <- foreach( # iterate over files asynchronously
      i             = seq_along(fnames),
      .inorder      = FALSE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      fname <- fnames[i]
      suppressWarnings({
        ff <- read.FCS(fname, ...)
      })
      d <- ff@exprs[, cols, drop = FALSE]
      
      mcls <- # metacluster label per cell
        meta[FlowSOM:::MapDataToCodes(
          codes   = codes[, intersect(colnames(codes), colnames(d))],
          newdata = d
        )[, 1]]
      
      for (j in seq_along(umcl)) {
        
        mcl       <- umcl[j]
        mask      <- meta_norm$Sample==fname & meta_norm$Metacluster==mcl
        if (any(mask)) {
          
          fname_mcl <- meta_norm$FileName[mask]
          ff_subset <- flowCore::read.FCS(fname_mcl, ...)
          ff@exprs[mcls==mcl, ] <- ff_subset@exprs
        }
      }
      
      ff <- ValidateFCS(exprs = ff@exprs, ff = ff)
      fname_out <-
        file.path(
          fpath_out,
          if (keep_full_paths) {
            paste0(prefix, gsub('[:/]', '_', fname))
          } else {
            paste0(prefix, basename(fname))
          }
        )
      suppressWarnings({
        flowCore::write.FCS(ff, fname_out)
      })
      return(invisible(NULL))
    }
    
    if (verbose) {
      close(pb)
    }
    stopCluster(clu)
    invisible(gc())
    rm(clu)
    invisible(gc())
    
    invisible(suppressWarnings({
      file.remove(meta_norm$FileName)  
    }))
    
  } else {
    
    ## Split samples by gates
    
    gating_fun <- model$gating
    gate_names <- attributes(gating_fun)$gate_names
    ng         <- length(gate_names)
    
    if (verbose) {
      message(
        'Splitting target samples into ', ng, ' gates:\n\t',
        paste0(gate_names, collapse = ', ')
      )
    }
    
    gate_split <-
      .SplitSamplesByGates(
        fnames     = fnames,
        batches    = batches,
        gating_fun = gating_fun,
        gate_names = gate_names,
        fpath_out  = tmp,
        cols       = cols,
        cores      = cores,
        verbose    = verbose,
        ...
      )
    
    ## Import multi-threading functions
    
    makeCluster    <- snow::makeCluster
    registerDoSNOW <- doSNOW::registerDoSNOW
    foreach        <- foreach::foreach
    `%dopar%`      <- foreach::`%dopar%`
    stopCluster    <- snow::stopCluster
    
    ## Apply gate-channel normalizations
    
    if (verbose) {
      message('Applying gate-channel normalizations by batch')
    }
    
    if (is.null(cols)) {
      cols <- GetPanels(fnames[1])[[1]]$Channel
    }
    
    packages <- c()
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = ng, style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    gate_norm <- foreach( # iterate over gates asynchronously
      i             = seq_along(gate_names),
      .combine      = rbind,
      .inorder      = TRUE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      gate <- gate_names[i]
      mask <-
        gate_split$Gate==gate & file.exists(gate_split$FileName)
      # ^ allow for non-existent split files if metacluster empty
      
      if (!any(mask)) {
        return(NULL)
      }
      
      fnames_gate  <- gate_split$FileName[mask]
      batches_gate <- gate_split$Batch[mask]
      samples_gate <- gate_split$Sample[mask]
      
      normMethod.normalize(
        model                 = model$clusterRes[[gate]],
        files                 = fnames_gate,
        labels                = batches_gate,
        outputDir             = fpath_out,
        prefix                = prefix,
        transformList         = NULL,
        transformList.reverse = NULL,
        removeOriginal        = FALSE,
        verbose               = FALSE,
        ...
      )
      
      return(
        data.frame(
          'Gate'        = gate,
          'FileName'    =
            file.path(fpath_out, paste0(prefix, basename(fnames_gate))),
          'Sample'      = samples_gate,
          'Batch'       = batches_gate
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
    
    invisible(suppressWarnings({
      file.remove(gate_split$FileName)  
    }))
    
    ## Stitch normalized split data back together
    
    if (verbose) {
      message('Stitching normalized split data back together')
    }
    
    packages <- c('flowCore')
    
    clu <- makeCluster(cores)
    registerDoSNOW(clu)
    
    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = length(fnames), style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }
    
    out <- foreach( # iterate over files asynchronously
      i             = seq_along(fnames),
      .inorder      = FALSE,
      .final        = NULL,
      .packages     = packages,
      .options.snow = opts
    ) %dopar% {
      
      fname <- fnames[i]
      suppressWarnings({
        ff <- read.FCS(fname, ...)
      })
      d <- ff@exprs[, cols, drop = FALSE]
      
      gates <- gating_fun(ff)
      
      for (j in seq_along(gate_names)) {
        
        gate <- gate_names[j]
        mask <- gate_norm$Sample==fname & gate_norm$Gate==gate
        if (any(mask)) {
          
          fname_gate <- gate_norm$FileName[mask]
          ff_subset <- flowCore::read.FCS(fname_gate, ...)
          ff@exprs[gates==gate, ] <- ff_subset@exprs
        }
      }
      
      ff <- ValidateFCS(exprs = ff@exprs, ff = ff)
      fname_out <-
        file.path(
          fpath_out,
          if (keep_full_paths) {
            paste0(prefix, gsub('[:/]', '_', fname))
          } else {
            paste0(prefix, basename(fname))
          }
        )
      suppressWarnings({
        flowCore::write.FCS(ff, fname_out)
      })
      return(invisible(NULL))
    }
    
    if (verbose) {
      close(pb)
    }
    stopCluster(clu)
    invisible(gc())
    rm(clu)
    invisible(gc())
    
    invisible(suppressWarnings({
      file.remove(gate_norm$FileName)  
    }))
  }
  
  return(invisible(NULL))
}
