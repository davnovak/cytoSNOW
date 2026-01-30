
#' Make a valid flowFrame
#'
#' Generates a `flowCore::flowFrame` based on an expression matrix and optional
#' template flowFrame and/or `desc` element.
#' Required slots of the flowFrame XML are filled in, which should prevent error
#' messages that can otherwise occur when the `exprs` slot of an existing
#' flowFrame gets rewritten.
#'
#' @param exprs numeric matrix with named columns. Row-wise cytometry
#' expression matrix
#' @param ff `flowCore::flowFrame` or `NULL`. Template flowFrame, on which the
#' new flowFrame should be based. Defaults to `NULL`
#' @param desc named string vector or `NULL`. Template `desc` vector (*e.g.*,
#' `ff@parameters@data$desc`) to use. Defaults to `NULL`
#'
#' @return `flowCore::flowFrame`
#'
#' @export
ValidateFCS <- function(
    exprs,
    ff     = NULL,
    desc   = NULL
) {
  
  if(!is.null(ff)) {
    
    mn <- flowCore::markernames(ff)
    ff <- flowCore::flowFrame(exprs, parameters = flowCore::parameters(ff))
    flowCore::markernames(ff) <- mn
  } else {
    
    mn <- `names<-`(colnames(exprs), colnames(exprs))
    ff <- flowCore::flowFrame(exprs)
  }
  
  params <- flowCore::parameters(ff)
  pd     <- NULL
  cols   <- as.vector(pd$name)
  idcs   <- match(cols, pd$name)
  
  if (any(is.na(idcs))) {
    stop('Invalid column specifier')
  }
  
  keyval <- list()
  for (channel_number in seq_len(ncol(exprs))) {
    
    channel_name <- colnames(exprs)[channel_number]
    if (is.null(desc)) {
      desc <- colnames(exprs)[channel_number]
    }
    channel_id    <- paste('$P', channel_number, sep = '')
    channel_range <- max(exprs[, channel_number]) + 1
    channel_min   <- min(0, min(exprs[, channel_number])-1)
    
    plist <- matrix(
      c(
        channel_name,
        desc[channel_number],
        channel_range,
        channel_min,
        channel_range-1
      )
    )
    rownames(plist) <- c('name', 'desc', 'range', 'minRange', 'maxRange')
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    
    keyval[[paste0('$P', channel_number, 'B')]] <- '32'
    keyval[[paste0('$P', channel_number, 'R')]] <- toString(channel_range)
    keyval[[paste0('$P', channel_number, 'E')]] <- '0,0'
    keyval[[paste0('$P', channel_number, 'N')]] <- channel_name
    keyval[[paste0('$P', channel_number, 'S')]] <- channel_name
  }
  
  params@data <- data.frame(pd)
  if(!is.null(ff)) {
    
    ff <- flowCore::flowFrame(exprs, parameters = params)
  } else {
    
    ff <- flowCore::flowFrame(exprs, parameters = params)
  }
  flowCore::keyword(ff) <- keyval
  
  if (!is.null(mn)) {
    flowCore::markernames(ff) <- mn
  }
  
  ff
}
