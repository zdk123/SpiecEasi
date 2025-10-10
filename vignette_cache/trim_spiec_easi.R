# Function to trim spiec.easi objects for caching
# Removes large intermediate data while keeping essential plotting/display info

trim_spiec_easi <- function(se) {
  if (is.null(se) || !inherits(se, "pulsar.refit")) {
    return(se)
  }
  
  # Keep only essential components
  trimmed <- list(
    refit = se$refit,           # Final selected network
    lambda = se$lambda,         # Lambda values
    fun = se$fun               # Function used
  )
  
  # Keep minimal select info for plotting and getOptInd()
  if (!is.null(se$select) && !is.null(se$select$stars)) {
    trimmed$select <- list(
      stars = list(
        summary = se$select$stars$summary,  # For plotting
        opt.index = 1,  # Always set to 1 since we only keep the optimal element
        merge = se$select$stars$merge  # For getOptMerge()
      )
    )
  }
  
  # Keep minimal est info - but preserve what's needed for getOptCov/getOptBeta
  if (!is.null(se$est)) {
    # Get the optimal index safely
    opt_idx <- if (!is.null(se$select$stars$opt.index)) se$select$stars$opt.index else 1
    
    trimmed$est <- list(
      method = se$est$method,
      lambda = se$est$lambda,
      sparsity = se$est$sparsity,
      df = se$est$df,
      sym = se$est$sym,
      # Keep the optimal path elements as the first (and only) element in lists
      path = if (!is.null(se$est$path) && length(se$est$path) >= opt_idx) list(se$est$path[[opt_idx]]) else NULL,
      beta = if (!is.null(se$est$beta) && length(se$est$beta) >= opt_idx) list(se$est$beta[[opt_idx]]) else NULL,
      cov = if (!is.null(se$est$cov) && length(se$est$cov) >= opt_idx) list(se$est$cov[[opt_idx]]) else NULL,
      # Keep data and resid for latent variable models
      data = se$est$data,
      resid = if (!is.null(se$est$resid) && length(se$est$resid) >= opt_idx) list(se$est$resid[[opt_idx]]) else NULL
    )
  }
  
  class(trimmed) <- class(se)
  return(trimmed)
}

# Function to trim multiple spiec.easi objects
trim_spiec_easi_list <- function(se_list) {
  if (is.list(se_list)) {
    lapply(se_list, trim_spiec_easi)
  } else {
    trim_spiec_easi(se_list)
  }
}
