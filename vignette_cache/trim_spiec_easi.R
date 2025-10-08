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
        opt.index = se$select$stars$opt.index,  # For getOptInd()
        merge = se$select$stars$merge  # For getOptMerge()
      )
    )
  }
  
  # Keep minimal est info - but preserve what's needed for getOptCov/getOptBeta
  if (!is.null(se$est)) {
    trimmed$est <- list(
      method = se$est$method,
      lambda = se$est$lambda,
      sparsity = se$est$sparsity,
      df = se$est$df,
      sym = se$est$sym,
      # Keep the optimal path elements for getOptCov/getOptBeta
      path = se$est$path[se$select$stars$opt.index],  # Only optimal path
      beta = se$est$beta[se$select$stars$opt.index],  # Only optimal beta
      cov = list(se$est$cov[[se$select$stars$opt.index]])  # Only optimal cov as list
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
