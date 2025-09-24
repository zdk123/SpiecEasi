#' @title SpiecEasi: Sparse Inverse Covariance Estimation for Ecological Association Inference
#' @description SpiecEasi provides methods for inferring ecological associations from compositional microbiome data using sparse inverse covariance estimation.
#' @details The package implements SPIEC-EASI (SParse InversE Covariance estimation for Ecological Association Inference) methods for microbiome network inference, including glasso and neighborhood selection approaches with stability selection for model selection.
#' @return This package provides functions for microbiome network inference and analysis
#' @importFrom graphics abline legend par plot
#' @importFrom methods as
#' @importFrom stats cor cov cov2cor lm median na.exclude na.omit optim pnorm ppoints qlnorm rmultinom rnorm rpois sd var
#' @keywords internal
"_PACKAGE"


#' @name AGP
#' @title American Gut Project
#' @description Round 1 and 2 community count datasets from the American Gut Project.
#' @format \enumerate{
#'  \item amgut1.filt: A matrix with 289 samples (rows) and 127 OTUs (cols).
#' \item amgut2.filt.phy: A phyloseq object
#' }
#' @source http://humanfoodproject.com/americangut/
#' @return List containing amgut1.filt matrix and amgut2.filt.phy phyloseq object
NULL


#' @rdname AGP
#' @usage data(amgut1.filt)
#' @name amgut1.filt
NULL

#' @rdname AGP
#' @usage data(amgut2.filt.phy)
#' @name amgut2.filt.phy
NULL

#' Human Microbiome Project 2
#'
#' Pre-filtered data from the intregrated human microbiome project.
#' @name hmp2
#' @docType data
#' @usage data(hmp2)
#' @format \enumerate{
#'    \item hmp216S: 16S data, phyloseq object 45 taxa and 47 samples.
#'    \item hmp2prot: protein data, A phyloseq object, 43 'taxa' and 47 samples.
#'  }
#' @source https://www.hmpdacc.org/ihmp/
#' @return List containing hmp216S and hmp2prot phyloseq objects
NULL

#' @name hmp216S
#' @rdname hmp2
NULL

#' @name hmp2prot
#' @rdname hmp2
NULL
