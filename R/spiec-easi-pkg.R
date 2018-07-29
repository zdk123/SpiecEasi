#' @importFrom graphics abline legend par plot
#' @importFrom methods as
#' @importFrom stats cor cov cov2cor lm median  na.exclude na.omit optim pnorm ppoints qlnorm rmultinom rnorm rpois sd var
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
NULL

#' @rdname hmp2
"hmp216S"

#' @rdname hmp2
"hmp2prot"
