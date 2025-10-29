# SpiecEasi NEWS

## Unreleased changes

### BUG FIXES

* Fixed BugReports URL in DESCRIPTION (removed double slash)
* Replaced sapply with vapply in R/fitdistr.R for improved type stability
* Fixed incomplete final line in .Rbuildignore file

### IMPROVEMENTS

* Added BiocManager installation instructions to main vignette
* Added runnable examples to man pages for getOptInd, robustPCA, stars.roc, and triu functions
* Cleaned up promotional language in NEWS entries
* Removed unused inst/extdata files to reduce package size

### INFRASTRUCTURE

* Added .registration = TRUE to useDynLib in NAMESPACE
* Updated package version to 0.99.2 for Bioconductor submission
* Verified vignette runtime compliance with Bioconductor guidelines

## Changes in version 1.99.2

### NEW FEATURES

* Added comprehensive Bioconductor package infrastructure
* Added automated GitHub Actions workflow for Bioconductor continuous integration
* Added test coverage reporting and automated testing
* Added Bioconductor-style vignettes and documentation
    - Added biocViews categorization: Software, Microbiome, Metagenomics, GraphAndNetwork, NetworkInference
    - Added comprehensive support documentation and issue templates
    - Implemented proper code of conduct and contribution guidelines

### INFRASTRUCTURE IMPROVEMENTS

* Added Bioconductor-specific GitHub Actions workflow (check-bioc.yml)
* Implemented automated R CMD check with Bioconductor standards
* Added BiocCheck integration for package validation
* Enhanced test coverage with automated reporting
* Added lifecycle badge indicating stable package status

### BUG FIXES

* Resolved compliance issues identified during BiocCheck
* Enhanced error handling and package robustness

## Changes in version 1.1.3

### NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

### SIGNIFICANT USER-VISIBLE CHANGES

* Adding other boilerplate for Bioconductor release

### BUG FIXES

* Your bug fixes. See more details at
<http://bioconductor.org/developers/package-guidelines/#news>.
