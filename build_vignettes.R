#!/usr/bin/env Rscript
# Build vignettes with optional force rebuild
# Usage: 
#   Rscript build_vignettes.R                    # Use cached data (fast)
#   Rscript build_vignettes.R --force            # Force rebuild (slow)
#   FORCE_VIGNETTE_REBUILD=TRUE Rscript build_vignettes.R  # Force rebuild via env var

args <- commandArgs(trailingOnly = TRUE)
force_rebuild <- "--force" %in% args || Sys.getenv("FORCE_VIGNETTE_REBUILD", "FALSE") == "TRUE"

if (force_rebuild) {
  message("🔄 FORCE REBUILD: Rebuilding all vignettes from scratch...")
  Sys.setenv(FORCE_VIGNETTE_REBUILD = "TRUE")
} else {
  message("⚡ FAST BUILD: Using cached vignette data...")
  Sys.setenv(FORCE_VIGNETTE_REBUILD = "FALSE")
}

message(paste("FORCE_VIGNETTE_REBUILD:", Sys.getenv("FORCE_VIGNETTE_REBUILD")))

# Build vignettes
message("Building vignettes...")
tools::buildVignettes(dir = ".", tangle = TRUE)
message("✅ Vignette build complete!")
