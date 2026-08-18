# Backward-compatible function loader for scripts that sourced the original file.
# The implementation is now split into focused, testable modules.
suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
})

source("R/io.R")
source("R/regions.R")
source("R/analysis.R")
source("R/primers.R")
