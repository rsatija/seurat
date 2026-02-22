#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(testthat)
  library(jsonlite)
})

source(file.path("tests", "testthat", "helper_core_golden_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else core_golden_root()

message("Generating core goldens in: ", output_dir)
invisible(generate_core_goldens(output_dir = output_dir))
message("Core golden generation complete")
