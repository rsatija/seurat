#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

read_benchmark_file <- function(path) {
  if (!file.exists(path)) {
    stop("Benchmark file not found: ", path, call. = FALSE)
  }
  fromJSON(path, simplifyVector = FALSE)
}

as_metric_map <- function(metrics) {
  stats::setNames(metrics, vapply(metrics, function(x) x$scenario_id, character(1)))
}

compare_benchmark_to_baseline <- function(
  current_path = file.path("tests", "benchmarks", "current", "core_v1_metrics.json"),
  baseline_path = file.path("tests", "benchmarks", "baseline", "core_v1_metrics.json"),
  threshold = 0.05
) {
  current <- read_benchmark_file(current_path)
  baseline <- read_benchmark_file(baseline_path)

  current_map <- as_metric_map(current$metrics)
  baseline_map <- as_metric_map(baseline$metrics)

  missing <- setdiff(names(baseline_map), names(current_map))
  if (length(missing) > 0) {
    stop("Current benchmark is missing scenarios: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  results <- vector("list", length(baseline_map))
  names(results) <- names(baseline_map)
  regressions <- character()

  i <- 1L
  for (scenario_id in names(baseline_map)) {
    b <- baseline_map[[scenario_id]]
    c <- current_map[[scenario_id]]

    runtime_ratio <- (c$median_elapsed_sec / b$median_elapsed_sec) - 1
    alloc_ratio <- (c$median_allocated_bytes / b$median_allocated_bytes) - 1

    status <- "ok"
    if (is.finite(runtime_ratio) && runtime_ratio > threshold) {
      status <- "regression"
      regressions <- c(
        regressions,
        sprintf("%s runtime regression %.4f > %.4f", scenario_id, runtime_ratio, threshold)
      )
    }
    if (is.finite(alloc_ratio) && alloc_ratio > threshold) {
      status <- "regression"
      regressions <- c(
        regressions,
        sprintf("%s allocation regression %.4f > %.4f", scenario_id, alloc_ratio, threshold)
      )
    }

    results[[i]] <- list(
      scenario_id = scenario_id,
      baseline_runtime_sec = b$median_elapsed_sec,
      current_runtime_sec = c$median_elapsed_sec,
      runtime_ratio = runtime_ratio,
      baseline_allocated_bytes = b$median_allocated_bytes,
      current_allocated_bytes = c$median_allocated_bytes,
      allocated_ratio = alloc_ratio,
      threshold = threshold,
      status = status
    )
    i <- i + 1L
  }

  list(
    ok = length(regressions) == 0,
    regressions = regressions,
    results = unname(results)
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  current_path <- if (length(args) >= 1) args[[1]] else file.path("tests", "benchmarks", "current", "core_v1_metrics.json")
  baseline_path <- if (length(args) >= 2) args[[2]] else file.path("tests", "benchmarks", "baseline", "core_v1_metrics.json")
  threshold <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.05

  comparison <- compare_benchmark_to_baseline(
    current_path = current_path,
    baseline_path = baseline_path,
    threshold = threshold
  )

  if (!comparison$ok) {
    stop("Benchmark comparison failed:\n", paste(comparison$regressions, collapse = "\n"), call. = FALSE)
  }

  message("Benchmark comparison passed for all scenarios")
}
