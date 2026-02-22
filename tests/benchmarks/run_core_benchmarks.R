#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

source(file.path("tests", "testthat", "helper_core_golden_utils.R"))
source(file.path("tests", "benchmarks", "compare_core_benchmarks.R"))

clone_object <- function(x) {
  unserialize(serialize(x, NULL))
}

parse_rprofmem <- function(path) {
  if (!file.exists(path)) {
    return(0)
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    return(0)
  }
  bytes <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", lines)))
  bytes <- bytes[is.finite(bytes)]
  if (length(bytes) == 0) {
    return(0)
  }
  sum(bytes)
}

measure_once <- function(f) {
  memfile <- tempfile("core-bench-mem")
  on.exit(unlink(memfile), add = TRUE)

  gc()
  Rprofmem(memfile)
  elapsed <- system.time(f())["elapsed"]
  Rprofmem(NULL)

  list(
    elapsed_sec = as.numeric(elapsed),
    allocated_bytes = as.numeric(parse_rprofmem(memfile))
  )
}

run_case <- function(scenario_id, fun, warmup, runs) {
  for (i in seq_len(warmup)) {
    invisible(measure_once(fun))
  }

  elapsed <- numeric(runs)
  allocated <- numeric(runs)
  for (i in seq_len(runs)) {
    m <- measure_once(fun)
    elapsed[[i]] <- m$elapsed_sec
    allocated[[i]] <- m$allocated_bytes
  }

  list(
    scenario_id = scenario_id,
    warmup_runs = warmup,
    measured_runs = runs,
    elapsed_sec = unname(elapsed),
    allocated_bytes = unname(allocated),
    median_elapsed_sec = as.numeric(stats::median(elapsed)),
    median_allocated_bytes = as.numeric(stats::median(allocated))
  )
}

build_benchmark_cases <- function() {
  pbmc <- core_fixture_pbmc_raw()

  base_norm <- CreateSeuratObject(counts = pbmc)
  base_norm <- NormalizeData(base_norm, verbose = FALSE)
  base_norm <- FindVariableFeatures(base_norm, verbose = FALSE)

  base_scaled <- ScaleData(clone_object(base_norm), verbose = FALSE)
  base_pca <- RunPCA(clone_object(base_scaled), npcs = 10, verbose = FALSE)
  base_neighbors <- FindNeighbors(clone_object(base_pca), dims = 1:10, k.param = 10, verbose = FALSE)

  transfer_fixture <- create_transfer_fixture()
  transfer_anchors <- FindTransferAnchors(
    reference = transfer_fixture$ref,
    query = transfer_fixture$query,
    k.filter = 50
  )

  integrate_fixture <- create_integrate_fixture()
  integrate_anchors <- suppressMessages(
    suppressWarnings(
      FindIntegrationAnchors(
        object.list = c(integrate_fixture$ref, integrate_fixture$query.list[[1]]),
        k.filter = NA,
        verbose = FALSE
      )
    )
  )

  list(
    list(
      id = "bench_scale_data",
      fn = function() {
        ScaleData(clone_object(base_norm), verbose = FALSE)
      }
    ),
    list(
      id = "bench_run_pca",
      fn = function() {
        RunPCA(clone_object(base_scaled), npcs = 10, verbose = FALSE)
      }
    ),
    list(
      id = "bench_run_umap",
      fn = function() {
        suppressWarnings(RunUMAP(clone_object(base_pca), dims = 1:10, seed.use = 42, verbose = FALSE))
      }
    ),
    list(
      id = "bench_find_neighbors",
      fn = function() {
        FindNeighbors(clone_object(base_pca), dims = 1:10, k.param = 10, verbose = FALSE)
      }
    ),
    list(
      id = "bench_find_clusters",
      fn = function() {
        FindClusters(clone_object(base_neighbors), random.seed = 42, verbose = FALSE)
      }
    ),
    list(
      id = "bench_find_transfer_anchors",
      fn = function() {
        FindTransferAnchors(
          reference = clone_object(transfer_fixture$ref),
          query = clone_object(transfer_fixture$query),
          k.filter = 50
        )
      }
    ),
    list(
      id = "bench_integrate_data",
      fn = function() {
        suppressWarnings(IntegrateData(anchorset = clone_object(integrate_anchors), k.weight = 50, verbose = FALSE))
      }
    ),
    list(
      id = "bench_transfer_data",
      fn = function() {
        TransferData(
          anchorset = clone_object(transfer_anchors),
          refdata = transfer_fixture$ref$RNA_snn_res.1,
          verbose = FALSE
        )
      }
    )
  )
}

run_core_benchmarks <- function(
  enforce = FALSE,
  output_path = file.path("tests", "benchmarks", "current", "core_v1_metrics.json"),
  baseline_path = file.path("tests", "benchmarks", "baseline", "core_v1_metrics.json"),
  threshold = 0.05,
  warmup = 1,
  runs = 7
) {
  cases <- build_benchmark_cases()
  metrics <- lapply(cases, function(case) {
    message("Running benchmark: ", case$id)
    run_case(case$id, case$fn, warmup = warmup, runs = runs)
  })

  payload <- list(
    schema_version = "core_v1",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    environment = list(
      r_version = as.character(getRversion()),
      platform = R.version$platform,
      seurat = as.character(packageVersion("Seurat")),
      seurat_object = as.character(packageVersion("SeuratObject"))
    ),
    protocol = list(
      warmup_runs = warmup,
      measured_runs = runs,
      threshold = threshold,
      metrics = c("median_elapsed_sec", "median_allocated_bytes")
    ),
    metrics = metrics
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_json(payload, output_path, pretty = TRUE, auto_unbox = TRUE, null = "null")

  comparison <- NULL
  if (enforce) {
    comparison <- compare_benchmark_to_baseline(
      current_path = output_path,
      baseline_path = baseline_path,
      threshold = threshold
    )
    if (!comparison$ok) {
      stop("Benchmark regression detected:\n", paste(comparison$regressions, collapse = "\n"), call. = FALSE)
    }
  }

  invisible(list(output_path = output_path, comparison = comparison))
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  output_path <- if (length(args) >= 1) args[[1]] else file.path("tests", "benchmarks", "current", "core_v1_metrics.json")
  baseline_path <- if (length(args) >= 2) args[[2]] else file.path("tests", "benchmarks", "baseline", "core_v1_metrics.json")
  threshold <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.05
  enforce <- if (length(args) >= 4) as.logical(args[[4]]) else FALSE

  run_core_benchmarks(
    enforce = enforce,
    output_path = output_path,
    baseline_path = baseline_path,
    threshold = threshold,
    warmup = 1,
    runs = 7
  )
  message("Core benchmark run complete")
}
