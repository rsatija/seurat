#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(Matrix)
  library(Seurat)
  library(SeuratData)
})

load_all('/Users/rahulsatija/seurat_rewrite2', reset = TRUE, quiet = TRUE, recompile = TRUE)

bmcite <- LoadData('bmcite')

message('Available assays: ', paste(Assays(bmcite), collapse = ', '))

run_stdassay_fvf <- function(assay_name, impl, reps = 10L, method = 'vst', nfeatures = 2000L, verbose = FALSE) {
  options(
    Seurat.FindVariableFeatures.StdAssay.impl = impl,
    Seurat.FindVariableFeatures.Assay.impl = impl,
    Seurat.FindVariableFeatures.V3Matrix.impl = impl,
    Seurat.FindVariableFeatures.disp.impl = impl,
    Seurat.FindVariableFeatures.mvp.impl = impl
  )

  times <- numeric(reps)
  baseline <- NULL

  for (i in seq_len(reps)) {
    obj <- bmcite[[assay_name]]
    t <- system.time({
      out <- suppressMessages(suppressWarnings(
        FindVariableFeatures(
          object = obj,
          selection.method = method,
          nfeatures = nfeatures,
          verbose = verbose
        )
      ))
    })['elapsed']
    times[i] <- as.numeric(t)

    if (is.null(baseline)) {
      baseline <- VariableFeatures(out)
    }
  }

  list(times = times, baseline = baseline, object = out)
}

run_matrix_fvf <- function(assay_name, impl, reps = 10L, method = 'vst', nfeatures = 2000L, verbose = FALSE) {
  options(
    Seurat.FindVariableFeatures.StdAssay.impl = impl,
    Seurat.FindVariableFeatures.Assay.impl = impl,
    Seurat.FindVariableFeatures.V3Matrix.impl = impl,
    Seurat.FindVariableFeatures.disp.impl = impl,
    Seurat.FindVariableFeatures.mvp.impl = impl
  )

  mat <- GetAssayData(bmcite[[assay_name]], slot = 'counts')
  times <- numeric(reps)
  baseline <- NULL

  for (i in seq_len(reps)) {
    t <- system.time({
      out <- suppressMessages(suppressWarnings(
        FindVariableFeatures(
          object = mat,
          selection.method = method,
          nfeatures = nfeatures,
          verbose = verbose
        )
      ))
    })['elapsed']
    times[i] <- as.numeric(t)

    if (is.null(baseline)) {
      baseline <- out
    }
  }

  list(times = times, baseline = baseline)
}

check_baseline <- function(a, b, type) {
  if (is.null(a) || is.null(b)) {
    return(list(same = FALSE, reason = 'missing baseline'))
  }
  if (is.character(a) && is.character(b)) {
    same <- identical(a, b)
    return(list(same = same, reason = if (same) 'exact' else 'differences in variable feature order/values'))
  }

  if (is.data.frame(a) && is.data.frame(b)) {
    same <- isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))
    return(list(
      same = same,
      reason = if (same) 'exact' else 'differences in HVF frame values/structure',
      max_abs_diff = if (same) 0 else max(abs(as.numeric(unlist(a)) - as.numeric(unlist(b)), na.rm = TRUE))
    ))
  }

  same <- isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))
  list(same = same, reason = if (same) 'equal' else 'differences', max_abs_diff = if (same) 0 else NA_real_)
}

summary_stats <- function(x) {
  c(
    n = length(x),
    min = min(x),
    q25 = as.numeric(quantile(x, 0.25)),
    median = median(x),
    mean = mean(x),
    q75 = as.numeric(quantile(x, 0.75)),
    max = max(x)
  )
}

assays_to_test <- intersect(c('RNA', 'ADT'), Assays(bmcite))
reps <- 8L

for (assay_name in assays_to_test) {
  message('\n===== Assay: ', assay_name, ' =====')

  # Run direct matrix path (V3Matrix method)
  legacy_mat <- run_matrix_fvf(assay_name = assay_name, impl = 'legacy', reps = reps)
  rewrite_mat <- run_matrix_fvf(assay_name = assay_name, impl = 'rewrite', reps = reps)

  ls_mat <- summary_stats(legacy_mat$times)
  rs_mat <- summary_stats(rewrite_mat$times)
  mat_speedups <- legacy_mat$times / rewrite_mat$times

  cmp_mat <- check_baseline(legacy_mat$baseline, rewrite_mat$baseline, 'matrix')
  message('Matrix-level V3Matrix path:')
  message('  Legacy: ', paste(format(round(legacy_mat$times, 4), nsmall = 4), collapse = ', '))
  message('  Rewrite: ', paste(format(round(rewrite_mat$times, 4), nsmall = 4), collapse = ', '))
  message('  legacy median = ', signif(ls_mat['median'], 4), ', rewrite median = ', signif(rs_mat['median'], 4))
  message('  median speedup (legacy/rewrite) = ', signif(median(mat_speedups), 4),
          '  mean speedup = ', signif(mean(mat_speedups), 4),
          '  p95 speedup = ', signif(quantile(mat_speedups, 0.95), 4))
  message('  output equality: ', cmp_mat$same, ' (', cmp_mat$reason, ')')
  if (!is.null(cmp_mat$max_abs_diff)) {
    message('  max_abs_diff: ', signif(cmp_mat$max_abs_diff, 6))
  }

  # Run StdAssay/Assay path (full object output with VariableFeatures())
  legacy_assay <- run_stdassay_fvf(assay_name = assay_name, impl = 'legacy', reps = reps)
  rewrite_assay <- run_stdassay_fvf(assay_name = assay_name, impl = 'rewrite', reps = reps)

  ls_as <- summary_stats(legacy_assay$times)
  rs_as <- summary_stats(rewrite_assay$times)
  ass_speedups <- legacy_assay$times / rewrite_assay$times

  cmp_as <- check_baseline(legacy_assay$baseline, rewrite_assay$baseline, 'vf')
  message('StdAssay/Assay path:')
  message('  Legacy: ', paste(format(round(legacy_assay$times, 4), nsmall = 4), collapse = ', '))
  message('  Rewrite: ', paste(format(round(rewrite_assay$times, 4), nsmall = 4), collapse = ', '))
  message('  legacy median = ', signif(ls_as['median'], 4), ', rewrite median = ', signif(rs_as['median'], 4))
  message('  median speedup (legacy/rewrite) = ', signif(median(ass_speedups), 4),
          '  mean speedup = ', signif(mean(ass_speedups), 4),
          '  p95 speedup = ', signif(quantile(ass_speedups, 0.95), 4))
  message('  variable feature list equality: ', cmp_as$same, ' (', cmp_as$reason, ')')
}
