#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(Matrix)
  library(Seurat)
  library(SeuratData)
})

# Load local package source with edits
load_all('/Users/rahulsatija/seurat_rewrite2', reset = TRUE, quiet = TRUE, recompile = TRUE)

message('Loading bmcite dataset...')
bmcite <- LoadData('bmcite')
adt_assay_name <- 'ADT'
if (!'ADT' %in% Assays(bmcite)) {
  stop('ADT assay not found in bmcite. Available assays: ', paste(Assays(bmcite), collapse = ', '))
}

run_norm <- function(obj, method, impl, margin = 2, scale.factor = 1e4, verbose = FALSE) {
  options(
    Seurat.NormalizeData.clr.impl = impl,
    Seurat.NormalizeData.rc.impl = impl
  )
  t <- system.time({
    out <- suppressMessages(suppressWarnings(
      NormalizeData(
        object = obj,
        assay = adt_assay_name,
        normalization.method = method,
        scale.factor = scale.factor,
        margin = margin,
        verbose = verbose
      )
    ))
  })['elapsed']
  expr <- out@assays[[adt_assay_name]]@data
  list(elapsed = as.numeric(t), expr = expr)
}

# one baseline legacy extraction for shape/naming checks
baseline <- run_norm(bmcite, method = 'CLR', impl = 'legacy', margin = 1)
base_dim <- dim(baseline$expr)
base_dimnames <- dimnames(baseline$expr)

methods <- c('CLR', 'RC')
replicates <- 10
set.seed(1)

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

for (m in methods) {
  message('\n=== Method: ', m, ' ===')
  margin <- if (identical(m, 'CLR')) 1L else 2L

  legacy_times <- numeric(replicates)
  rewrite_times <- numeric(replicates)

  for (i in seq_len(replicates)) {
    legacy <- run_norm(bmcite, method = m, impl = 'legacy', margin = margin)
    rewrite <- run_norm(bmcite, method = m, impl = 'rewrite', margin = margin)

    legacy_times[i] <- legacy$elapsed
    rewrite_times[i] <- rewrite$elapsed

    if (i == 1) {
      legacy_mat <- legacy$expr
      rewrite_mat <- rewrite$expr
      same_shape <- all(dim(legacy_mat) == base_dim)
      same_dimnames <- identical(dimnames(legacy_mat), base_dimnames)
      max_abs_diff <- if (same_shape) max(abs(as.matrix(legacy_mat) - as.matrix(rewrite_mat)), na.rm = TRUE) else NA_real_
      all_equal <- same_shape && isTRUE(all.equal(as.matrix(legacy_mat), as.matrix(rewrite_mat), tolerance = 1e-10, check.attributes = FALSE))
    }

    message(sprintf('Rep %02d: legacy=%.4fs rewrite=%.4fs', i, legacy$elapsed, rewrite$elapsed))
  }

  message('Legacy timings: ', paste(format(round(legacy_times, 4), nsmall = 4), collapse = ', '))
  message('Rewrite timings: ', paste(format(round(rewrite_times, 4), nsmall = 4), collapse = ', '))

  ls <- summary_stats(legacy_times)
  rs <- summary_stats(rewrite_times)
  speedups <- legacy_times / rewrite_times

  message('legacy summary: min=', signif(ls['min'],4), ', q25=', signif(ls['q25'],4), ', med=', signif(ls['median'],4), ', mean=', signif(ls['mean'],4), ', q75=', signif(ls['q75'],4), ', max=', signif(ls['max'],4))
  message('rewrite summary: min=', signif(rs['min'],4), ', q25=', signif(rs['q25'],4), ', med=', signif(rs['median'],4), ', mean=', signif(rs['mean'],4), ', q75=', signif(rs['q75'],4), ', max=', signif(rs['max'],4))
  message('median speedup (legacy/rewrite)=', signif(median(speedups), 4),
          '  mean speedup=', signif(mean(speedups), 4),
          '  p95 speedup=', signif(quantile(speedups, 0.95), 4))
  message('shape match: ', same_shape, '  dimnames match: ', same_dimnames, '  max_abs_diff: ', max_abs_diff, '  all_equal: ', all_equal)
}
