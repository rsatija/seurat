#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(Matrix)
  library(Seurat)
  library(SeuratData)
})

# Load local package source with our edits
load_all('/Users/rahulsatija/seurat_rewrite2', reset = TRUE, quiet = TRUE, recompile = TRUE)

message('Loading bmcite dataset...')
bmcite <- LoadData('bmcite')

adt_assay_name <- 'ADT'
if (!'ADT' %in% Assays(bmcite)) {
  stop('ADT assay not found in bmcite. Available assays: ', paste(Assays(bmcite), collapse = ', '))
}

run_norm <- function(method, impl, verbose = FALSE, margin = 2, scale.factor = 1e4) {
  options(
    Seurat.NormalizeData.clr.impl = impl,
    Seurat.NormalizeData.rc.impl = impl
  )
  t <- system.time({
    out <- suppressMessages(suppressWarnings(
      NormalizeData(
        object = bmcite,
        assay = adt_assay_name,
        normalization.method = method,
        scale.factor = scale.factor,
        margin = margin,
        verbose = verbose
      )
    ))
  })['elapsed']
  expr <- out@assays[[adt_assay_name]]@data
  list(
    elapsed = as.numeric(t),
    expr = expr
  )
}

methods <- c('CLR', 'RC')
results <- list()

for (m in methods) {
  message('\n### Method: ', m)

  # margin choice for CLR; RC is column-wise by design, margin is ignored
  mgn <- if (identical(m, 'CLR')) 1L else 2L

  legacy <- run_norm(method = m, impl = 'legacy', verbose = FALSE, margin = mgn)
  rewrite <- run_norm(method = m, impl = 'rewrite', verbose = FALSE, margin = mgn)

  legacy_mat <- as.matrix(legacy$expr)
  rewrite_mat <- as.matrix(rewrite$expr)

  same_shape <- all(dim(legacy_mat) == dim(rewrite_mat))
  same_dimnames <- identical(dimnames(legacy_mat), dimnames(rewrite_mat))
  max_abs_diff <- if (same_shape) max(abs(legacy_mat - rewrite_mat), na.rm = TRUE) else NA_real_

  results[[m]] <- list(
    legacy_seconds = legacy$elapsed,
    rewrite_seconds = rewrite$elapsed,
    speedup = ifelse(rewrite$elapsed > 0, legacy$elapsed / rewrite$elapsed, NA_real_),
    same_shape = same_shape,
    same_dimnames = same_dimnames,
    max_abs_diff = max_abs_diff,
    all_equal = same_shape && isTRUE(all.equal(legacy_mat, rewrite_mat, tolerance = 1e-10, check.attributes = FALSE))
  )

  message('legacy: ', signif(results[[m]]$legacy_seconds, 4), 's')
  message('rewrite: ', signif(results[[m]]$rewrite_seconds, 4), 's')
  message('speedup: ', ifelse(is.finite(results[[m]]$speedup), signif(results[[m]]$speedup, 4), NA_character_))
  message('max abs diff: ', signif(results[[m]]$max_abs_diff, 4))
  message('all.equal: ', results[[m]]$all_equal)
}

for (m in methods) {
  r <- results[[m]]
  writeLines(sprintf(
    '%s: legacy=%.4fs rewrite=%.4fs speedup=%.3fx same_shape=%s same_dimnames=%s max_abs_diff=%g all_equal=%s',
    m,
    r$legacy_seconds,
    r$rewrite_seconds,
    r$speedup,
    r$same_shape,
    r$same_dimnames,
    r$max_abs_diff,
    r$all_equal
  ))
}
