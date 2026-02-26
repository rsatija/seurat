library(devtools)
library(Matrix)
library(Seurat)
library(SeuratData)

load_all('.', reset = TRUE, quiet = TRUE, recompile = TRUE)

bmcite <- LoadData('bmcite')
rna <- GetAssayData(bmcite, assay = 'RNA', layer = 'counts')

sparse_equal <- function(x, y) {
  if (!identical(dim(x), dim(y))) return(FALSE)
  if (!identical(class(x), class(y))) return(FALSE)
  if (!identical(slot(x, 'p'), slot(y, 'p'))) return(FALSE)
  if (!identical(slot(x, 'i'), slot(y, 'i'))) return(FALSE)
  max(abs(slot(x, 'x') - slot(y, 'x'))) < 1e-10
}

run_norm <- function(data, impl, reps = 10L, scale.factor = 1e4) {
  options(Seurat.NormalizeData.sparse_normalize.impl = impl)
  times <- numeric(reps)
  baseline <- NULL

  for (i in seq_len(reps)) {
    t <- system.time({
      out <- suppressWarnings(suppressMessages(
        NormalizeData(
          object = data,
          normalization.method = 'LogNormalize',
          cmargin = 2L,
          scale.factor = scale.factor,
          verbose = FALSE
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

legacy <- run_norm(rna, 'legacy', reps = 10L)
rewrite <- run_norm(rna, 'rewrite', reps = 10L)

same <- sparse_equal(legacy$baseline, rewrite$baseline)
max_abs_diff <- if (same) 0 else max(abs(slot(legacy$baseline, 'x') - slot(rewrite$baseline, 'x')))

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

ls <- summary_stats(legacy$times)
rs <- summary_stats(rewrite$times)
speedups <- legacy$times / rewrite$times

cat('Legacy timings:', paste(format(round(legacy$times, 4), nsmall = 4), collapse = ', '), '\n')
cat('Rewrite timings:', paste(format(round(rewrite$times, 4), nsmall = 4), collapse = ', '), '\n')
cat('legacy summary: min=', signif(ls['min'],4), ', q25=', signif(ls['q25'],4), ', med=', signif(ls['median'],4), ', mean=', signif(ls['mean'],4), ', q75=', signif(ls['q75'],4), ', max=', signif(ls['max'],4), '\n')
cat('rewrite summary: min=', signif(rs['min'],4), ', q25=', signif(rs['q25'],4), ', med=', signif(rs['median'],4), ', mean=', signif(rs['mean'],4), ', q75=', signif(rs['q75'],4), ', max=', signif(rs['max'],4), '\n')
cat('median speedup (legacy/rewrite)=', signif(median(speedups), 4),
    '  mean speedup=', signif(mean(speedups), 4),
    '  p95 speedup=', signif(quantile(speedups, 0.95), 4), '\n')
cat('exact match (sparse structure and values):', same, '\n')
cat('max_abs_diff:', max_abs_diff, '\n')

