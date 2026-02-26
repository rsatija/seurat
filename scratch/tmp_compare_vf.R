library(devtools)
library(Seurat)
library(SeuratData)

load_all('/Users/rahulsatija/seurat_rewrite2', reset=TRUE, quiet=TRUE, recompile=TRUE)

bmcite <- LoadData('bmcite')
o <- bmcite[['RNA']]

legacy_opts <- list(
  Seurat.FindVariableFeatures.StdAssay.impl='legacy',
  Seurat.FindVariableFeatures.Assay.impl='legacy',
  Seurat.FindVariableFeatures.V3Matrix.impl='legacy',
  Seurat.FindVariableFeatures.disp.impl='legacy',
  Seurat.FindVariableFeatures.mvp.impl='legacy'
)
rewrite_opts <- list(
  Seurat.FindVariableFeatures.StdAssay.impl='rewrite',
  Seurat.FindVariableFeatures.Assay.impl='rewrite',
  Seurat.FindVariableFeatures.V3Matrix.impl='rewrite',
  Seurat.FindVariableFeatures.disp.impl='rewrite',
  Seurat.FindVariableFeatures.mvp.impl='rewrite'
)

options(legacy_opts)
l <- FindVariableFeatures(o, selection.method='vst', nfeatures=2000, verbose=FALSE)
options(rewrite_opts)
r <- FindVariableFeatures(o, selection.method='vst', nfeatures=2000, verbose=FALSE)

lv <- VariableFeatures(l)
rv <- VariableFeatures(r)

cat('same_exact=', identical(lv, rv), '\n')
cat('set_equal=', setequal(lv, rv), ' intersection=', length(intersect(lv, rv)), '\n')

diff_idx <- which(lv != rv)
cat('diff_count=', length(diff_idx), '\n')
if (length(diff_idx) > 0) {
  n_show <- min(20, length(diff_idx))
  sel <- diff_idx[seq_len(n_show)]
  cat('first positions with diffs=', paste(sel, collapse = ','), '\n')
  cat('legacy at positions=', paste(lv[sel], collapse = ','), '\n')
  cat('rewrite at positions=', paste(rv[sel], collapse = ','), '\n')
}
