library(devtools)
library(Seurat)
library(SeuratData)

load_all('/Users/rahulsatija/seurat_rewrite2', reset=TRUE, quiet=TRUE, recompile=TRUE)

bmcite <- LoadData('bmcite')
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
l_assay <- FindVariableFeatures(bmcite[['RNA']], selection.method='vst', nfeatures=2000, verbose=FALSE)
options(rewrite_opts)
r_assay <- FindVariableFeatures(bmcite[['RNA']], selection.method='vst', nfeatures=2000, verbose=FALSE)

l_hvf <- l_assay[['vst.variable_features']]
r_hvf <- r_assay[['vst.variable_features']]

# Align by rownames
common <- intersect(rownames(l_hvf), rownames(r_hvf))
# Compare whether legacy top scores for differing positions are tied
cat('legacy top 2050 score range:\n')
cat(paste(format(head(na.omit(l_hvf[1:2050,'vst.variance.standardized']), digits=12), collapse=' , '), '\n')
cat('---\n')

l_top <- l_hvf[1:1300, c('vst.variance.standardized')]
r_top <- r_hvf[1:1300, c('vst.variance.standardized')]

# report first 12 differences in ranks
diff_idx <- which(l_top != r_top)
cat('diff_idx_in_top=', length(diff_idx), '\n')
if (length(diff_idx) > 0) {
  idx <- diff_idx[1:min(20, length(diff_idx))]
  cat('first diff rank positions=', paste(idx, collapse=','), '\n')
  for (i in idx) {
    cat('rank', i, ':', rownames(l_top)[i], l_top[i,1], '|||', rownames(r_top)[i], r_top[i,1], '\n')
  }
}

# direct pairwise score compare for swapped features
features <- c('RP11-301G19.1','RP11-879F14.2','LINC01474','LRRC4C','TTC36','GALNT8','CABP1','CTC-297N7.1')
cat('selected feature scores\n')
for (f in features) {
  if (f %in% rownames(l_hvf) && f %in% rownames(r_hvf)) {
    cat(f, as.numeric(l_hvf[f,'vst.variance.standardized']), as.numeric(r_hvf[f,'vst.variance.standardized']), '\n')
  }
}
