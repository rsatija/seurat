# Core Golden Baseline (`core_v1`)

This directory contains serialized golden outputs for the 10 core baseline scenarios used by the Step 1 lock.

- Manifest: `manifest.json`
- Scenario files: `<scenario_id>.rds`
- Source of truth for scenario definitions: `tests/testthat/helper_core_golden_utils.R`

## Generate / Verify

Generate all goldens:

```sh
Rscript tests/scripts/generate_core_goldens.R
```

Verify against current implementation:

```sh
Rscript -e "library(testthat); library(Seurat); testthat::test_file('tests/testthat/test_core_goldens.R')"
```

## Common Runner Wrapper

Every scenario is executed through `run_core_scenario(spec)` with:

```r
if (!is.null(spec$rng_version)) suppressWarnings(RNGversion(spec$rng_version))
set.seed(spec$seed)
payload <- spec$runner()
```

## Scenario Matrix and Commands

### 1) `core_preprocess_log`

```r
pbmc <- core_fixture_pbmc_raw()
object <- CreateSeuratObject(counts = pbmc)
object <- NormalizeData(object, verbose = FALSE, scale.factor = 1e6)
object <- FindVariableFeatures(object, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
seurat_signature(object)
```

### 2) `core_preprocess_clr_rc`

```r
pbmc <- core_fixture_pbmc_raw()
clr <- NormalizeData(pbmc, normalization.method = "CLR", verbose = FALSE)
rc <- NormalizeData(pbmc, normalization.method = "RC", verbose = FALSE)
list(
  clr = matrix_signature(clr),
  rc = matrix_signature(rc)
)
```

### 3) `core_pca_umap_neighbors_clusters`

```r
pbmc <- core_fixture_pbmc_raw()
assay <- CreateAssay5Object(pbmc)
object <- CreateSeuratObject(assay)
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, npcs = 10, verbose = FALSE)
object <- FindNeighbors(object, k.param = 10, dims = 1:10, verbose = FALSE)
object <- FindClusters(object, random.seed = 42, verbose = FALSE)
object <- RunUMAP(object, dims = 1:10, seed.use = 42, verbose = FALSE)
list(
  seurat = seurat_signature(object),
  clusters = list(
    values = as.character(object$seurat_clusters),
    levels = levels(object$seurat_clusters)
  )
)
```

### 4) `core_sctransform_v1`

```r
pbmc <- core_fixture_pbmc_raw()
object <- CreateSeuratObject(counts = pbmc)
object <- NormalizeData(object, verbose = FALSE)
object <- suppressWarnings(
  SCTransform(object, verbose = FALSE, vst.flavor = "v1", seed.use = 1448145)
)
list(
  sct_assay = assay_signature(object[["SCT"]]),
  sct_results = SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
)
```

### 5) `core_transfer_anchors_default`

`create_transfer_fixture()` setup commands:

```r
set.seed(42)
ref <- suppressWarnings(UpdateSeuratObject(pbmc_small))
query <- CreateSeuratObject(
  counts = as.sparse(
    GetAssayData(ref[["RNA"]], layer = "counts") +
      rpois(n = ncol(ref), lambda = 1)
  )
)
query <- NormalizeData(query, verbose = FALSE)
query <- FindVariableFeatures(query, verbose = FALSE, nfeatures = 100)
ref <- FindVariableFeatures(ref, verbose = FALSE, nfeatures = 100)
```

Scenario commands:

```r
anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50)
anchorset_signature(anchors)
```

### 6) `core_integrate_anchors_twoobj`

`create_integrate_fixture()` setup commands:

```r
set.seed(42)
ref <- suppressWarnings(UpdateSeuratObject(pbmc_small))
ref <- FindVariableFeatures(ref, verbose = FALSE, nfeatures = 100)
query <- CreateSeuratObject(
  counts = as.sparse(
    GetAssayData(ref[["RNA"]], layer = "counts") +
      rpois(n = ncol(ref), lambda = 1)
  )
)
query2 <- CreateSeuratObject(
  counts = as.sparse(
    LayerData(ref[["RNA"]], layer = "counts")[, 1:40] +
      rpois(n = ncol(ref), lambda = 1)
  )
)
query.list <- list(query, query2)
query.list <- lapply(query.list, NormalizeData, verbose = FALSE)
query.list <- lapply(query.list, FindVariableFeatures, verbose = FALSE, nfeatures = 100)
query.list <- lapply(query.list, ScaleData, verbose = FALSE)
query.list <- suppressWarnings(lapply(query.list, RunPCA, verbose = FALSE, npcs = 20))
```

Scenario commands:

```r
anchors <- suppressMessages(
  suppressWarnings(
    FindIntegrationAnchors(
      object.list = c(ref, query.list[[1]]),
      k.filter = NA,
      verbose = FALSE
    )
  )
)
list(
  class = class(anchors),
  anchors = matrix_signature(anchors@anchors),
  anchor_features = anchors@anchor.features,
  reference_cells = anchors@reference.cells,
  query_cells = anchors@query.cells,
  reference_objects = anchors@reference.objects,
  offsets = anchors@offsets
)
```

### 7) `core_integrate_data_twoobj`

Uses the same `create_integrate_fixture()` setup as scenario 6.

Scenario commands:

```r
anchors <- suppressMessages(
  suppressWarnings(
    FindIntegrationAnchors(
      object.list = c(ref, query.list[[1]]),
      k.filter = NA,
      verbose = FALSE
    )
  )
)
integrated <- suppressWarnings(IntegrateData(anchorset = anchors, k.weight = 50, verbose = FALSE))
integration_object_signature(integrated)
```

### 8) `core_integrate_data_threeobj`

Uses the same `create_integrate_fixture()` setup as scenario 6.

Scenario commands:

```r
anchors <- suppressMessages(
  suppressWarnings(
    FindIntegrationAnchors(
      object.list = c(ref, query.list),
      k.filter = NA,
      verbose = FALSE
    )
  )
)
integrated <- suppressWarnings(IntegrateData(anchorset = anchors, k.weight = 25, verbose = FALSE))
integration_object_signature(integrated)
```

### 9) `core_transfer_data_default`

Uses the same `create_transfer_fixture()` setup as scenario 5.

Scenario commands:

```r
anchors <- FindTransferAnchors(reference = ref, query = query, k.filter = 50)
preds.standard <- TransferData(anchorset = anchors, refdata = ref$RNA_snn_res.1, verbose = FALSE)
pred.assay <- TransferData(anchorset = anchors, refdata = GetAssayData(ref[["RNA"]]), verbose = FALSE)
transfer_output_signature(preds.standard, pred.assay)
```

### 10) `core_mapquery_default`

Uses the same `create_transfer_fixture()` setup as scenario 5.

Scenario commands:

```r
ref <- ScaleData(ref, verbose = FALSE)
ref <- suppressWarnings(RunPCA(ref, npcs = 30, verbose = FALSE))
ref <- suppressWarnings(
  RunUMAP(ref, reduction = "pca", dims = 1:30, return.model = TRUE, seed.use = 42, verbose = FALSE)
)
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  dims = 1:30,
  reference.reduction = "pca",
  k.filter = 50
)
mapped <- suppressWarnings(
  MapQuery(
    anchorset = anchors,
    query = query,
    reference = ref,
    refdata = list(predicted.id = ref$RNA_snn_res.1),
    reference.reduction = "pca",
    reduction.model = "umap",
    verbose = FALSE
  )
)
predicted_cols <- grep("^predicted\\.", colnames(mapped[[]]), value = TRUE)
list(
  seurat = seurat_signature(mapped),
  metadata_columns = colnames(mapped[[]]),
  predicted_metadata = if (length(predicted_cols) > 0) mapped[[]][, predicted_cols, drop = FALSE] else NULL
)
```
