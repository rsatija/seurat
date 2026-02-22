find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(current, "NAMESPACE")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop("Unable to locate repository root from: ", start, call. = FALSE)
}

core_golden_root <- function() {
  file.path(find_repo_root(), "tests", "testdata", "golden", "core_v1")
}

core_manifest_path <- function(golden_dir = core_golden_root()) {
  file.path(golden_dir, "manifest.json")
}

core_fixture_pbmc_raw <- function() {
  pbmc.file <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")
  if (!nzchar(pbmc.file)) {
    stop("Unable to locate Seurat extdata/pbmc_raw.txt", call. = FALSE)
  }
  as.sparse(as.matrix(read.table(pbmc.file, sep = "\t", row.names = 1)))
}

md5_of_file <- function(path) {
  unname(tools::md5sum(path))
}

md5_of_object <- function(object) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  md5_of_file(tmp)
}

matrix_signature <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(list(
      class = class(x),
      dim = dim(x),
      dimnames = dimnames(x),
      i = x@i,
      p = x@p,
      x = x@x
    ))
  }
  if (is.matrix(x)) {
    return(list(
      class = class(x),
      dim = dim(x),
      dimnames = dimnames(x),
      x = unclass(x)
    ))
  }
  if (is.data.frame(x)) {
    return(list(
      class = class(x),
      dim = c(nrow(x), ncol(x)),
      dimnames = list(rownames(x), colnames(x)),
      x = x
    ))
  }
  stop("Unsupported matrix type for signature: ", paste(class(x), collapse = ", "), call. = FALSE)
}

assay_signature <- function(assay) {
  layers <- tryCatch(Layers(assay), error = function(e) character())
  layer_payload <- list()
  if (length(layers) > 0) {
    for (layer_name in layers) {
      layer_payload[[layer_name]] <- matrix_signature(LayerData(assay, layer = layer_name))
    }
  }
  list(
    class = class(assay),
    layers = layers,
    layer_payload = layer_payload,
    var_features = VariableFeatures(assay),
    meta_features = if ("meta.features" %in% slotNames(assay)) assay@meta.features else NULL
  )
}

dimreduc_signature <- function(dr) {
  list(
    class = class(dr),
    key = Key(dr),
    embeddings = matrix_signature(Embeddings(dr)),
    loadings = matrix_signature(Loadings(dr)),
    stdev = Stdev(dr),
    misc = dr@misc
  )
}

seurat_signature <- function(object) {
  assay_names <- Assays(object)
  assays <- setNames(lapply(assay_names, function(a) assay_signature(object[[a]])), assay_names)

  reduction_names <- Reductions(object)
  reductions <- setNames(lapply(reduction_names, function(r) dimreduc_signature(object[[r]])), reduction_names)

  graph_names <- Graphs(object)
  graphs <- setNames(lapply(graph_names, function(g) matrix_signature(object[[g]])), graph_names)

  list(
    class = class(object),
    default_assay = DefaultAssay(object),
    assays = assay_names,
    assay_payload = assays,
    reductions = reduction_names,
    reduction_payload = reductions,
    graphs = graph_names,
    graph_payload = graphs,
    cells = Cells(object),
    meta_data = object[[]],
    active_ident = as.character(Idents(object)),
    active_ident_levels = levels(Idents(object)),
    commands = names(object@commands)
  )
}

anchorset_signature <- function(anchors) {
  list(
    class = class(anchors),
    anchors = matrix_signature(anchors@anchors),
    anchor_features = anchors@anchor.features,
    reference_cells = anchors@reference.cells,
    query_cells = anchors@query.cells,
    reference_objects = anchors@reference.objects,
    offsets = anchors@offsets,
    neighbors = anchors@neighbors,
    object_1 = seurat_signature(anchors@object.list[[1]])
  )
}

integration_object_signature <- function(object) {
  sig <- seurat_signature(object)
  integration_tool <- Tool(object = object, slot = "Integration")
  tool_slots <- slotNames(integration_tool)
  sig$integration_tool <- list(
    class = class(integration_tool),
    slots = tool_slots,
    sample_tree = if ("sample.tree" %in% tool_slots) integration_tool@sample.tree else NULL,
    preserve_order = if ("preserve.order" %in% tool_slots) integration_tool@preserve.order else NULL,
    weight = if ("weight" %in% tool_slots) integration_tool@weight else NULL
  )
  sig
}

transfer_output_signature <- function(preds_standard, pred_assay) {
  list(
    preds_standard = preds_standard,
    pred_assay = assay_signature(pred_assay)
  )
}

create_transfer_fixture <- function() {
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
  list(ref = ref, query = query)
}

create_integrate_fixture <- function() {
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

  list(ref = ref, query.list = query.list)
}

core_fixture_fingerprints <- function() {
  pbmc_path <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")
  list(
    pbmc_raw_txt_md5 = md5_of_file(pbmc_path),
    pbmc_small_md5 = md5_of_object(pbmc_small)
  )
}

scenario_meta <- function(id, seed, rng_version = NULL) {
  list(
    scenario_id = id,
    seed = seed,
    rng_version = rng_version,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_fingerprint = list(
      seurat = as.character(packageVersion("Seurat")),
      seurat_object = as.character(packageVersion("SeuratObject")),
      testthat = as.character(utils::packageVersion("testthat"))
    ),
    session_info = capture.output(utils::sessionInfo()),
    fixture_fingerprints = core_fixture_fingerprints()
  )
}

scenario_core_preprocess_log <- function() {
  pbmc <- core_fixture_pbmc_raw()
  object <- CreateSeuratObject(counts = pbmc)
  object <- NormalizeData(object, verbose = FALSE, scale.factor = 1e6)
  object <- FindVariableFeatures(object, verbose = FALSE)
  object <- ScaleData(object, verbose = FALSE)
  seurat_signature(object)
}

scenario_core_preprocess_clr_rc <- function() {
  pbmc <- core_fixture_pbmc_raw()
  clr <- NormalizeData(pbmc, normalization.method = "CLR", verbose = FALSE)
  rc <- NormalizeData(pbmc, normalization.method = "RC", verbose = FALSE)
  list(
    clr = matrix_signature(clr),
    rc = matrix_signature(rc)
  )
}

scenario_core_pca_umap_neighbors_clusters <- function() {
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
}

scenario_core_sctransform_v1 <- function() {
  pbmc <- core_fixture_pbmc_raw()
  object <- CreateSeuratObject(counts = pbmc)
  object <- NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(SCTransform(object, verbose = FALSE, vst.flavor = "v1", seed.use = 1448145))

  list(
    sct_assay = assay_signature(object[["SCT"]]),
    sct_results = SCTResults(object = object, assay = "SCT", slot = "feature.attributes")
  )
}

scenario_core_transfer_anchors_default <- function() {
  fixture <- create_transfer_fixture()
  anchors <- FindTransferAnchors(reference = fixture$ref, query = fixture$query, k.filter = 50)
  anchorset_signature(anchors)
}

scenario_core_integrate_anchors_twoobj <- function() {
  fixture <- create_integrate_fixture()
  anchors <- suppressMessages(
    suppressWarnings(
      FindIntegrationAnchors(
        object.list = c(fixture$ref, fixture$query.list[[1]]),
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
}

scenario_core_integrate_data_twoobj <- function() {
  fixture <- create_integrate_fixture()
  anchors <- suppressMessages(
    suppressWarnings(
      FindIntegrationAnchors(
        object.list = c(fixture$ref, fixture$query.list[[1]]),
        k.filter = NA,
        verbose = FALSE
      )
    )
  )
  integrated <- suppressWarnings(IntegrateData(anchorset = anchors, k.weight = 50, verbose = FALSE))
  integration_object_signature(integrated)
}

scenario_core_integrate_data_threeobj <- function() {
  fixture <- create_integrate_fixture()
  anchors <- suppressMessages(
    suppressWarnings(
      FindIntegrationAnchors(
        object.list = c(fixture$ref, fixture$query.list),
        k.filter = NA,
        verbose = FALSE
      )
    )
  )
  integrated <- suppressWarnings(IntegrateData(anchorset = anchors, k.weight = 25, verbose = FALSE))
  integration_object_signature(integrated)
}

scenario_core_transfer_data_default <- function() {
  fixture <- create_transfer_fixture()
  anchors <- FindTransferAnchors(reference = fixture$ref, query = fixture$query, k.filter = 50)
  preds.standard <- TransferData(anchorset = anchors, refdata = fixture$ref$RNA_snn_res.1, verbose = FALSE)
  pred.assay <- TransferData(anchorset = anchors, refdata = GetAssayData(fixture$ref[["RNA"]]), verbose = FALSE)
  transfer_output_signature(preds.standard, pred.assay)
}

scenario_core_mapquery_default <- function() {
  fixture <- create_transfer_fixture()
  ref <- fixture$ref
  query <- fixture$query

  ref <- ScaleData(ref, verbose = FALSE)
  ref <- suppressWarnings(RunPCA(ref, npcs = 30, verbose = FALSE))
  ref <- suppressWarnings(RunUMAP(ref, reduction = "pca", dims = 1:30, return.model = TRUE, seed.use = 42, verbose = FALSE))

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
}

core_scenario_specs <- function() {
  list(
    list(id = "core_preprocess_log", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_preprocess_log),
    list(id = "core_preprocess_clr_rc", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_preprocess_clr_rc),
    list(id = "core_pca_umap_neighbors_clusters", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-5), runner = scenario_core_pca_umap_neighbors_clusters),
    list(id = "core_sctransform_v1", seed = 1448145L, rng_version = "3.5.0", tolerance = list(numeric = 1e-6), runner = scenario_core_sctransform_v1),
    list(id = "core_transfer_anchors_default", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-7), runner = scenario_core_transfer_anchors_default),
    list(id = "core_integrate_anchors_twoobj", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-7), runner = scenario_core_integrate_anchors_twoobj),
    list(id = "core_integrate_data_twoobj", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_integrate_data_twoobj),
    list(id = "core_integrate_data_threeobj", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_integrate_data_threeobj),
    list(id = "core_transfer_data_default", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_transfer_data_default),
    list(id = "core_mapquery_default", seed = 42L, rng_version = NULL, tolerance = list(numeric = 1e-6), runner = scenario_core_mapquery_default)
  )
}

get_core_scenario_specs <- function() {
  core_scenario_specs()
}

run_core_scenario <- function(spec) {
  stopifnot(is.list(spec), !is.null(spec$id), !is.null(spec$runner), is.function(spec$runner))

  old_rng_kind <- RNGkind()
  if (!is.null(spec$rng_version)) {
    suppressWarnings(RNGversion(spec$rng_version))
  }
  on.exit({
    if (!is.null(spec$rng_version)) {
      do.call(RNGkind, as.list(old_rng_kind))
    }
  }, add = TRUE)

  set.seed(spec$seed)
  payload <- spec$runner()

  list(
    meta = scenario_meta(id = spec$id, seed = spec$seed, rng_version = spec$rng_version),
    payload = payload,
    tolerance = spec$tolerance
  )
}

write_manifest <- function(manifest, manifest_path) {
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

generate_core_goldens <- function(output_dir = core_golden_root()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- get_core_scenario_specs()

  scenario_entries <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    result <- run_core_scenario(spec)
    out_file <- file.path(output_dir, paste0(spec$id, ".rds"))
    saveRDS(result, out_file, version = 3)

    scenario_entries[[i]] <- list(
      scenario_id = spec$id,
      file = basename(out_file),
      md5 = md5_of_file(out_file),
      seed = spec$seed,
      rng_version = spec$rng_version,
      tolerance = spec$tolerance
    )
  }

  manifest <- list(
    schema_version = "core_v1",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    generator = list(
      r_version = as.character(getRversion()),
      seurat = as.character(packageVersion("Seurat")),
      seurat_object = as.character(packageVersion("SeuratObject"))
    ),
    scenarios = scenario_entries
  )

  write_manifest(manifest, core_manifest_path(output_dir))
  invisible(manifest)
}

validate_core_manifest <- function(golden_dir = core_golden_root()) {
  manifest_file <- core_manifest_path(golden_dir)
  if (!file.exists(manifest_file)) {
    stop("Core golden manifest not found: ", manifest_file, call. = FALSE)
  }

  manifest <- jsonlite::fromJSON(manifest_file, simplifyVector = FALSE)
  if (is.null(manifest$scenarios) || length(manifest$scenarios) == 0) {
    stop("Core golden manifest has no scenarios", call. = FALSE)
  }

  declared_ids <- vapply(manifest$scenarios, function(x) x$scenario_id, character(1))
  expected_ids <- vapply(get_core_scenario_specs(), function(x) x$id, character(1))
  if (!identical(sort(declared_ids), sort(expected_ids))) {
    stop("Scenario IDs in manifest do not match expected core scenario matrix", call. = FALSE)
  }

  for (entry in manifest$scenarios) {
    file_path <- file.path(golden_dir, entry$file)
    if (!file.exists(file_path)) {
      stop("Missing golden file: ", entry$file, call. = FALSE)
    }

    observed <- md5_of_file(file_path)
    if (!identical(observed, entry$md5)) {
      stop(
        "Manifest hash mismatch for ",
        entry$scenario_id,
        ": expected ",
        entry$md5,
        ", observed ",
        observed,
        call. = FALSE
      )
    }
  }

  manifest
}

collect_compare_errors <- function(expected, observed, path = "$", tolerance = 0) {
  errors <- character()

  if (is.null(expected) || is.null(observed)) {
    if (!identical(expected, observed)) {
      return(sprintf("%s: NULL mismatch", path))
    }
    return(errors)
  }

  if (!identical(class(expected), class(observed))) {
    return(sprintf(
      "%s: class mismatch expected [%s] observed [%s]",
      path,
      paste(class(expected), collapse = ","),
      paste(class(observed), collapse = ",")
    ))
  }

  if (is.list(expected) && !is.data.frame(expected)) {
    if (!identical(names(expected), names(observed))) {
      return(sprintf("%s: list names mismatch", path))
    }
    if (!identical(length(expected), length(observed))) {
      return(sprintf("%s: list length mismatch", path))
    }
    for (nm in names(expected)) {
      child <- collect_compare_errors(expected[[nm]], observed[[nm]], paste0(path, "$", nm), tolerance)
      errors <- c(errors, child)
    }
    return(errors)
  }

  if (is.data.frame(expected)) {
    if (!identical(colnames(expected), colnames(observed))) {
      return(sprintf("%s: data.frame column names mismatch", path))
    }
    if (!identical(rownames(expected), rownames(observed))) {
      return(sprintf("%s: data.frame row names mismatch", path))
    }
    if (!identical(vapply(expected, class, character(1)), vapply(observed, class, character(1)))) {
      return(sprintf("%s: data.frame column class mismatch", path))
    }
    for (nm in colnames(expected)) {
      child <- collect_compare_errors(expected[[nm]], observed[[nm]], paste0(path, "$", nm), tolerance)
      errors <- c(errors, child)
    }
    return(errors)
  }

  if (is.matrix(expected) && is.numeric(expected)) {
    if (!identical(dim(expected), dim(observed))) {
      return(sprintf("%s: matrix dim mismatch", path))
    }
    if (!identical(dimnames(expected), dimnames(observed))) {
      return(sprintf("%s: matrix dimnames mismatch", path))
    }
    if (length(expected) == 0) {
      return(errors)
    }

    delta <- abs(expected - observed)
    delta[is.na(delta)] <- 0
    if (any(delta > tolerance)) {
      return(sprintf("%s: numeric matrix mismatch max_delta=%g tolerance=%g", path, max(delta), tolerance))
    }
    return(errors)
  }

  if (is.numeric(expected)) {
    if (!identical(length(expected), length(observed))) {
      return(sprintf("%s: numeric length mismatch", path))
    }
    if (!identical(names(expected), names(observed))) {
      return(sprintf("%s: numeric names mismatch", path))
    }
    if (length(expected) == 0) {
      return(errors)
    }

    delta <- abs(expected - observed)
    delta[is.na(delta)] <- 0
    if (any(delta > tolerance)) {
      return(sprintf("%s: numeric mismatch max_delta=%g tolerance=%g", path, max(delta), tolerance))
    }
    return(errors)
  }

  if (!identical(expected, observed)) {
    return(sprintf("%s: value mismatch", path))
  }

  errors
}

compare_core_payload <- function(expected, observed, tolerance = 0) {
  errs <- collect_compare_errors(expected, observed, tolerance = tolerance)
  list(ok = length(errs) == 0, errors = errs)
}

verify_core_goldens <- function(golden_dir = core_golden_root()) {
  manifest <- validate_core_manifest(golden_dir)
  entries <- setNames(manifest$scenarios, vapply(manifest$scenarios, function(x) x$scenario_id, character(1)))

  all_errors <- character()
  for (spec in get_core_scenario_specs()) {
    entry <- entries[[spec$id]]
    golden_file <- file.path(golden_dir, entry$file)
    expected <- readRDS(golden_file)
    observed <- run_core_scenario(spec)

    tol <- expected$tolerance$numeric
    if (is.null(tol)) {
      tol <- 0
    }

    payload_cmp <- compare_core_payload(expected$payload, observed$payload, tolerance = tol)
    if (!payload_cmp$ok) {
      all_errors <- c(
        all_errors,
        paste0("[", spec$id, "] ", payload_cmp$errors)
      )
    }

    meta_fields <- c("scenario_id", "seed", "rng_version", "fixture_fingerprints")
    for (field in meta_fields) {
      if (!identical(expected$meta[[field]], observed$meta[[field]])) {
        all_errors <- c(all_errors, paste0("[", spec$id, "] meta mismatch for field `", field, "`"))
      }
    }
  }

  if (length(all_errors) > 0) {
    stop(
      "Core golden verification failed:\n",
      paste(all_errors, collapse = "\n"),
      call. = FALSE
    )
  }

  TRUE
}
