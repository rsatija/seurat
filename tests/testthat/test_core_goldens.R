context("core-goldens")

test_that("core golden manifest is present and valid", {
  expect_no_error(validate_core_manifest())
})

test_that("core scenarios match serialized goldens", {
  expect_true(verify_core_goldens())
})

test_that("determinism: same seed reproduces core_preprocess_log", {
  spec <- get_core_scenario_specs()[[1]]
  run1 <- run_core_scenario(spec)
  run2 <- run_core_scenario(spec)
  cmp <- compare_core_payload(
    expected = run1$payload,
    observed = run2$payload,
    tolerance = run1$tolerance$numeric
  )
  expect_true(cmp$ok, info = paste(cmp$errors, collapse = "\n"))
})

test_that("determinism: explicit RNG version reproduces core_sctransform_v1", {
  specs <- setNames(get_core_scenario_specs(), vapply(get_core_scenario_specs(), function(x) x$id, character(1)))
  spec <- specs[["core_sctransform_v1"]]
  run1 <- run_core_scenario(spec)
  run2 <- run_core_scenario(spec)
  cmp <- compare_core_payload(
    expected = run1$payload,
    observed = run2$payload,
    tolerance = run1$tolerance$numeric
  )
  expect_true(cmp$ok, info = paste(cmp$errors, collapse = "\n"))
})

test_that("negative: missing golden file fails manifest validation", {
  src <- core_golden_root()
  tmp <- tempfile("core-golden-missing")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  copied_dir <- file.path(tmp, "core_v1")
  dir.create(copied_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(list.files(src, full.names = TRUE), copied_dir, recursive = TRUE)
  expect_true(all(copied))
  manifest <- jsonlite::fromJSON(core_manifest_path(copied_dir), simplifyVector = FALSE)
  file.remove(file.path(copied_dir, manifest$scenarios[[1]]$file))

  expect_error(validate_core_manifest(copied_dir), "Missing golden file")
})

test_that("negative: manifest hash mismatch fails validation", {
  src <- core_golden_root()
  tmp <- tempfile("core-golden-hash")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  copied_dir <- file.path(tmp, "core_v1")
  dir.create(copied_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(list.files(src, full.names = TRUE), copied_dir, recursive = TRUE)
  expect_true(all(copied))
  manifest_path <- core_manifest_path(copied_dir)
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  manifest$scenarios[[1]]$md5 <- "00000000000000000000000000000000"
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null")

  expect_error(validate_core_manifest(copied_dir), "Manifest hash mismatch")
})
