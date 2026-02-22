# Step 1 Core Baseline Review

## Summary of changes
- Added a core golden harness with deterministic scenario execution, payload serialization, manifest generation, hash validation, and deep payload comparison.
- Added full serialized golden outputs for 10 core scenarios in `tests/testdata/golden/core_v1/` plus manifest metadata in `tests/testdata/golden/core_v1/manifest.json`.
- Added benchmark harness scripts with runtime/allocation measurement and threshold comparison:
  - `tests/benchmarks/run_core_benchmarks.R`
  - `tests/benchmarks/compare_core_benchmarks.R`
  - baseline metrics at `tests/benchmarks/baseline/core_v1_metrics.json`
- Added tests covering manifest validation, golden verification, determinism checks, and negative-path checks in `tests/testthat/test_core_goldens.R`.
- Added generation entrypoint `tests/scripts/generate_core_goldens.R`.
- Updated packaging ignore rules in `.Rbuildignore` for `diffs/` and `review/` artifacts.

## Compatibility notes
- No exported API names, signatures, defaults, or return types were changed.
- No production behavior rewrite was introduced.
- `CreateSeuratObject` remains reexport-driven and is baseline-validated only through scenarios.
- Baseline scope is v5-focused and limited to the selected core pipelines for this step.

## Unified diff reference
- Full unified patch artifact: `diffs/step1_core_baseline.patch`
- Excerpt:

```diff
diff --git a/.Rbuildignore b/.Rbuildignore
index d203ce14..146edcd4 100644
--- a/.Rbuildignore
+++ b/.Rbuildignore
@@ -13,3 +13,6 @@ CODE_OF_CONDUCT.md
 ^index\.md$
 ^vignettes$
 ^LICENSE\.md$
+
+^diffs$
+^review$
```

## Golden verification and benchmark checks run
- `Rscript tests/scripts/generate_core_goldens.R` (passed)
- `Rscript -e "library(testthat); library(Seurat); testthat::test_file('tests/testthat/test_core_goldens.R')"` (passed; non-failing UMAP method warning observed)
- `Rscript tests/benchmarks/run_core_benchmarks.R tests/benchmarks/baseline/core_v1_metrics.json` (passed)
- `Rscript tests/benchmarks/compare_core_benchmarks.R tests/benchmarks/baseline/core_v1_metrics.json tests/benchmarks/baseline/core_v1_metrics.json 0.05` (passed)

## Notes
- Benchmark gating remains manual artifact review mode in this step; no CI hard-fail benchmark gate was added.
- Golden regeneration has now been performed intentionally for this baseline lock step.
