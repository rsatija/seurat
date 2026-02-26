# NormalizeData Rewrite Changes Summary

## Scope

This summary covers changes made for the NormalizeData rewrite effort in
`/Users/rahulsatija/seurat_rewrite2`, with a focus on speed, reproducibility,
and clean dispatch structure.

## Files changed

- `/Users/rahulsatija/seurat_rewrite2/src/data_manipulation.cpp`
- `/Users/rahulsatija/seurat_rewrite2/src/data_manipulation.h`
- `/Users/rahulsatija/seurat_rewrite2/src/RcppExports.cpp`
- `/Users/rahulsatija/seurat_rewrite2/R/RcppExports.R`
- `/Users/rahulsatija/seurat_rewrite2/R/preprocessing.R`
- `/Users/rahulsatija/seurat_rewrite2/R/preprocessing5.R`

## What changed in `NormalizeData` path

### 1) Added/expanded C++ normalization kernels

In `/Users/rahulsatija/seurat_rewrite2/src/data_manipulation.cpp`:

- `LogNorm()` now performs sparse column scaling and log1p in C++.
- `CLRNorm()` now performs CLR in C++ for margin 1 (row-wise) and margin 2 (column-wise).
- `RelativeCountsNorm()` was added for RC in C++.

#### Rationale

- Move heavy element-wise work off R loops into compiled C++.
- Preserve sparse storage and avoid dense conversion.
- Reduce temporary allocation and extra traversal operations in old R paths.

#### Implementation notes

- `LogNorm`:
  - Replaced `data.transpose() * ones` style sum approach with raw column traversal using sparse internals.
- `CLRNorm`:
 - Initial rewrite used transpose + two-pass per-column strategy.
 - Improved to avoid transpose for row-wise CLR by accumulating row statistics directly from non-zero entries.
- `RelativeCountsNorm`:
  - Uses direct traversal of each column’s non-zero block and rescales entries.

### 2) Exposed rewrite dispatch helpers

In `/Users/rahulsatija/seurat_rewrite2/R/preprocessing.R` and `/Users/rahulsatija/seurat_rewrite2/R/preprocessing5.R`:

- `.CLRNormalize()`:
  - Added implementation selector controlled by:
    - `options("Seurat.NormalizeData.clr.impl")` (`legacy|rewrite`)
  - Default remains legacy with optional rewrite path for sparse `dgCMatrix`.
  - Uses C++ `CLRNorm()` when rewrite is selected.

- `.RCNormalize()`:
  - Initially implemented dispatch via `options("Seurat.NormalizeData.rc.impl")`.
  - Finalized to always call legacy path (decision: rewrite is currently slower).
  - Keeps the rewrite wrapper present but not used in normal execution.

### 3) V3Matrix normalization path updated

In `NormalizeData.V3Matrix` (`/Users/rahulsatija/seurat_rewrite2/R/preprocessing.R`):

- Rewired `CLR`/`RC` branches to use internal helpers (`.CLRNormalize`, `.RCNormalize`)
  rather than custom raw `apply` closures.
- Preserved existing parallel semantics for `LogNormalize` and non-rewrite cases.
- Added logic to disable V3Matrix parallel chunking when CLR/RC rewrite is active
  (avoid orchestration overhead that dominated small/medium workloads).

### 4) Legacy fallback behavior and reproducibility

- Non-`dgCMatrix` inputs to rewrite-sensitive functions fall back to legacy
  implementations in both CLR/RC wrappers.
- Existing behavior for names and output classes was preserved in R-level wrappers
  through explicit dimname assignment where rewrite returns a matrix object from C++.

## Why RC rewrite was turned off

Repeated A/B benchmarking on the `bmcite` ADT assay (`NormalizeData` methods CLR/RC)
showed rewrite speed was not consistently better for RC:

- Initial rewrite RC was slower than legacy across most runs.
- Multiple optimization passes narrowed overhead, but final median remained slower:
  - RC legacy median remained around ~`0.0645s`
  - RC rewrite median remained around `0.088s` on this benchmark (medians).

Decision taken:
- Keep RC on legacy permanently for now.
- Leave C++ RC path in code for potential future experimentation but do not dispatch to it.

## Performance findings (bmcite ADT repeated benchmark)

Benchmark script:
`/Users/rahulsatija/seurat_rewrite2/scratch/abtest_clr_rc_bmcite_repeat.R`

10-run repeated `bmcite` ADT comparisons were used.

- CLR:
  - Legacy and rewrite reached parity after optimizations.
  - Final observed medians were effectively tied (around ~`0.102s` each), with rewrite slightly better in some runs and very small numerical differences (`max_abs_diff ~1e-15`, `all_equal = TRUE`).
- RC:
  - Legacy remained faster than rewrite in final runs.
  - Final rewrite disabled by policy.

## Additional bmcite LogNorm benchmark (RNA assay)

Benchmark script:
`/Users/rahulsatija/seurat_rewrite2/scratch/abtest_lognorm_bmcite_rna.R`

This benchmark evaluates `NormalizeData(..., normalization.method = "LogNormalize")` on
the `RNA` assay `counts` matrix from `bmcite` using
`GetAssayData(assay = "RNA", layer = "counts")`.

- Legacy `sparse_normalize`:
  - `legacy summary: min = 1.125, q25 = 1.161, med = 1.246, mean = 1.357, q75 = 1.365, max = 2.064`
- Rewrite `sparse_normalize`:
  - `rewrite summary: min = 1.062, q25 = 1.105, med = 1.201, mean = 1.312, q75 = 1.240, max = 2.212`
- speedup (legacy/rewrite):
  - median `1.037`, mean `1.045`, p95 `1.174`
- correctness:
  - exact structure and value match for sparse output (`exact match: TRUE`, `max_abs_diff: 0`)

Raw timing vectors (10 repeats each):
- Legacy: `1.2780, 1.3840, 1.3080, 2.0640, 1.7120, 1.2140, 1.1280, 1.2090, 1.1250, 1.1450`
- Rewrite: `1.0620, 1.2160, 1.2480, 2.2120, 1.6600, 1.2170, 1.0920, 1.1460, 1.0800, 1.1870`

## Net outcome

- Rewrite now provides value primarily for CLR and general sparse log-normalization pathways.
- RC remains legacy for performance and predictability.
- The code now has explicit dispatch points and cleaner internal structure for future A/B toggles while preserving existing outputs.

## Suggested next direction (optional)

1. Keep RC legacy and continue C++ work only if:
   - the input size shifts to where C++ benefits are clear, or
   - a thread-parallel C++ RC kernel is added.
2. Add stronger benchmark harness output (wall-time, memory, and worker count) to the markdown log path for each run.
