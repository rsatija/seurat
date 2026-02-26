# FindVariableFeatures Rewrite Summary

## Scope
This summarizes changes related to `FindVariableFeatures` and its internal selectors (`VST`, `DISP`, `MVP`) across:

- `R/preprocessing5.R`
- `R/preprocessing.R`

It includes implementation rewrites, dispatch changes, and documentation updates focused on speed, parity, and maintainability.

## 1) Dispatch architecture introduced for A/B testing

### Added implementation switches
- `FindVariableFeatures.V3Matrix`
- `FindVariableFeatures.Assay`
- `DISP`
- `MVP`
- `FindVariableFeatures.StdAssay` (earlier in session)

Each now selects legacy/rewrite behavior from options (e.g. `Seurat.FindVariableFeatures.Assay.impl`) using:
- `legacy` path: existing behavior
- `rewrite` path: optimized implementation

#### Rationale
- Enables deterministic A/B benchmarking without API changes.
- Lets teams toggle behavior per process without code edits.
- Creates a controlled migration path for speed gains while preserving exact outputs.

#### Performance expectation
- Removes repeated switching overhead by enabling per-run opt-in.
- Speeds iteration during development/validation because rewrite can be benchmarked quickly against legacy by toggling options.

---

## 2) Shared top-k selection primitive

### New helper
- `.TopKIndices(x, k, decreasing = TRUE)` in `R/preprocessing5.R`

Used by:
- `DISP_legacy`
- `DISP_rewrite`
- `VST.IterableMatrix`
- `VST.dgCMatrix`
- `VST.matrix`/`.VST` path via internal calls
- `FindVariableFeatures.Assay_rewrite` (vst branch)

#### Rationale
- Duplicate top-k code previously repeated in multiple selectors and methods.
- Needed both speed and deterministic behavior when ties occur at selection boundary.

#### Behavior contract
- Uses fast partial selection (`sort.int(..., partial = ...)`) when safe.
- Detects ties at the boundary and falls back to full deterministic ordering to preserve legacy-like selection semantics.
- Handles zero/NA/inf-safe boundary conditions and empty requests consistently.

#### Performance expectation
- Avoids full vector ordering for common “small `nselect`” cases.
- Reduces memory churn from repeated `order()` allocations for large feature spaces.
- Keeps tie-case exactness only when necessary.

---

## 3) DISP rewrite and cleanup

### Changes
- `DISP` dispatch implemented with legacy + rewrite branches.
- `DISP_legacy` now uses `.TopKIndices` for feature ranking.
- `DISP_rewrite` added as an optimized path with same output schema.

#### Rationale
- Centralize and unify variance-dispersion ranking behavior.
- Eliminate duplicated ranking blocks while preserving outputs.

#### Performance expectation
- Faster selection of top dispersions when only top features are needed.
- Reduced transient vector sorting overhead on large matrices.

---

## 4) MVP rewrite with reduced sorting surface

### Changes
- Added `MVP` implementation switch (legacy vs rewrite).
- `MVP_legacy` retains historical pipeline.
- `MVP_rewrite` now:
  - filters on `mean.cutoff`/`dispersion.cutoff` first,
  - ranks only selected rows by dispersion,
  - writes ranks on selected indices only.

#### Rationale
- Previous behavior sorted full `hvf.info` prior to filtering for all cases.
- Reordering to “filter then rank” avoids unnecessary global work and makes logic clearer.

#### Performance expectation
- Saves sorting work when the filtered set is much smaller than total features.
- Cuts intermediate indexing overhead in large datasets with strong cutoffs.

---

## 5) VST path cleanup and parity alignment

### Changes
- `VST.IterableMatrix`:
  - replaced `head(order(...), n)` selection with `.TopKIndices`-based branch.
- `VST.dgCMatrix`:
  - same top-k branch refactor.
- `.VST`:
  - removed redundant top-k branch logic and now uses a single `.TopKIndices` fast path plus clear full-sort fallback for broad selections.
- Legacy/rewrite DISP usage in `.VST` now aligned around shared helper.

#### Rationale
- Ensure identical intent across multiple VST code paths.
- Remove redundancy that increased cognitive load and potential divergence.
- Preserve behavior when all or nearly all features are requested.

#### Performance expectation
- Faster top feature extraction for common default requests (`nselect` far smaller than total features).
- Predictable performance by avoiding duplicate partial sort checks.

---

## 6) FindVariableFeatures.Assay rewrite modernization

### Changes
- `FindVariableFeatures.Assay_rewrite` now avoids full sort for VST top-feature extraction:
  - uses `.TopKIndices(x = hvf.info$vst.variance.standardized, k = nfeatures)`
- Added explicit legacy/rewrite wrapper functions and preserved legacy wrapper contract.

#### Rationale
- Reduce work in the most frequently used user-facing path.
- Keep downstream behavior (`VariableFeatures`, `<method>.variable`) identical while optimizing ranking.

#### Performance expectation
- Faster selection for standard `nfeatures=2000` defaults on large feature sets.
- No additional object materialization for ordering full feature vectors.

---

## 7) Documentation and readability updates

### Changes
- Added concise function-level comments for each rewrite-relevant function in:
  - `R/preprocessing5.R`
  - `R/preprocessing.R`
- Added inline comments at key decision points:
  - sparse vs dense branches,
  - top-k path vs full sort,
  - tie-detection fallback,
  - legacy vs rewrite behavior boundaries.

#### Rationale
- Meet rewrite-goal requirement: clean, understandable code for long-term maintainability.
- Make behavior decisions visible to new contributors.

#### Performance expectation
- No runtime change; reduces maintenance and onboarding cost while lowering regression risk.

---

## Net effect on result reproducibility and speed

- Rewrites were built to preserve output contracts of legacy paths:
  - same returned columns, same flag semantics (`variable`, `rank`, `VariableFeatures`),
  - same downstream object updates.
- Performance gains come from reduced sort work, partial selection, and reduced redundant processing.
- Tie-aware fallback ensures correctness and stable semantics remain aligned with existing expectations.

## Notes on benchmark status

- A/B benchmark scripts were run in prior sessions for this rewrite stream and earlier iterations.
- The code is now structured for immediate benchmark reuse via the implementation options above.
- For a fresh, auditable benchmark cycle, run a dedicated benchmark script that toggles:
  - `options(Seurat.FindVariableFeatures.*.impl = "legacy")`
  - `options(Seurat.FindVariableFeatures.*.impl = "rewrite")`
  - and records wall time + exactness checks on the same matrix/method combinations.
