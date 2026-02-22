# Purpose: Seurat Rewrite

We are incrementally rewriting Seurat for:
- Speed: reduce runtime and allocations in core pipelines.
- Robustness: improve input validation and edge-case handling.
- Maintainability: clearer internal boundaries, consistent conventions, and strong comments.
- Backwards compatibility: exported APIs and outputs must remain identical to legacy Seurat.

Approach:
- Freeze current behavior with golden datasets + pipelines.
- Add legacy/rewrite dual implementations with dispatch via reserved arg `.__seurat_engine` passed through `...`.
- Maintain identical outputs and object structure; differences must be detectable and explained.
- Every rewrite produces review artifacts and passes golden verification and benchmark checks.

# Rules (MUST follow)

## Broader purpose (always optimize for this)
We are rebuilding Seurat step-by-step to be faster, more robust, and better commented, while preserving backwards compatibility.
The key requirement is that users see identical outputs and object structure for the same inputs and parameters.

## Non-negotiables
1) Backwards compatibility:
   - Do NOT change exported function names, signatures, defaults, or return types/structure.
   - Do NOT change object slots/fields, dimnames, ordering, or factor levels unless explicitly approved.
2) Identical outputs:
   - All rewritten behavior must match legacy outputs on golden datasets/pipelines (strict or explicitly defined tolerances).
3) No global options:
   - Engine selection is ONLY via reserved `...` argument `.__seurat_engine` (default "legacy").
4) Don’t regenerate goldens:
   - Never modify golden expected outputs unless explicitly instructed with a regeneration ticket.

## Dispatch pattern (required)
For any routed function <name>:
- Implement:
  - <name>_legacy  (legacy code moved verbatim)
  - <name>_rewrite (new code, well-commented)
  - <name>         (wrapper/dispatcher)
- The wrapper must:
  1) extract `.__seurat_engine` from `...` (default "legacy")
  2) remove it from dots
  3) dispatch to <name>_legacy or <name>_rewrite
- Any function calling a routed function must forward `.__seurat_engine`.

## Reserved args policy (strict)
- Reserved arg: `.__seurat_engine` only. Allowed values: "legacy", "rewrite".
- Always strip reserved args before calling external packages (do not leak `.__seurat_engine` to deps).

## Commenting standard (required for *_rewrite)
Each <name>_rewrite must include a header block:
- Purpose
- Inputs/Outputs
- Determinism/RNG notes
- Compatibility notes (what must match legacy)
Plus:
- Section comments for major phases
- Any line that exists to match legacy must include `# COMPAT: ...`

## Review artifacts required
Any rewrite must produce/update:
- diffs/<component>.patch
- review/<component>.md
The review doc must include:
- short summary of changes
- explicit compatibility notes
- side-by-side or unified diff
- which golden tests were run/passed

## Tests required
All changes must add/maintain tests:
- golden tests for outputs
- unit tests for edge cases
- no new warnings/errors unless explicitly approved
