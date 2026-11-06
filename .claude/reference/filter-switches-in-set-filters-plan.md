# The plausibility and feasibility switches in `set_filters()`

Date created: 2026-11-05
Date last edited: 2026-11-06

> **IMPLEMENTED 2026-11-06.** See "Outcome" at the end.

(The previous plan in this file, `.poly_clip()` defined once, is implemented and archived as
`.claude/reference/poly-clip-single-definition-plan.md`. On approval, this plan is archived as
`.claude/reference/filter-switches-in-set-filters-plan.md`.)

## Context

### Where the hull and the training ranges filter the local boundaries today

All of it is in `bl_find_local_cf()` (`R/local_cf.R`). For each component pair, after the rotation,
each opposing-class contour segment on the rotated plane passes four filters in this order:

| # | Filter | Thesis term | Code | Switch today |
|---|---|---|---|---|
| 1 | Hull of the training data + target on this plane | plausibility (p.116-117) | built at `:745` (`.build_hull_polygon(..., outlie = 1)`), applied in 2D at `:803` (`.poly_clip()`) | `bl_find_local_cf(hull = TRUE)` |
| 2 | Per-feature constraints | actionability (p.116) | `:818`, `.apply_actionability()` on the back-projected `Bx` | `set_filters()` |
| 3 | Training ranges | feasibility (p.116-117) | `:826`, `get_filter_logical_vector(Bx, train_ranges)` | **none** (only by setting `bl_result$train_ranges <- NULL`) |
| 4 | Model re-score: the point really is on the other side | not a constraint | `:841` | always on |

The nearest surviving vertex is the pair's candidate, and the pair with the smallest distance wins.
After the search, two off-plane points are checked against the training ranges again, not the hull:
the `"fixed"` restore (`.restore_fixed()`, `:1145`) and the `"integer"` rounding
(`.round_integers()`, `:1268`).

Outside the local search: `bl_find_boundary()` (global) clips to the global biplot's hull and only
*flags* `train_ranges` after selection; `bl_find_sparse_cf()` checks neither; `bl_select_target()`
warns when a target is outside `train_ranges`.

### The proposal and my view

Move both switches into `set_filters()`, so one object holds the thesis's three constraints
(actionable, feasible, plausible; p.116), and `set_filters(tgt, hull = FALSE, train_ranges = FALSE)`
gives the thesis's "unconstrained solution" (p.115): the nearest point on the model's local
boundary, even outside the training data. With per-feature constraints given as well, only those
apply. **I think it is better:** the training ranges get a switch they lack today, the result
records which constraints applied, and the API follows the paper's framing. **Costs:**
`set_filters = NULL` must keep both rules on (see "Why `set_filters = NULL` must keep both rules
on" below); a feature
literally named `hull` or `train_ranges` cannot be constrained; about 20 calls change.

**Decisions (user, 2026-11-05):** the names are `hull` and `train_ranges`, documented as
plausibility and feasibility; `bl_find_local_cf()`'s `hull` argument is removed (0.2.0 is untagged
and the argument is one day old).

### Why `set_filters = NULL` must keep both rules on

1. **What NULL does today.** `bl_find_local_cf(bl_results, tgt)` passes no `set_filters`, so
   `set_filters = NULL`. The hull still applies (its own default, `hull = TRUE`) and so do the
   training ranges (always on). A call without filters gives the counterfactual within the hull and
   the ranges, not the unconstrained one.
2. **What changes with the move.** The switches will live inside the `bl_filters` object. A
   `NULL` has no fields, so the code must decide what `NULL` means:
   - **NULL = both on** (planned): every existing call gives exactly what it gives today.
   - **NULL = both off:** every call without `set_filters` silently loses both rules. Its
     counterfactual changes wherever the nearest boundary point lies outside the plane's hull or
     outside a training range.
3. **The calls that pass NULL today** would all be affected:
   - the loan harness's free and external runs (`03_...R:630`, `:702`), whose numbers the walkthrough
     quotes;
   - the vignette's equivalents (`:566`, `:609`) and the README's first example (`:118`);
   - the Pima `run_local()` in `00_setup.R:112` and `run_local_d()` in `04_...R:68`, whenever a run has
     no filters;
   - `05_full-space-gap.R:52`, the A1 experiment. Its reported ratio (median 1.19) is computed with
     the hull on, and its check "every search point lies inside its plane's hull" would fail.
4. **How large the change is.** For the hull alone, the measured effect (`local-convex-hull.md` §4)
   is none on the loan harness, but on Pima 11 of 79 test targets move. One example is the TP
   target: without the hull it goes from pair (3,4) at distance 1.84 to pair (2,3) at 1.50, with
   *Age* raised to 59.6, outside the training data. Dropping the ranges too can move more.
5. **An inconsistency.** `set_filters(tgt)`, a filter object with no constraints, fills in
   `hull = TRUE, train_ranges = TRUE`. If `NULL` meant both off, then "no filter object" and "an
   empty filter object" would give different counterfactuals for the same target.
6. **The thesis.** The method's counterfactual is the constrained one: it "must be selected from a
   local decision boundary that is within the convex hull of the training data" (p.117). Dropping
   the hull is the stated exception, "if an unconstrained solution is acceptable" (p.115). So the
   default stays constrained, and the unconstrained counterfactual is an explicit opt-out:
   `set_filters(tgt, hull = FALSE, train_ranges = FALSE)`.

The roxygen for `bl_find_local_cf(set_filters = )` says this: "`NULL`: no actionability
constraints; the hull and the training ranges still apply (the defaults of `set_filters()`)".

## Design

### `set_filters(bl_target, ..., hull = TRUE, train_ranges = TRUE)` (`R/local_cf.R:414`)

- After `...`, so both must be named exactly; existing calls are unchanged.
- Each must be one `TRUE`/`FALSE`, else `stop(..., call. = FALSE)`: "'hull' must be TRUE or FALSE.
  A feature named 'hull' cannot be constrained through set_filters()." (Same for `train_ranges`.)
- Returns `list(constraints, bl_target, hull, train_ranges)`, class `"bl_filters"`. `constraints`
  stays per-feature only, so `.fixed_vars()`, `.integer_vars()`, `.constraint_labels()`,
  `.cf_compare_table()` and `bl_find_sparse_cf()` are untouched.
- Roxygen: `@param hull` (plausibility: the convex hull of the training data and the target on each
  rotated plane; thesis §4.3.1, p.117) and `@param train_ranges` (feasibility: every feature within
  its training range; p.116-117), a short section "Actionability, feasibility and plausibility",
  and the sentence that both `FALSE` with no per-feature constraint gives the unconstrained
  counterfactual (p.115).
- `print.bl_filters()`: two more lines, `Plausibility (hull)` and `Feasibility (train_ranges)`, each
  `on`/`off`.

### New private helper `.filter_switches(set_filters)` (`R/local_cf.R`, next to `.fixed_vars()`)

Returns `list(hull = , train_ranges = )`. `NULL` gives both `TRUE` (the thesis default, so default
results are unchanged); a `bl_filters` built by hand or saved before this change (no fields) gives
`TRUE` for a missing field.

### `bl_find_local_cf()`

- Remove the `hull` formal and its check; `sw <- .filter_switches(set_filters)`; `hull <- sw$hull`;
  the hull code is otherwise unchanged.
- The ranges filter at `:826` runs only when `sw$train_ranges` (and `is.list(train_ranges)`, as now).
- `.restore_fixed()` and `.round_integers()` gain an argument `train_ranges = bl_result$train_ranges`
  (default keeps their direct tests working); `bl_find_local_cf()` passes `NULL` when the switch is
  off, so their existing `is.list()` checks skip the range test.
- Result: `hull` (now from the switch) and a new logical `train_ranges` field beside it.
- No-solution message: "(hull = FALSE drops this rule)" becomes "(set_filters(..., hull = FALSE)
  drops this rule)", plus the same for the training ranges when that rule is on.
- `print.bl_local_result()`: one line, `Not applied    : ...`, printed **only** when a switch is
  off, so the default output (and every quoted harness number) is unchanged.
- `.warn_target_range()`: "The local search looks for counterfactuals within that range" becomes
  "By default the local search ... (set_filters(..., train_ranges = FALSE) drops this rule)".

### Not changed

`bl_find_boundary()` (global hull and range flag), `bl_find_sparse_cf()`, `bl_shapley()`, filter 4
and the opposing-class contour choice (correctness, not constraints), CLAUDE.md §4's rule against
the global hull in the local search, the version (0.2.0).

## Callers to update (`hull = X` in `bl_find_local_cf()` -> `set_filters(tgt, hull = X)`)

- `development/03_loan_status_Boundary_Logic.R`: Step 10 `set_filters()` (`:546`) gains `hull = TRUE,
  train_ranges = TRUE` with a comment naming plausibility and feasibility; Step 11 (`:562`) drops
  `hull`. The free (`:630`) and external (`:702`) calls stay on the defaults.
- `vignettes/Boundary_Logic_loan_default_workflow.Rmd` (`:485-499`); its `:562` text "without
  actionability constraints gives the geometrically nearest boundary point" becomes "the nearest
  boundary point within the hull and the training ranges; `set_filters(tgt, hull = FALSE,
  train_ranges = FALSE)` gives the unconstrained one". `README.md` (`:348-351`).
- Rotation scripts `documentation/biplot-rotation/17_1`, `17_2`, `17_3`: about 15 calls,
  `hull = FALSE` -> `set_filters = set_filters(<that call's target>, hull = FALSE)`.
- Pima (local): `scripts/paper1-pima/04_method-change-effects.R` (the wrapper `run_local_d()` `:66`,
  which builds `set_filters` with the `hull` switch even when a run has no filters; M6 at `:444`,
  `:463`), `05_full-space-gap.R:52` (`hull = hull` -> `set_filters = set_filters(tgt, hull =
  hull)`); text in `review.qmd` (`:162`).

## Documentation

- `man/`: regenerated (`set_filters.Rd`, `bl_find_local_cf.Rd`).
- `NEWS.md` 0.2.0 "Breaking changes": rewrite the hull bullet for `set_filters(hull = )`; add a bullet
  for the new `train_ranges` switch and the unconstrained counterfactual.
- `2 implementation_summary.txt`: the signature rows of `set_filters()` and `bl_find_local_cf()`, the
  local-search steps (`:381-406`), the filters table (`:432`), the `bl_filters` and `bl_local_result`
  fields.
- `1 Foundation intro documents.txt` (`:253-257`, `:527`): the mechanics wording only.
- `CLAUDE.md`: §2 (local-convex-hull line; one archive line), §4 (`:103`), §6 (filter order line and
  "the rotation scripts pin it").
- `documentation/local-convex-hull.md` (the switch row and §3-§5), `fixed-constraint-restore.md`
  (range check "when `train_ranges` is on"), `biplot-rotation/biplot_rotation_and_pair_selection.md`
  (`:62` and the pinning notes).
- `documentation/methodology-changes/methodology-changes.qmd` (`:225`, `:486`): wording; re-render.
- Walkthrough `development/reference/review_section9_to_15.md`: Parameters tables (`set_filters()`
  gains two rows, `bl_find_local_cf()` loses `hull`), filters 2-3 text, the result-field table, the
  call tree (`:1115`), and the Step 10 `print(flt)` output taken from a real harness run.
- Not changed: archived plans, `future-improvements-review.md` (history), memory files.

## Tests (`tests/testthat/test-local_cf.R`)

1. `set_filters()` stores `hull = TRUE`, `train_ranges = TRUE` by default; `constraints` holds only
   the features.
2. `hull = "increase"`, `NA`, `c(TRUE, FALSE)` (and the same for `train_ranges`) stop with
   "must be TRUE or FALSE".
3. `.filter_switches()`: `NULL` -> both `TRUE`; a hand-built `bl_filters` without the fields -> both
   `TRUE`; `set_filters(tgt, hull = FALSE)` -> `hull` `FALSE`.
4. `print.bl_filters()` shows both switches (`expect_output`).
5. On the iris fixture: `set_filters(tgt, hull = FALSE)` gives `$hull` `FALSE` and `hull_polygon`
   `NULL`; the default gives `$hull` and `$train_ranges` `TRUE`.
6. Unconstrained is never further: with both switches off, the selector distance is <= the default's
   (removing filters only adds candidate vertices); `$train_ranges` is `FALSE`.
7. `.restore_fixed()` with `train_ranges = NULL` uses a restored point outside the ranges (the
   counterpart of the existing test at `:128`).

## Error cases

- A feature named `hull` or `train_ranges`: the validation error above (documented in `@param`).
- An old call `bl_find_local_cf(hull = ...)`: R's "unused argument" error.
- `bl_result$train_ranges` `NULL` with the switch on: skipped, as today.
- Any default-path output differs from `HEAD` (below): stop and report; do not adjust tests.

## Verification

1. **Baseline first:** export `HEAD` to the scratchpad (`git archive HEAD DESCRIPTION NAMESPACE
   LICENSE R inst development documentation/biplot-rotation .Rbuildignore | tar -x`); copy the Pima
   `output/` (`table1.csv`, `numbers.rds`, `method-change-effects.csv`) aside.
2. `devtools::document()`; `devtools::test()`: 473 + new expectations, 0 FAIL, 4 WARN.
3. **Loan harness** against both roots (scratchpad driver per CLAUDE.md §5); diff the logs: the only
   differences are the two new `print(flt)` lines.
4. **Rotation scripts** `17_1`, `17_2`, `17_3` against both roots: identical apart from timings.
5. **Pima:** `run_all.R`, `tests.R` (13/13), `04_method-change-effects.R` (6/6 self-checks); the three
   output files equal the saved copies (`run_log` timestamps aside). Then render `review.qmd` and the
   methodology note (`--to docx`, `QUARTO_R` set).
6. One new-behaviour check printed from a driver: the harness target with both switches off, its
   pair and distance against the default's.
7. `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`: 0 / 0 / 0.
8. Update the walkthrough's quoted Step 10 output from step 3's log; archive this plan; update
   `progress.md`. Stop and report; no commit or push without an instruction.

## Outcome (2026-11-06)

Implemented as planned. Additions during implementation:
- `.check_switch()` (private) validates the two switches; `.filter_switches()` reads them with
  `[[` (exact names).
- The harness's Step 15 title and the vignette's "Alternative" section now say "without
  actionability constraints" (they apply the hull and the ranges); a harness comment shows the fully
  unconstrained call.
- README: the "Actionability constraints" section and the function table name the switches.
- `local-convex-hull.md` §6: where the rules run, why NULL keeps both on, and the measured effect of
  the switches on the harness targets.
- The methodology note's history table: the 5 Nov row now cites `c13fe44`; a 6 Nov row was added.

Results:
- `devtools::test()`: 514 PASS / 0 FAIL / 4 WARN / 0 SKIP (473 + 41 new expectations).
- Loan harness against `HEAD`: identical apart from timestamps, the two new `print(flt)` lines and
  the reworded target-range warning.
- Rotation scripts: `17_1`, `17_2` identical; `17_3` timings only.
- Pima: `tests.R` 13/13; `04` 6/6 self-checks; `05` both checks. `table1.csv`, `full_space_gap.csv`
  and `.rds` identical; `numbers.rds` identical apart from `meta` and the no-solution message of TP
  runs 2-4 (the new wording); `method-change-effects.csv` differs only in four `meta` rows. No figure
  or review text quotes the message. `review.html` and `methodology-changes.docx` re-rendered.
- New behaviour (harness, no actionability constraints): the Step 10 target is unchanged under every
  setting (pair (3,4), sqrt(d2_M) 0.5554). The Step 16 external applicant (credit score 420, below the
  training range): default pair (1,2), 1.7026; `train_ranges = FALSE` pair (1,4), 0.8521; both off
  pair (1,4), 0.4750.
- `R CMD check`: 0 / 0 / 0 (run twice; the second after the last vignette edit).
- Not committed.
