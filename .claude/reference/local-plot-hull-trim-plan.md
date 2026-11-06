# The local biplot shows the hull, as the global biplot does

Date created: 2026-11-06
Date last edited: 2026-11-06

> **IMPLEMENTED 2026-11-06.** See "Outcome" at the end.

(The previous plans in this file are implemented and archived as
`.claude/reference/poly-clip-single-definition-plan.md` and
`.claude/reference/filter-switches-in-set-filters-plan.md`. On approval, this plan is archived as
`.claude/reference/local-plot-hull-trim-plan.md`.)

## Context

With the hull switch on (`set_filters(hull = TRUE)`, the default), the local search keeps only the
boundary points inside the convex hull of the training data and the target on each rotated plane.
The local biplot (`plot.bl_local_result()`, and `plot.bl_sparse_result()` through it) does not show
this: it colours the whole prediction grid, so a reader cannot see which boundary lines the search
could choose from. The user asked for the Pima local figures in `review.html` to show the hull
filter, looking like Figure 5(a), the global biplot; and for Figure 7(a), now a text panel ("no
counterfactual"), to be a real plot whose caption says why no counterfactual exists.

**What Figure 5(a) does.** `bl_build_grid(calc_hull = TRUE)` keeps only the grid points inside the
global hull (`R/biplot_grid.R:239`, `!is.na(sp::over(G_df, polygon))`), so the colours stop at the
hull and the area outside is white. The contour lines are drawn in full. No outline, no legend.

**Why 7(a) has no counterfactual.** The TP target's run 4 (Pregnancies fixed, Age cannot decrease)
finds none with the hull. Without the hull it finds one on pair (1,2): Pregnancies 7 (held), Age 51
→ 59.5, √d²_M 1.99 (`method-change-effects.csv`, M6). Inside the hull every TP boundary point lowers
Age (`local-convex-hull.md` §4), which run 4 forbids.

**Decisions (user, 2026-11-06):**
- the local plot looks like Figure 5(a);
- no legend: the captions explain it;
- Figure 7(a) shows run 4 searched without the hull, on its plane, with the colours trimmed to that
  plane's hull: the cross marks the point the hull excludes, and the caption explains.

## Design

### 1. `bl_find_local_cf()` stores the winning plane's hull also when the hull filter is off (`R/local_cf.R`)

After the pair loop, when a solution exists and `hull` is `FALSE`, build the hull for the winning pair
the same way the loop does: `.build_hull_polygon(rbind(Z_train_rot, Z_target), min_val, max_val,
outlie = 1)` from `best_result`, and store its vertices in `hull_polygon`. One extra off-screen hull
per call, only with the switch off; the search and every number are unchanged. `hull` still says
whether the rule was applied. Reason: the plot must not compute geometry (CLAUDE.md §4,
"never entangle plotting code with computation code"), and 7(a) needs the hull of a hull-off search.

### 2. New private helper `.grid_in_hull(Z, hull_coords)` (`R/hull_utils.R`, next to `.points_in_polygon()`)

Logical vector, `TRUE` inside the hull: rebuilds the `SpatialPolygons` object from the stored vertices
(as `.build_hull_polygon()` builds it) and calls the existing `.points_in_polygon()` (`sp::over`), the
global grid's rule. `NULL` vertices give all `TRUE`.

### 3. `plot.bl_local_result()` (`R/local_cf.R:1756`)

- New argument `trim_to_hull = NULL`, after `no_grid`. `NULL` means `isTRUE(x$hull)` (trim when the
  search applied the hull); `TRUE` trims also for a hull-off search (7(a)); `FALSE` never trims.
- When trimming and `x$hull_polygon` is not `NULL`, Step 2 draws only the grid points inside the hull:
  `keep <- .grid_in_hull(x_grid, x$hull_polygon)`. Everything else (contours, axes, target, cross,
  arrow) is drawn as now.
- Display only; an object built before `hull_polygon` existed gets the full grid.
- Roxygen: `@param trim_to_hull`; one sentence in the description.
- `plot.bl_sparse_result()` passes `...` through (figures 9(a), 9(b)); one roxygen line.

### 4. Pima scripts (`scripts/paper1-pima/`, local)

- `00_setup.R`: `run_local()` gains `hull = TRUE`, built like `04`'s `run_local_d()`
  (`flt <- if (length(filters) || !hull) do.call(set_filters, c(list(tgt), filters, list(hull = hull)))`).
  The default leaves every run as now.
- `03_local.R`: for each target whose run 4 has no counterfactual, run 4 again with `hull = FALSE`;
  print it in the run log; store `extract_numbers()` of it plus `cross_in_hull` and `search_in_hull`
  (`.grid_in_hull()` of `B_z` and of `fixed_check$search_point$B_z` against its `hull_polygon`) as
  `numbers$run4_no_hull[[tg]]`; store the object in `local.rds` as `run4_no_hull`. Table 1, `runs` and
  every existing number are untouched.
- `fig-07.R`: when run 4 has no counterfactual and `run4_no_hull` has one, plot it with
  `trim_to_hull = TRUE`, title "TP target, Pregnancies fixed, Age increase: no counterfactual inside
  the hull"; the text panel stays only as a fallback (neither has a solution). Comment block updated.
- `review.qmd`:
  - **7(a) caption, from `numbers$run4_no_hull$TP`:** no counterfactual. Searched without the hull,
    run 4 finds one on pair (1,2), the cross: it holds Pregnancies and raises Age from 51 to 59.5. It
    lies outside the convex hull of the training data on its plane (the white area), so the hull
    removes it. Inside the hull every boundary point lowers Age.
  - **If `cross_in_hull` is `TRUE`:** the caption says instead that the search point lies outside the
    hull, and that restoring Pregnancies moves its projection.
  - **Captions of 5(b), 7(b), 7(c), 9(a), 9(b):** one shared sentence: colours are drawn only inside
    the convex hull of the training data and the target on this plane, as in Figure 5(a); boundary
    lines in the white area lie outside it, so the search could not choose them.

### Not changed

The search result's numbers, `plot_biplotEZ()`, the global grid, Table 1, the methodology note and its
measurements (04 is not re-run: it reads no hull polygon), the harness code.

## Documents

- `NEWS.md` 0.2.0: the `trim_to_hull` argument; `hull_polygon` stored also with the hull off (edit the
  hull bullet's "`hull_polygon` (the winning plane's hull)").
- `R/local_cf.R` roxygen `@return hull_polygon`: stored for the winning plane in either case.
- `2 implementation_summary.txt`: `hull_polygon` (the "NULL when the hull switch is off" sentence) and
  the `plot.bl_local_result()` paragraph.
- `development/reference/review_section9_to_15.md`: the `hull_polygon` result row; Step 12's layer 2
  and a Parameters row for `trim_to_hull`; Step 14's pass-through sentence.
- `vignettes/Boundary_Logic_loan_default_workflow.Rmd` Step 14: one sentence.
- `documentation/local-convex-hull.md` §3: the stored hull and the plot.
- `progress.md`: figures 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) changed (for the Paper 1 plan).

## Tests

- `test-hull_utils.R`: `.grid_in_hull()` on the unit square (inside, outside, `NULL` vertices).
- `test-local_cf.R`:
  - the existing `expect_null(no_hull$hull_polygon)` becomes: `hull` `FALSE`, `hull_polygon` stored;
  - on the iris fixture (`pdf(NULL)`), `plot()` runs with `trim_to_hull` `NULL`, `TRUE` and `FALSE`,
    for a default and a hull-off search, and with `hull_polygon` removed (an older object);
  - `.grid_in_hull(loc$Zgrid, loc$hull_polygon)` keeps between 0 and `m^2` points (exclusive).

## Error cases

- Malformed `hull_polygon` (fewer than 3 vertices): `sp::Polygon()` errors; not reachable from
  `bl_find_local_cf()`. No extra handling.
- A restored, rounded or sparse cross can fall in the white area (off-plane points; the hull applies on
  the plane only, `local-convex-hull.md` §3). Correct; mentioned in a caption only where it happens.
- `fig-07.R` with neither result solved: the text panel, as now.

## Verification

1. `devtools::document()`; `devtools::test()`: 514 + new, 0 FAIL, 4 WARN.
2. Rotation scripts `17_1`-`17_3` (the hull-off callers) and the loan harness against the current
   logs: identical apart from timestamps and timings.
3. Pima: `run_all.R`; `table1.csv` identical; `numbers.rds` identical apart from `meta` and the new
   `run4_no_hull`; `tests.R` 13/13. View fig-05a next to 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) by eye
   (stepped colour edge, lines on white outside, 7(a)'s cross position); then `quarto render
   review.qmd` and read the 7(a) caption.
4. `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`: 0 / 0 / 0.
5. Archive this plan; update `progress.md`. Stop and report; no commit or push without an
   instruction.

## Outcome (2026-11-06)

Implemented as planned. Small additions: `.grid_in_hull()` returns an unnamed vector (`sp::over()`
names its result); `fig-07.R`'s 7(a) title reads "TP target, Pregnancies fixed, Age increase: / no
counterfactual inside the hull (pair (1,2), hull off)"; `03_local.R` also runs the TN target's run 4
without the hull (it has none either: the fixed restore fails), stored but not plotted.

- `devtools::test()`: 527 PASS / 0 FAIL / 4 WARN / 0 SKIP (514 + 13).
- Loan harness and rotation scripts `17_1`-`17_3`: identical to the previous logs apart from
  timestamps and `17_3`'s timings.
- Pima: `table1.csv` identical; `numbers.rds` identical apart from `meta` and the new `run4_no_hull`;
  `tests.R` 13/13. TP run 4 without the hull: pair (1,2), Pregnancies 7, Age 59.52, √d²_M 1.99; its
  cross and its search point both lie outside the plane's hull (`cross_in_hull`, `search_in_hull`
  FALSE). Figures 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) checked by eye against 5(a); `review.html`
  re-rendered with the new captions.
- `R CMD check`: 0 / 0 / 0.
- Not committed.
