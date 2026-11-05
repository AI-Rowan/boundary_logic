# Plan: from open issues to future improvements, and one design-decisions document

> **IMPLEMENTED (steps 0-5), 2026-11-05.** Archived from `.claude/plans/cryptic-snuggling-bee.md`.
> The hull on each rotated plane (B10), A1, A3, A8, A9, B3, B9 and B24 are in the working tree, not
> committed. Results: `devtools::test()` 465 PASS / 0 FAIL; Pima `tests.R` 13/13; `devtools::check()`
> 0/0/0. Phases 4 and 5 (the design-decisions document, retiring `known-issues.md`) wait for the
> supervisor's decisions. The decisions themselves are in
> `documentation/methodology-changes/future-improvements-review.md`.

Date created: 2026-11-05
Date last edited: 2026-11-05

## Context

Every open item from Paper 1's preparation (the 9 known issues, about 24 other items, and the
resolved ones) has been through the reply form
`documentation/methodology-changes/future-improvements-review.md` (three rounds, 2026-11-05).

**Round 4 (2026-11-05), two changes:**
- **No test plan.** The package only has to work on your computer, so the tests of computed values
  (B21) and the CVA identity test (B5) are not done now. B2's paper sentence notes the risk.
- **The convex hull on the rotated plane (B10) is done now, not as future work.** You decided the
  methodology note's §3.6 yourself: the local search returns to the thesis rule (p.117). The
  supervisor has not seen the methodology note yet (corrected 2026-11-05); the note presents §3.6
  as restored.

Your design choices for B10 (AskUserQuestion, 2026-11-05):
- the hull is built from the training data plus the target, as Appendix A;
- it is switchable, `hull = TRUE` by default;
- the rotation scripts pin `hull = FALSE`;
- there is no enclosure rule.

**Consequence:** B10 changes Paper 1's Phase 3 numbers (Table 1, figures 5–9, the review page,
the methodology note). Phases 1 and 2 are unchanged. No commit or push without your instruction.

**Measured during step 1 (2026-11-05), and your decisions on it:**
- `hull = FALSE` reproduces the earlier Pima results exactly. With `hull = TRUE`:
  - the TN target is unchanged;
  - the TP target moves to pair (3,4) (1.50 to 1.84), and runs 3, 4 and 5 have no counterfactual,
    since inside the hull every TP counterfactual lowers *Age*;
  - over the 79 test rows, 11 move and none is lost.
- Your decisions: keep the hull, and keep figure 7(a) on TP. When run 4 has no counterfactual,
  `fig-07.R` draws a panel that says so, instead of skipping it.
- Add **run 7** (*Pregnancies* fixed, *Glucose* `"decrease"`) to `runs_spec` in `00_setup.R` and
  to `table1_cols` in `03_local.R`.

**Additions made during the implementation (2026-11-05):**
- `1 Foundation intro documents.txt` stated that the local search applies no hull "intentionally". It
  takes precedence on methodological questions, so its three passages now say the global polygon is
  not used and that a hull is rebuilt per rotated plane.
- The Pima `issues.csv` (local): issues 5, 11 and 12 rewritten for the hull. Issue 11's "older
  patient" no longer holds, since the TP counterfactual now lowers *Age*.
- The methodology note: every TP-dependent passage revised (§3.4, §3.6, §3.7, §3.9, §3.10, §3.12,
  Table 4, the interpretation points), and its claims check updated.
- `devtools::check()` runs on a copy of the package files, because the note's `.docx` was open in Word
  and R CMD build could not copy it.

---

## The final list

### 1. Already implemented

These were fixed before or during the Task D list (`documentation/known-issues.md` §2). They
become design decisions in the final document.

| # | Issue | Decision taken |
|---|---|---|
| C1 | The nearest boundary vertex was chosen with a random tie-break, so results drifted between calls | The exact nearest vertex (`ties.method = "first"`), 2026-10-30 |
| C2 | The cutoff was not passed to the grid | It is passed |
| C3 | `rotate_deg` could turn the local biplots | Removed from the local plots (an error if passed) |
| C4 | Models from `xgboost::xgboost()` were not accepted | Accepted |
| C5 | An empty CVA outcome group broke the projection | Empty groups are dropped, with a message |
| C6 | The PCA pair search used the Mahalanobis distance, against the thesis | `distance = "auto"`: Euclidean under PCA, Mahalanobis under CVA, 2026-11-05 |
| C7 | `used_in_sparse` was wrong for fixed variables | Fixed |
| C8 | Stale prose and line references in the notes | The documentation sweep, 2026-11-04 |

### 2. Implemented now

| # | Issue | Decision | Moves Paper 1's numbers? | In the paper? |
|---|---|---|---|---|
| **B10** | The local search does not keep the counterfactual within the convex hull of the training data (thesis p.117) | A hull rebuilt on each rotated plane from the training data and the target; `hull = TRUE` by default | **Yes: Phase 3** | Yes: §5, as part of the method |
| A1 | Paper 1's distances are upper bounds of unknown tightness | A script-only experiment on the 79 test rows (local search only), against a constrained full-space search, with the hull checks | No; it adds a result | Yes: §8 |
| A3 | `plot.bl_projection()` uses a different colour scheme | One shared colour map, on the thesis scheme | No | No |
| A8 | The README constraint example finds no counterfactual | Paper 1's run 4 constraints, with the printed result and one sentence on `blocking_constraint` | No | No |
| A9 | The loan harness title says "prior defaulters" | "Loan default (no prior default) -- XGB, CVA biplot" | No | No |
| B3 | "Eigenvector pair" suggests the plane is the eigenplane | The agreed paper sentence; the help pages, messages and walkthrough reworded | No | Yes: §5 |
| B9 | A target at the centroid makes the Procrustes matrix zero | A warning only | No | No |
| B24 | Housekeeping before the release | Version 0.2.0, NEWS, final thesis page citations, pkgdown; merge and tag on your instruction | No | No |

### 3. Future improvements

**In the paper** (§9's four paragraphs from G1, unless another section is named):

| # | Issue | Paper placement |
|---|---|---|
| A2 | Only the winning pair's plane and point are kept | §5, sentence (ii) (also covers B12) |
| A4 | The fixed tolerance (0.5 model units) cannot be changed | Constraints; sentence (c) |
| A5 | A failed fixed restore gives no counterfactual | Constraints |
| A6 | The rotation and pair search include the fixed variables | Constraints |
| A7 | Integer rounding cannot let the free variables compensate | Constraints (its numbers re-checked after the re-run) |
| B1 | The local plane is not unique; B4 (frame convention) goes with it | The local plane |
| B2 | Cross-platform reproducibility is untested | The local plane, revised sentence below |
| B6 | Feasibility is a filter on contour vertices, not a constraint in the search | Constraints |
| B7 | The reported counterfactual lies off the plane | §5: Table 1's two distances |
| B12 | The fixed restore is not retried on the next-best pair | Through A2's sentence |
| B13 | Integer rounding at other step sizes | Constraints (a phrase) |
| B14 | The per-variable importance overstates the full Mahalanobis distance by 27% | Attribution |
| B15 | No aggregated Shapley values | Attribution |
| B16 | CVA needs a nonsingular W (p > n) | Scope |
| B17 | Binary classifiers only | Scope |
| B18 | Only affine transformations on the axes | Scope |

**Revised B2 sentence:** "The package was developed and run on one platform (R 4.6.0 on Windows).
The Procrustes rotation determines the local plane only up to its second direction, which the
linear algebra routine settles, so the results are reproducible on that platform, but another
platform may select a different plane, and hence a different counterfactual; the minimal rotation
proposed as future work would make the plane unique." The supplement includes `sessionInfo()`.

**In the design document only:**

| # | Issue | Note |
|---|---|---|
| B5 | No unit test pins the CVA identity | Round 4: future, with B21 |
| B8 | No `pairs =` argument | |
| B11 | No feasibility flag in the global search | 0 of 313 and 0 of 79 on the seven-variable model |
| B21 | Few tests of computed values | Round 4: future, with B23 (only if the package goes beyond your computer) |
| B22 | No method vignette | |
| B23 | Not prepared for CRAN | |

### 4. No action

B19 (non-interactive arrow layer) and B20 (label direction per variable).

### 5. Limitations the paper states (not future work)

From the rotation note §4.9:
- a 2D biplot is a slice of the decision surface;
- rotating to fit one target worsens the other points' representation;
- the search finds the nearest boundary within the planes it draws (A1 measures this);
- the CVA class means leave the plane;
- "optimal" refers only to the smallest distance.

---

## Implementation (after approval)

### Step 0: update the reply form

In `future-improvements-review.md`, add a Round 4 note and update:
- B10 becomes "now", with your four design choices, and its status Resolved;
- B21 and B5 become future, design document only;
- remove the tests plan from the Round 3 table;
- the revised B2 sentence;
- G1 (B10 moves from §8 to §5);
- G2 gains B21.

### Step 1: the hull on each rotated plane (B10)

**Read first:**
- `R/local_cf.R:626-910` (the pair loop and the no-solution messages);
- `R/hull_utils.R` (`.build_hull_polygon()`, `.poly_clip()`);
- `R/boundary.R:250-300` (how the global search clips to its hull);
- `R/biplot_grid.R:176` (the global grid's `outlie = 1`);
- the methodology note §3.6.

**Design** (in `bl_find_local_cf()`):
- **The argument.** A new argument `hull = TRUE`, after the existing arguments, so positional calls
  are unaffected. It is validated as a single TRUE or FALSE (`call. = FALSE`).
- **The hull.** Per pair, after `Z_train_rot`, `Z_target`, `min_val` and `max_val` (`:700-708`):
  `hull_poly <- .build_hull_polygon(rbind(Z_train_rot, Z_target), min_val, max_val, outlie = 1)`
  when `hull` is TRUE.
  - The full hull (fraction 1) matches the global grid.
  - The helper already restores the caller's graphics device.
- **The filter.** In the segment loop, right after `Mi` is built (`:757`), and before the
  actionability, `train_ranges` and model filters: `Mi <- .poly_clip(Mi, hull_poly)`, then `next`
  if it is empty.
  - This is the same vertex filter the global search uses.
  - `.poly_clip()` in `hull_utils.R` is the authoritative copy.
  - The global grid's Z-space polygon is still never used, as CLAUDE.md §4 requires.
- **The result.** Two new fields: `hull` (the flag) and `hull_polygon` (the winning pair's polygon
  coordinates). No field is removed.
- **The messages.** When `hull` is TRUE, the no-solution messages (`:893-903`) name the hull, and
  say that `hull = FALSE` drops the rule (text only).
- **Off the plane.** The restore, the rounding and the sparse counterfactual work on off-plane
  points. The hull applies to the search point on the plane only, as `train_ranges` already does
  to the in-plane candidates. The B7 sentence covers this.

**Documentation:**
- the roxygen: `@param hull`, and the `@details` text at `:457-458`, rewritten to say the global
  polygon is not used and that a hull is rebuilt per plane (thesis p.117, Appendix A pp.207, 209);
- a new technical note, `documentation/local-convex-hull.md`: why (the thesis text, Appendix A,
  and the branch's test-data hull), what the package does, the switch, and the measured effect
  (M6, below);
- `2 implementation_summary.txt`: the signature and the two new fields;
- NEWS.md, under breaking changes;
- CLAUDE.md:
  - the §4 rule, reworded: the *global* polygon is never applied; a hull rebuilt per rotated plane
    is applied;
  - §6 ("Only `train_ranges` and `set_filters` constrain the local search");
  - the line on running the Pima example (adds `05`);
- **The harness example shows the argument (your comment):** the Step 11 call in
  `development/03_loan_status_Boundary_Logic.R:559` gains `hull = TRUE` explicitly, with a comment
  (the thesis plausibility rule; `hull = FALSE` drops it). The free (`:626`) and external (`:698`)
  calls keep the default.
- `development/reference/review_section9_to_15.md`:
  - the Step 11 call block (`:339`) gains `hull = TRUE`;
  - its Parameters table (`:348`) gains a row: `hull` | `TRUE` | `TRUE` | keeps only boundary points
    inside the convex hull of the training data and the target on each rotated plane (thesis
    p.117); `FALSE` drops the rule;
  - the stage-by-stage logic gains the hull filter, before the actionability filter;
  - the "convex hull is NOT applied" passage (`:1178`) is rewritten;
  - the numbers and console quotes are updated after the re-run.
- The vignette `vignettes/Boundary_Logic_loan_default_workflow.Rmd` mirrors the harness: its Step 13
  call (`:493`) gains `hull = TRUE`, and its text (`:483`) names the hull filter.
- `README.md:326` ("Local search does not apply the convex hull") is rewritten to describe the hull
  on each rotated plane, and `hull = FALSE`.

**Error cases:**
- the hull removes every candidate on every pair: no solution, and the message names the hull;
- a target outside the training cloud: it is part of its own hull;
- existing tests that change because they depend on the search point: pin `hull = FALSE` only where
  the test checks something else (the restore, the rounding). No new tests.

### Step 2: small package fixes (A3, A8, A9, B3, B9)

- **A3:**
  - a private constant `.confusion_cols <- c(TP = "red", TN = "blue", FP = "purple", FN = "orange")`
    in `R/utils.R`;
  - used at `R/projection.R:320` (and its roxygen at `:376`), `R/boundary_plot.R:170`,
    `R/project_points.R:129` and `R/local_cf.R:1682,1705`;
  - the expectation at `tests/testthat/test-projection.R:72` is updated.
- **A8:** at `README.md:188`, with the hull on:
  - the run 4 constraints on the same patient. If no counterfactual is found, drop *Glucose*, then
    try the next predicted-diabetic patient;
  - run through a scratchpad driver;
  - paste the printed table, plus one sentence on `blocking_constraint`.
- **A9:** the title at `development/03_loan_status_Boundary_Logic.R:251`, and the walkthrough quote
  at `development/reference/review_section4_to_6.md:358`.
- **B3:** "eigenvector pair" becomes "component pair":
  - in `bl_find_local_cf()`'s roxygen (`:441-572`) and its messages (`:678`, `:895`, `:900`);
  - in `plot.bl_local_result()`'s roxygen (`:1520`), with "optimal rotation" becoming "Procrustes
    rotation" (`:1628`);
  - one `@details` sentence with the agreed meaning;
  - the same wording in `review_section9_to_15.md`.
- **B9:**
  - after `x_target_st` (`R/local_cf.R:669`), warn once (`call. = FALSE`) when
    `sqrt(sum(x_target_st^2)) < sqrt(.Machine$double.eps)`;
  - first confirm today's behaviour in a scratchpad script (expected: every pair gives plane
    (1,2));
  - `.bl_rotate()` is untouched.
- **Also:** NEWS.md, the implementation summary, the rotation note §4.8, and `known-issues.md`
  (1.3, 1.8 and 1.9 move to its §2 table).

### Step 3: check, re-run and update the numbers

1. Run `devtools::document()` and `devtools::test()`, and report PASS/FAIL.
2. **Loan harness:** re-run it by the CLAUDE.md driver, then update the numbers and console quotes
   in `review_section9_to_15.md`.
3. **Pima:**
   - run `run_all.R`, then `tests.R`, then render `review.qmd`;
   - record which cells of Table 1 and which figures change;
   - update the review page's text where it no longer holds.
4. **Methodology effects:**
   - add **M6** to `04_method-change-effects.R`: `hull = TRUE` against `FALSE` for both targets and
     the runs, as the other switches;
   - run it, and update the methodology note's §3.6: thesis, package (now restored), reason,
     effect, and "decided by Adriaan Rowan, 2026-11-05";
   - update its History and Table 5, and render to `.docx`;
   - the `claims` chunk must pass. Phase 4 still does the full restructure.
5. **Rotation scripts:** add `hull = FALSE` to the `bl_find_local_cf()` calls in `17_1`, `17_2`
   and `17_3`, beside their `distance = "mahalanobis"` pin. Add one line in the rotation note's
   Appendix B. Re-run `17_3` to confirm its output is unchanged.
6. **Vignettes:** grep them for `bl_find_local_cf()`, and update any text that quotes a number
   that moves.
7. **The agreed paper sentences:** re-check every number they quote against the new run (A5's
   empty runs, A7's 0.14 to 0.42, B7). If a number moves, update the form, and I tell you.

### Step 4: the approximation-gap experiment (A1)

**New script:** `scripts/paper1-pima/05_full-space-gap.R`, run after `run_all.R`.
1. `source("00_setup.R")`; `bl_results <- fresh_result()`; the raw GAM from `reduced.rds`.
2. For each of the 79 test rows: `run_local()` with no filters. Record the search point and
   `sqrt(dist_mahalanobis)`, with `hull = TRUE` (the method) and also `hull = FALSE`, so the hull's
   effect over many rows is counted: points that move or vanish.
3. Full-space problem:
   - minimise (x − x₀)′W⁻¹(x − x₀) subject to f(x) = 0.5 and the training ranges as bounds;
   - f is the raw GAM probability (`mgcv::predict.gam(type = "response")`), not the floored
     3-dp prediction;
   - solved with `nloptr::slsqp()`, started from the package points and from 50 random-direction
     boundary crossings (`uniroot`; fixed seed);
   - the best converged solution with |f − 0.5| < 1e-6, within the bounds.
4. Ratios (the package's distance / the full-space distance) for `hull = TRUE` and for `FALSE`:
   the median, quartiles and range; the share at or below 1.05; by outcome group; and the TP and
   TN targets.
5. The hull check for the full-space point (A1.6(b)): whether it lies inside the four-dimensional
   convex hull of the training data (`geometry::convhulln`/`inhulln`;
   `install.packages("geometry")`). The ratio is also given for the plausible rows. A1.6(a) is
   covered by step 2's counts.
6. Output: `output/full_space_gap.csv` (one row per target) and `output/full_space_gap.rds` (the
   summary), not `numbers.rds`, which `run_all.R` rewrites. A new section in `review.qmd`, skipped
   if the file is missing.

**Error cases:**
- a target with no package solution: excluded and counted;
- every start fails: NA, counted;
- a ratio below 1 − 1e-6: flagged;
- a target outside the training ranges: the starts are clamped.

**Checks:** with `hull = FALSE`, the TP and TN rows reproduce the old search distances (1.5018,
0.1405); with `hull = TRUE`, they match the new Table 1. If the median ratio is well above 1, I
bring it to you before any paper sentence uses it.

### Step 5: release housekeeping (B24)

- `DESCRIPTION`: Version 0.2.0. NEWS.md: "development version" becomes "boundarylogic 0.2.0".
- Thesis citations:
  - the page citations in the rotation note (about 30), `boundary-distance-and-robustness.md` (3)
    and `2 implementation_summary.txt` (1) move to the final thesis
    (`../Explain XAI code/reference-materials/1 FINAL_ADRIAAN ROWAN_27022026.pdf`; printed page =
    PDF page − 16), each found by `pdftotext -layout` on both PDFs;
  - the methodology note's sentence on the draft is changed in Phase 4.
- `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`: expect 0/0/0.
- Rebuild pkgdown by the CLAUDE.md procedure (`clean_site`, `build_site`, unlink
  `docs/CLAUDE.html`/`progress.html` and their `.md`, delete the helper script).
- Merge and tag: only on your instruction, after the supervisor's Gate 3 decisions.

### After the implementation

- Copy this plan to `.claude/reference/` with an IMPLEMENTED banner.
- Stop with the PASS/FAIL count and the list of Paper 1 numbers that moved. No commit.

### Later (unchanged)

- **Phase 4:** after the supervisor's decisions, evolve `methodology-changes.qmd` into the single
  design-decisions document, with the future section built from list 3 above.
- **Phase 5:** retire `known-issues.md` (with your confirmation) and repoint its references.

## Verification

- B10:
  - `devtools::test()` passes;
  - with `hull = FALSE`, the Pima runs reproduce today's `table1.csv` exactly;
  - with `hull = TRUE`, every reported search point lies inside its plane's hull (checked in
    `05`);
  - `tests.R` passes.
- Step 2: the README output comes from a real run; the B9 warning fires once at the centroid.
- Step 3: the methodology note renders and its claims pass; `17_3`'s output is unchanged.
- Step 5: `devtools::check()` gives 0/0/0; pkgdown builds; the moved citations are checked against
  the final PDF's text.
