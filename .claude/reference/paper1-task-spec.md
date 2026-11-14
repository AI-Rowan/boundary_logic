> **ARCHIVED 2026-11-07.** The Paper 1 task specification (Tasks A-E), written by the Paper 1
> project (its original is the paper repo's `.claude/plans/t2-package-fix-prompt.md`). Tasks A-E are
> done in this repo; the open items (Gate 3, the `v0.2.0` tag, the hand-over to the paper writer, the
> paper-repo task T4) are in `progress.md`, "Outstanding items".

# Prompt for the Boundary_logic session: fixes, the Paper 1 example, and the method record

Date created: 2026-10-29
Date last edited: 2026-10-29

Written by the Paper 1 project (`Explain XAI code/reference-materials/Paper_nr1/Exploring_boundaries/`)
from its plan `t2-package-fix-prompt.md`. Paste everything below the line into a Claude Code session
opened in `Boundary_logic/`.

---

## 0. Context and ground rules

I am writing a paper for the *Journal of Computational and Graphical Statistics* (JCGS):
**"Decision-Boundary Biplots: Global and Local Counterfactual Interpretation of Machine Learning
Models"** (Paper 1). Every figure and number in it will come from this package. JCGS reviewers run
the code, so the package must behave predictably and be pinned to a tagged commit. This session has
five tasks:
- **A:** fix a graphics-device side effect;
- **B:** fix a target-class lookup;
- **C:** write and run the paper's example code here, so that any fault it exposes is fixed at the
  source;
- **D:** record the issues we defer;
- **E:** write a note for my supervisor on every methodology change since the thesis, plus possible
  future changes.

**Rules:**
- **This package's `CLAUDE.md` wins** wherever it conflicts with this prompt. Follow it throughout:
  - plan first, in `.claude/plans/`;
  - stay on `method_developments`;
  - R source in ASCII only;
  - verification scripts as `.R` files run with `Rscript`, never `Rscript -e`;
  - `bl_results` and `bl_shapley_values` as object names;
  - update `2 implementation_summary.txt` when an S3 field changes;
  - update the `development/reference/` walkthrough for any covered source file you change;
  - no new `.md` files in the package root (pkgdown publishes them).
- **Leave the existing uncommitted work alone.** At the time of writing (2026-10-29, `c6ca185`) the
  working tree has:
  - modified `2 implementation_summary.txt`, `CLAUDE.md` and `progress.md`;
  - a rename of `documentation/rotation_plane_underdetermination_note.md` to `development/`;
  - an untracked `documentation/biplot-rotation/`.

  Where your edits touch the same files, add them alongside without reverting anything, and list
  the overlap in your report. Compare against `git status` at the start before deleting anything.
- **Never commit, push or tag** unless I say so.
- **Environment on this machine:**
  - **Rscript:** run it through PowerShell, `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" file.R`.
    It segfaults under Git Bash.
  - **Quarto:** use RStudio's binary,
    `C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe`. Its knitr engine needs the
    `rmarkdown` package.
- **If a line reference below has moved,** find the function by name, and say so.

**Write one plan covering A to E,** then work through it with three review gates. At each gate,
stop and wait for me:
1. **Gate 1:** after Tasks A and B pass `devtools::document()` and `devtools::test()`. Report the
   PASS/FAIL counts.
2. **Gate 2:** after Task C's first run. I review the outputs and decide, issue by issue, whether
   to fix now or defer.
3. **Gate 3:** after the drafts of Tasks D and E.

---

## Task A: the graphics-device side effect, in two helpers

**The problem.** Two internal helpers each open a temporary `png()` and close it with a bare
`dev.off()`. After closing a device, R activates the *next* open device, not the one that was
active before. So when two or more devices are open, the exported functions that reach these
helpers silently switch the caller's active device. Later plots then land on the wrong device.

| helper | location | reached from |
|---|---|---|
| `.build_hull_polygon()` | `R/hull_utils.R:28` | `bl_prepare_data()` → `bl_filter_outliers()` (`R/outlier_filter.R:89`); `bl_build_result()` → `bl_build_grid()` (`R/biplot_grid.R:212`) |
| `.bl_get_plot_bounds()` | `R/biplot_grid.R:31` | `bl_build_result()` → `bl_build_grid()` (`R/biplot_grid.R:194`) |

The hull is built even at `hull_fraction = 1`, so turning outlier removal off does not avoid it.

**Reproduction (measured on 0.1.0):**

```r
library(boundarylogic)
loan_raw <- read.csv(system.file("extdata", "loan_data.csv", package = "boundarylogic"))
loan <- loan_raw[loan_raw$previous_loan_defaults_on_file == "No", c(10, 12, 14)]
pdf(NULL); A <- dev.cur()   # the caller's device
pdf(NULL); B <- dev.cur()   # a second device, active
bl_dat <- bl_prepare_data(data = loan, class_col = "loan_status",
                          feature_cols = c("loan_percent_income", "credit_score"),
                          train_fraction = 0.01, seed = 121L, hull_fraction = 1)
dev.cur() == B   # FALSE: the active device is now A
graphics.off()
```

**The fix, in both helpers.** It was verified on a patched copy of `.build_hull_polygon()`: the hull
coordinates were identical, and the active device was kept.

```r
  prev_dev <- grDevices::dev.cur()
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(tmp_file, width = 600L, height = 600L)   # keep each helper's own size
  hull_dev <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(hull_dev)
    if (prev_dev > 1L) grDevices::dev.set(prev_dev)
    unlink(tmp_file)
  }, add = TRUE)
```

Two details matter:
- **`dev.off(hull_dev)`** closes the temporary device and no other.
- **The `prev_dev > 1L` guard:** `dev.cur()` is 1 when nothing is open, and `dev.set(1)` would open
  a new device.

`.bl_get_plot_bounds()` imports only `png` and `dev.off` from `grDevices`. Keep its roxygen
`@importFrom` lines consistent, and single-line.

**Tests,** in a new `tests/testthat/test-graphics_device.R`. Use `pdf(NULL)` devices only, never
`windows()` or `x11()`, and close every device each test opens.
- With two devices open, each helper leaves `dev.cur()` unchanged and `dev.list()` as it was:
  ```r
  set.seed(1); Z <- matrix(rnorm(200), ncol = 2)
  grDevices::pdf(NULL); d1 <- grDevices::dev.cur()
  grDevices::pdf(NULL); d2 <- grDevices::dev.cur()
  on.exit(grDevices::graphics.off(), add = TRUE)
  .build_hull_polygon(Z, -4, 4, outlie = 0.9)
  expect_identical(grDevices::dev.cur(), d2)
  expect_identical(unname(grDevices::dev.list()), unname(c(d1, d2)))
  ```
- With no device open, the helpers leave none open (`expect_null(grDevices::dev.list())`).
- End to end: `bl_prepare_data()` and `bl_build_result()`, run with two devices open, keep
  `dev.cur()`. Use a small iris fixture, as in `tests/testthat/test-result.R`.

**NEWS entry:**

> * `bl_prepare_data()`, `bl_filter_outliers()`, `bl_build_grid()` and `bl_build_result()` no longer
>   change the active graphics device when more than one device is open. The internal helpers
>   `.build_hull_polygon()` and `.bl_get_plot_bounds()` now restore the caller's device after
>   closing their temporary `png()`.

**Walkthroughs:** update `development/reference/review_section4_to_6.md` (`biplot_grid.R`) and
`review_section7_to_8.md` (`hull_utils.R`) where they describe these lines.

---

## Task B: `bl_select_target(target = i, data = df)` colours the target by the wrong row's class

**The problem.**
- `bl_select_target()` accepts a row index into any data frame, but stores only `row_id`
  (`R/local_cf.R:134`; lookup at `:156-164`; structure at `:183-194`).
- `plot.bl_local_result()` (`R/local_cf.R:799`) then reads the true class as
  `bl_result$test_data[[class]][row_id]` (`:884-894`).
- Whenever `data` is not `test_data`, that is another observation. The target is drawn in the wrong
  confusion colour, with no warning; the label still shows `i`, so the plot looks right.
- **Measured on the loan data:** 221 training rows, `train_fraction = 0.01`. `test_data[i]`'s class
  differed from `train_data[i]`'s in **120 of 221** rows.

**Step 1, a worked example:** a note-and-script bundle in `documentation/select-target-class/`,
following the `documentation/biplot-rotation/` convention (a note, plus the script that reproduces
every number in it). Use Pima with the Paper 1 set-up (Task C, "Set-up"). Report:
- the class mismatch between `train_data[i]` and `test_data[i]`, counted over every training row;
- one figure that shows a wrongly coloured target.

If Pima shows no mismatch, say so; the loan numbers above are then the evidence.

**Step 2, the fix:**
- In `bl_select_target()`, resolve the class where the row is chosen. Store it as a new field,
  `class_obs`; it is `NA_integer_` when the data has no class column. Use
  `resp <- setdiff(names(src), var_names)`, the same rule the plot uses now, where `src` is the
  one-row source:
  - `data[row_id, , drop = FALSE]` for an index;
  - `target` itself for a data frame.
- In `plot.bl_local_result()`, replace the lookup with
  ```r
  target_true_class <- bl_target$class_obs
  if (is.null(target_true_class)) target_true_class <- NA_integer_   # objects built before the fix
  ```
- **The results:**
  - an index into `test_data` behaves as before;
  - an index into any other frame uses the right observation's class;
  - a one-row data-frame target that includes the class column now gets its confusion colour too;
  - older `bl_target` objects plot by prediction colour.

**Because `bl_target` gains a field:** update `2 implementation_summary.txt` (the S3 field list)
and `development/reference/review_section9_to_15.md`.

**Tests:** the package's first Phase 3 test file, `tests/testthat/test-local_cf.R`. Use a small
fixture in the style of `test-result.R`, choosing one whose `train_data[1]` and `test_data[1]`
differ in class. Test that:
- `bl_select_target(bl_results, target = 1L, data = bl_results$train_data)$class_obs` equals the
  class of `train_data[1, ]`;
- a one-row data-frame target, `train_data[1, ]`, gives the same `class_obs`;
- `bl_select_target(bl_results, target = 1L)` (the default, `test_data`) is unchanged.

**NEWS entry:**

> * `plot.bl_local_result()` now colours the target by its own true class. Previously it looked the
>   class up in `test_data` by row number, which gave another observation's class whenever
>   `bl_select_target()` indexed a different data frame. `bl_target` gains a `class_obs` field; a
>   one-row data frame target that includes the class column now gets its confusion colour too.

**→ Gate 1.**

---

## Task C: the Paper 1 example code, run and reviewed here

**Purpose:**
- Produce every output Paper 1 needs, using the package as it stands after A and B.
- I review the outputs here and decide what to fix in the package.
- Once the package is updated, committed and reinstalled, the same scripts are copied into the
  paper's repo and re-run against the installed package.

So the scripts must run unchanged in both places, apart from the one line that loads the package.

**Location:** `scripts/paper1-pima/`, with outputs in `scripts/paper1-pima/output/`. `scripts/` is
git-ignored and `.Rbuildignore`d, and this code must stay local: it must not reach the public repo.

### Set-up (fixed for the paper; do not change it without asking)

- **Data:** `inst/extdata/pima_diabetes.csv` (768 rows; class column `Outcome`).
  - Drop every row with a zero in Insulin, SkinThickness, Glucose or BMI, and assert that
    **392 rows** remain.
  - Drop BloodPressure through `bl_prepare_data(feature_cols = …)`, never only in the model,
    leaving p = 7: Pregnancies, Glucose, SkinThickness, Insulin, BMI, DiabetesPedigreeFunction and
    Age.
- **Split:** `bl_prepare_data(data, class_col = "Outcome", feature_cols = <7 features>,
  train_fraction = 0.8, seed = 121L, hull_fraction = 1)`.
- **Model:**
  - `mgcv::gam(class ~ s(Pregnancies) + s(Glucose) + … + s(Age), family = binomial)` on
    `train_data`, in raw clinical units, with the default basis and no interactions.
  - Wrap it with `bl_wrap_model()`: `model_type = "GAM"` if that path predicts correctly, otherwise
    `"custom"` with a `predict_fn` as this package's `CLAUDE.md` §9 describes. Say which you used.
  - Cutoff 0.5.
- **Projection:** `bl_build_result(bl_dat, bl_mod, method = "CVA")`. With a model passed, this is
  CVA on the four prediction outcomes (TP, TN, FP, FN), so K = 3.
- **Local search:** `max_pairs = 21L` (the default since 2026-10-31; it was `10L`) and `distance = "mahalanobis"` (the default).

### Files

| file | does |
|---|---|
| `00_setup.R` | Libraries, paths and the seed. **Package loading is the one line that differs between repos:** `devtools::load_all("<package root>")` here, `library(boundarylogic)` in the paper repo. Records the package version and the git SHA (or "installed") in the outputs. Defines `fresh_result()`, which rebuilds `bl_results` from the saved data and model objects. |
| `01_data-model.R` | the 392-row filter (asserted), `bl_prepare_data()`, the GAM, `bl_wrap_model()`; saves the objects to `output/` |
| `02_global.R` | `bl_results`; `bl_find_boundary(bl_results, data = train_data)` and `(…, data = test_data)`; global importance with `plot(bl_bnd, type = "boxplot")` and `bl_robustness()`; figures 3 and 4 |
| `03_local.R` | target choice, `run_local()`, all local runs; figures 5 to 9 and Table 1 |
| `fig-03.R` … `fig-09.R` | one base-R figure each, written as PDF and PNG to `output/`. Separate files, so the paper can include them with Quarto's `#| file:`. |
| `tests.R` | the assertions below; stops on the first failure |
| `review.qmd` | an HTML review page for me: the figures, Table 1, the assertion results, and the issues-found table |
| `output/numbers.rds` | every number the paper will quote: data sizes, confusion counts, global importance, target IDs and values, and for each local run `best_pair`, √d²_M, `B_x`, `B_pred`, `solution_found`, the Shapley values and the sparse counterfactual |

### Local runs (`03_local.R`)

- **Targets:** one test-set TP and one test-set TN (true class against prediction at 0.5). Take the
  one whose global counterfactual distance (from `bl_find_boundary(…, data = test_data)`) is nearest
  its group's median. This avoids hand-picking. Say which distance field you used.
- **Select targets from `test_data`.**
- **`run_local(target, filters, max_pairs)`** builds a **fresh** `bl_results` on every call, then
  runs:
  - `bl_select_target()`;
  - `set_filters()` (when there are filters);
  - `bl_find_local_cf()`;
  - `bl_shapley()`;
  - `bl_find_sparse_cf()`.

  Results depend on call history, and `best_pair` is not always reproducible. So quote numbers
  only from a first call on a fresh `bl_results`. Never "check" one by repeating the call.
- **Runs, for each target:**
  1. unconstrained;
  2. both cannot decrease: `set_filters(tgt, Pregnancies = "increase", Age = "increase")`;
  3. both fixed: `Pregnancies = "fixed", Age = "fixed"`;
  4. **mixed, the paper's main constrained example:** `Pregnancies = "fixed", Age = "increase"`;
  5. `max_pairs = 1L` (pins pair (1,2)), for figure 6(b).

### Figures and Table 1

| # | content | package objects |
|---|---|---|
| 3 | the global CVA biplot: four outcome groups, the approximated decision surface, the local decision boundaries | `plot(bl_results)` |
| 4 | (a) the training points with their nearest boundary points on the global biplot; (b) per-variable distance importance (boxplot) | `bl_find_boundary()`, `plot(bl_bnd, type = "boxplot")` |
| 5 | the TP target: (a) on the global plane; (b) on its Rotated Component biplot | `bl_results$biplot_grid$ct`, the target's `z_obs`; `bl_local$Vr_rot`, `bl_local$ct_local` |
| 6 | (a) a dot chart of √d²_M for each of the 21 pairs, both targets, the shortest marked; (b) the TP target on pair (1,2) beside its winning pair, only if they differ | `all_distances_mahalanobis` (stored squared); run 5 |
| 7 | the TP target, unconstrained (run 1) beside mixed (run 4), each on its own winning plane, never overlaid on one plane | `plot(bl_local)` |
| 8 | Shapley values towards the counterfactual (Supports/Contradicts), both targets | `plot(bl_shapley_values)` |
| 9 | the sparse counterfactual on the local biplot, both targets | `plot(bl_sparse)` |
| Table 1 | for each target: the observed values against the counterfactuals of runs 1 to 4, and the sparse counterfactual; `solution_found` and √d²_M per run | |

- Figures 1 and 2 are drawn schematics. They are made in the paper repo, not here.
- **Plot the package's own objects.** Never re-derive a plane or a projection by hand.
- **If a figure needs something the package cannot draw non-interactively** (4(a) may be such a
  case: boundary arrows exist only in the interactive `bl_pick_point()`), do not hand-roll it. Log
  it as an issue, with a proposed package change, for gate 2.

### Assertions (`tests.R`)

- 392 rows, p = 7, and all four outcome groups non-empty.
- Every counterfactual's `B_pred` is on the other side of 0.5 from the target's prediction.
- The Shapley values sum to `B_pred − pred_prob`, within numerical tolerance.
- **Constraints hold:**
  - "fixed" variables are within ±0.5 of their observed value;
  - "increase" variables are ≥ their observed value.
- In the sparse counterfactual, each fixed variable is back at its observed value.
- Every counterfactual lies within `train_ranges`.
- **Reproducibility:** run the whole pipeline twice, each time in a fresh R session. The two
  `numbers.rds` files are identical, apart from the version fields.

### The review page and the issues-found table

`review.qmd` renders to `review.html`. It shows each figure, Table 1, a pass/fail line for each
assertion, and a table of every issue found:

| column | content |
|---|---|
| symptom | what went wrong or looked wrong |
| where | file:line in the package |
| effect on the paper | which figure or number it touches |
| proposed fix | the change, and whether it changes results |
| decision | left blank, for me: fix now / defer |

**→ Gate 2.**
- For each issue I mark "fix now": a short plan, the fix, a test and a NEWS entry, then re-run Task C
  from a fresh session.
- Each issue I mark "defer" goes to Task D.

---

## Task D: record the deferred issues (no code changes)

**The note:** `documentation/known-issues.md`, with one entry per issue: symptom, file:line, impact
on results, and suggested fix. Also add **one bullet** to `CLAUDE.md` §9 ("Deferred / Future Work")
pointing to the note. §9 already holds the rotation under-determination, so the `best_pair` entry
refers to that bullet rather than repeating it.

**The issues to record:**
1. **`best_pair` instability:** `best_pair` is not always reproducible across identical calls.
2. **Call-history dependence:** the result depends on earlier calls (`bl_find_local_cf()`,
   `plot(bl_results)`, `bl_find_boundary()`, `bl_surrogate()`) in the same session.
3. **The approximation-gap experiment has not been done:** how far the 2D boundary lies from the
   p-dimensional one.
4. **Only the winning pair's plane is stored.** The other pairs keep only `all_distances*`.
5. **Plot colours are inconsistent:** `plot.bl_projection` uses a different scheme from the other
   biplots.
6. **The ±0.5 tolerance for "fixed" is hard-coded** (in `local_cf.R`).
7. **The cutoff is not passed to `bl_build_grid`.**
8. **`rotate_deg` is not applied to the sparse counterfactual's cross.**
9. **Stale documentation:**
   - roxygen: `dist_z`, the sparse colours, the Shapley `seed`;
   - the README;
   - the implementation summary's description of the split as "stratified".
10. **Every issue from Task C that I chose to defer.**

---

## Task E: methodology changes since the thesis, for my supervisor's approval

My supervisor must approve every change to the method made since the thesis before the paper
relies on it. Write a note that my supervisor can comment on in Word.

**Deliverable:** `documentation/methodology-changes/methodology-changes.qmd`, rendered to `.docx`
with the Quarto above. The rendered `.docx` stays next to the `.qmd`. Record in the note the
package commit at which it was measured.

**Baseline, "the thesis":**
- **The final thesis text:** `../Explain XAI code/reference-materials/1 FINAL_ADRIAAN ROWAN_27022026.pdf`
  (Stellenbosch, March 2026). Its chapters are available as markdown in
  `../Explain XAI code/reference-materials/Paper_nr1/thesis-md/`, chapters 07 to 11. Do **not** use the
  30 Sept 2025 draft in `documentation/`.
- **The thesis code:** the branch `origin/original_PhD_code`. Read it with `git show` and
  `git diff`, and never check it out.
- **If the thesis text and the thesis code disagree,** report both. The text defines "the thesis
  version".

**Sources for what changed:**
- `git log method_developments` (the package was started on 2026-01-15);
- `NEWS.md`;
- the design decisions in `2 implementation_summary.txt`;
- the notes in `documentation/`;
- the plan archives in `.claude/reference/`.

### Part 1: changes made

The starting list. **Verify each one against the sources, correct it where it is wrong, and add any
others you find:**
1. **Variable importance:** standardised by the total SD → standardised by the within-group SD,
   √diag(**W**).
2. **The metric W:** the pooled covariance of the four outcome groups, divided by (n − g).
3. **The pair search:** all p(p−1)/2 pairs (21 for p = 7) → the first `max_pairs = 10` in canonical
   order (1,2), (1,3), (2,3), (1,4), …, so components 6 and 7 are never reached when p = 7. The Paper 1
   plan recorded that a thesis example selected pair (1,7); check this. (2026-10-31: the default is now
   `max_pairs = 21L`, so the package tries all 21 pairs for p = 7; `10L` is the earlier default.)
4. **Components beyond K:** components 4 to 7 are an arbitrary completion, because four-group CVA
   has only K = 3 canonical axes.
5. **The feasibility filter:** the convex hull applied in the local search → no hull after the
   rotation; only `train_ranges` and `set_filters` constrain it.
6. **A new constraint type:** "cannot decrease" only → the added "fixed" type, which may move ±0.5
   during the search and then reverts in the sparse counterfactual.
7. **The contour band:** b ± 0.001, with predictions floor-rounded to 3 decimals.
8. **Re-scoring:** contour vertices are re-scored by the model.
9. **The `.bl_rotate()` fix** changed the results for pairs other than (1,2).

**For each change, give:**
- the thesis version, with its section;
- the package version, with the function;
- the reason for the change;
- **the effect on the Pima results.**
  - Measure it with Task C's set-up wherever the package can switch back: `distance = "euclidean"`
    for 1 and the pair selector; `max_pairs = 21L` against `10L` for 3.
  - Write the measurement script into `scripts/paper1-pima/` (local), and have it write a small
    `method-change-effects.csv` into `documentation/methodology-changes/` for the `.qmd` to read.
  - Where a change cannot be switched back, describe its effect and say plainly that it was not
    measured.
- **a decision box** for my supervisor: approve / query / reject, plus a comment.

**Exclude** what does not change the method:
- API renames, plotting and the list of models;
- the exact Shapley computed over 2^p subsets instead of p! orderings (the values are the same);
- two-group CVA (the thesis also used four outcome groups).

Add a one-line appendix of the bug fixes from Tasks A, B and C, so the record is complete.

### Part 2: possible future changes

For each one, give: what it is, why, its expected effect on the results, the effort, and your
recommendation. Include:
- the minimal Givens rotation with a plane-family search
  (`documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §6, a proposal);
- storing every pair's plane;
- the approximation-gap experiment;
- the full-Mahalanobis Shapley decomposition in Phase 2 (`CLAUDE.md` §9);
- aggregated Shapley values (global summaries of the local attributions);
- a configurable "fixed" tolerance;
- a CVA biplot for a singular **W** (p > n), by the generalised SVD of Ganey and Gardner-Lubbe (2026,
  *Statistical Methods & Applications* 35, 825–847, R package `wideRhino`);
- the LDA multi-class route (`CLAUDE.md` §9);
- per-variable standardisation (`CLAUDE.md` §9);
- the Task C issues I deferred.

**→ Gate 3.**

---

## A recommendation, then your report

**The recommendation.** Once I have approved and committed A to E, push `method_developments` and
tag a release. The paper can then pin it with `devtools::install_github("AI-Rowan/boundary_logic@<tag>")`.
Push and tag only when I say so.

**Report back with:**
- the `devtools::test()` PASS/FAIL counts;
- the files changed, and any overlap with the uncommitted work listed in §0;
- the issues found in Task C, with my decision on each;
- the paths of `scripts/paper1-pima/review.html` and `documentation/methodology-changes/methodology-changes.docx`;
- the commit SHA, once I commit;
- the reinstall command for this machine: `devtools::install("<package root>")`, or the
  `install_github` line once tagged.

## Error cases

| case | what to do |
|---|---|
| a line reference has moved | find the function by name, and say so |
| this prompt conflicts with the package's `CLAUDE.md` | the package's file wins; say where |
| the uncommitted work touches a file you edit | edit alongside it, revert nothing, and list the overlap |
| a device test misbehaves without a screen | use `pdf(NULL)` only |
| an outcome group is empty, so CVA fails | stop with a message; do not change the set-up |
| a constrained run returns `solution_found = FALSE` | report it in Table 1; do not re-pick the target |
| `best_pair` differs between two fresh runs | record it as an issue; never hide it by re-running |
| `s(Pregnancies)` fails (too few distinct values for the default basis) | lower its `k`, and flag it as a set-up change for my approval |
| a Task E change cannot be measured on Pima | describe it, and say it was not measured |
| the thesis text and `origin/original_PhD_code` disagree | report both; the text defines the thesis version |
