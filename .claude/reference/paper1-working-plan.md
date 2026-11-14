> **ARCHIVED 2026-11-07.** The working plan for Paper 1's package tasks (Part 0, Tasks A-E, with
> their implementation logs), 2026-10-29 to 2026-11-05. Every task is done; Task C was superseded by
> `paper1-pima-example-plan.md`. The open items (Gate 3, the tag, `method_developments`) are in
> `progress.md`, "Outstanding items".

# Paper 1 package changes: commit, branch, then Tasks A-E

Date created: 2026-10-29
Date last edited: 2026-11-04

## Implementation log

- **Part 0 (done 2026-10-29):** commit `b7aee65` was pushed to `origin/method_developments`.
  `paper_nr1_changes` was created at `b7aee65` and pushed with upstream. CLAUDE.md §10 records
  the branch; that edit is uncommitted.
- **Tasks A and B (done 2026-10-29). Changes from the plan:**
  - `class_obs` coercion is wrapped in `suppressWarnings()`, so a non-numeric label becomes `NA`
    without a coercion warning.
  - The bundle script also writes `select_target_class_results.txt`, the source of every number
    in the note.
  - The pre-fix run at `b7aee65` drew both panels orange, which confirms the bug. Its figure is kept
    only in the scratchpad.
  - The 9c passage in `review_section9_to_15.md` was wrong: `plot_biplotEZ()` draws `target_point`
    yellow, not in a confusion colour. It was corrected.
  - `review_section4_to_6.md` §4b held a broken sentence; it was rewritten.
  - The six line references in `review_section9_to_15.md` point into `project_points.R` and
    `shapley.R`, so the `local_cf.R` edits do not move them.
  - The `shapley.R` line references there were already stale. They are fixed in the Gate 1
    follow-up (P3).
  - Tests: 158 PASS / 0 FAIL / 4 WARN. The warnings are the same CVA 2-class notices as before.
  - The harness `03_loan_status_Boundary_Logic.R` ran end to end in 0.9 min, SHAP comparison
    included.
  - The harness warns 6 times that `person_age` is not found, from `label_offset_var` and
    `ticks_var` in the reduced-model plots, because `person_age` was pruned in `f540607`. The root
    cause is in the package helpers; it is fixed in the Gate 1 follow-up (P2).

- **Gate 1 follow-up, P1-P3 (done 2026-10-29). Changes from the plan:**
  - The length-mismatch warning fires only when at least one listed variable is valid. This keeps
    `label_offset_var = 0L` with a vector of offsets silent.
  - The harness comment names the "Prune least-important variables" block. Step 4b was the wrong
    label, because `vars_to_drop` is set before Step 4b.
  - `review_section9_to_15.md` 9c gained a sentence on the pairing rule.
  - The older `set_filters(person_age = "fixed")` examples in `review_section9_to_15.md` (:164,
    :198, :731) are outside P1-P3. They are logged for the Task D stale-docs item.
  - Tests: 176 PASS / 0 FAIL / 4 WARN, of which 18 are the new `test-plot_helpers.R`.
  - The harness re-ran clean in 1.0 min, with no warnings.

- **Gate 1 follow-up 2, the S1-S6 sweep (done 2026-10-29). Changes from the plan:**
  - **Small additions:**
    - `plot.bl_boundary()` printed "25%%" literally in its note. The `cat()` string is fixed, with a
      NEWS entry.
    - A stale code comment in `predict_utils.R` claimed "GBM fitted by `bl_fit_model()`"; it is
      fixed.
  - **Walkthrough restructuring:**
    - `review_section9_to_15.md`: `bl_select_target()` moved from 9b into Step 10 (as 10a,
      matching the script); 9c became 9b, and `set_filters` became 10b.
    - `review_section4_to_6.md`: the subsections became 6.1-6.6, and §4b was folded into Step 4.
  - **Real output:** all console samples come from harness runs (`harness_log3.txt`).
  - **The harness script:** the header and comments are fixed (Step 3 population, CVA,
    `bl_pick_point` with `data =`, "four variables"). The unused `shap_top5` / `dist_top5`
    variables are left in place (no code change).
  - **Not stale, so left alone:**
    - the pre-existing non-ASCII characters in `R/` comments (32; allowed by `R CMD check`);
    - the harness title string "Loan default (prior defaulters)", which is a code string.
  - **Results:**
    - tests: 176 / 0 / 4;
    - the README chunks all run;
    - the README Pima constraint example returns `solution_found = FALSE` for its first predicted-1
      target (no error);
    - the harness ran clean.
- **Committed 2026-10-29** as `f2a076b` on `paper_nr1_changes` (not pushed), with CLAUDE.md and
  progress.md updated for the handoff.
- **Next:** the user reviews the three `review_section*.md` walkthroughs; Task C starts only on
  their approval.
- **Sessions 2-7** (Task C, the Gate 2 fixes, Task D) are logged in `progress.md` and in the
  archived plans in `.claude/reference/`.
- **Task E (drafted 2026-11-04, session 8; plan
  `.claude/reference/methodology-changes-taskE-plan.md`). Changes from this plan:**
  - **Structure (user, at plan review):** a history of changes, then the core methodology (one
    decision box per change), then the code improvements summarised (one box for the section). This
    replaces the one-line bug-fix appendix.
  - **Three baselines:** the thesis text, its Appendix A, and the branch, which was committed on
    2026-03-09, after submission. They differ, for example on the tie-break and on the data behind the
    local hull.
  - **The pair selector** is reported as the CVA identity, confirmed on Pima (ratio n-g = 309 for
    every pair). It is not measured under PCA.
  - **`max_pairs`** is measured on the 4-feature model only (10 = 21; all 6 pairs against pair (1,2)).
  - **Measurements:** `scripts/paper1-pima/04_method-change-effects.R` writes
    `documentation/methodology-changes/method-change-effects.csv`.
    - It has five self-checks. In one respect it departs from the plan: each configuration shares one
      fresh `bl_results` across its runs, as `03_local.R` does, instead of building one per run.
  - The note is `methodology-changes.qmd`, rendered to `.docx`.

## Context

Paper 1 (JCGS, "Decision-Boundary Biplots") will take every figure and number from this package, and
reviewers will run the code. `paper1-task-spec.md` (written by the Paper 1 project) specifies
five tasks:
- **A:** fix a graphics-device side effect.
- **B:** fix a target-class lookup.
- **C:** run the paper's Pima example here.
- **D:** record the deferred issues.
- **E:** write a methodology-changes note for the supervisor.

The user has added two steps that change the prompt's ground rules:
1. Commit and push the current uncommitted documentation work on `method_developments` first. The
   prompt's "leave the uncommitted work alone" rule is therefore moot.
2. Do all of A-E on a new branch, `paper_nr1_changes`, not on `method_developments`. This overrides
   the prompt's "stay on `method_developments`" and CLAUDE.md §10's "start every session on
   `method_developments`".

Everything else in the prompt stands, including its three review gates. Where the prompt and the
package's CLAUDE.md conflict, CLAUDE.md wins.

**Pre-checks run during planning:**
- **Line references:** every line reference in the prompt is still exact. None has moved.
- **Branches:** `origin/original_PhD_code` exists. `paper_nr1_changes` does not exist, locally or on
  the remote. `method_developments` and `origin/main` both sit at `c6ca185`. There are no tags.
- **Pima data:** the zero filter leaves exactly **392** rows (262 class 0, 130 class 1).
  Pregnancies has 17 distinct values, so the default `s()` basis (k = 10) fits.
- **Model and tools:** `bl_wrap_model()` already accepts `model_type = "GAM"` (`R/model_fit.R:219`).
  `rmarkdown`, `knitr` and `mgcv` are installed. There is no renv. Quarto is at the RStudio path.

---

## Part 0: commit, push, branch (the only git writes in this plan)

1. **Move the prompt file out of the root.** Per the user's decision:
   `mv paper1-task-spec.md `.claude/reference/paper1-task-spec.md`. `.claude/plans/` is
   gitignored, so the file is never committed or published.
2. **Stage explicit paths only, never `git add -A` on the whole tree:**
   ```bash
   git add -- "2 implementation_summary.txt" CLAUDE.md progress.md documentation/biplot-rotation
   git add -A -- development/rotation_plane_underdetermination_note.md
   ```
   The second line stages the deletion. Against HEAD, the net effect is that
   `documentation/rotation_plane_underdetermination_note.md` is deleted. Its content now lives in
   §4.1 of `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`.
3. **Check what is staged with `git status`.** Expect:
   - 3 modified files;
   - 1 deleted file;
   - 6 new files under `documentation/biplot-rotation/`;
   - nothing untracked except gitignored paths.
4. **Commit**, using a heredoc for the message:
   ```
   Add biplot-rotation note bundle; drop superseded rotation-plane note

   Adds documentation/biplot-rotation/: the technical note on .bl_rotate() and
   the eigenvector-pair search, plus the four scripts that reproduce its
   measurements. Removes documentation/rotation_plane_underdetermination_note.md,
   whose content the new note's section 4.1 now carries, and points CLAUDE.md,
   the implementation summary and progress.md at the new note. Documentation
   only; no R/ change.

   Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
   ```
5. **Push** with `git push origin method_developments`. `main` is not touched.
6. **Create and push the branch:** `git checkout -b paper_nr1_changes`, then
   `git push -u origin paper_nr1_changes`. The second command follows the user's decision.
7. **Record the branch in CLAUDE.md §10, on the new branch, uncommitted.** Add a `paper_nr1_changes`
   bullet: the Paper 1 fixes, branched from `method_developments` at `<sha>` on 2026-10-29. Also
   qualify the "always start on `method_developments`" rule: while Paper 1 work is open, start
   sessions on `paper_nr1_changes`. Without this, the next session's `git checkout
   method_developments` would leave the work.

From here on, **no commit, push, merge or tag** happens without an explicit instruction.

**Environment conventions for every step below:**
- **R:** run R as `.R` files only. Invoke them through the PowerShell tool:
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <file>.R`. The prompt reports that Rscript crashes
  under Git Bash, and CLAUDE.md bans `Rscript -e`. No `.ps1` files are written. Temporary driver
  scripts go in the scratchpad.
- **Quarto:** use `C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe`.
- **R source:** ASCII only.
- **Object names in scripts:** `bl_results` and `bl_shapley_values`.

---

## Files to read first

- **Task A:**
  - `R/hull_utils.R:22-56`: `.build_hull_polygon()`. `png()` is at :31; the bare `dev.off()` is
    at :35.
  - `R/biplot_grid.R:26-51`: `.bl_get_plot_bounds()`. `png()` is at :33; `dev.off()` is at :35.
  - The callers: `R/biplot_grid.R:194` and `:211-213`, and `R/outlier_filter.R:89`.
- **Task B:**
  - `R/local_cf.R:120-195`: `bl_select_target()` and its roxygen.
  - `R/local_cf.R:799-1010`: `plot.bl_local_result()`. The class lookup is at :884-894; the
    colouring is at :895-905.
  - `R/shapley.R:236`: `row_id` is passed through there.
- **Fixture pattern:** `tests/testthat/test-result.R`, and the `build_shap_pipeline()` helper in
  `tests/testthat/test-set_scaling.R`.
- **Doc files to update:** `NEWS.md` (the development-version sections), `2 implementation_summary.txt`
  (§3 at :79 and §4.8 at :322-326), and the three files in `development/reference/`.
- **Harness script:** `development/03_loan_status_Boundary_Logic.R`, especially Steps 4-6
  (:152-260), Steps 9-10 (:438-495) and Step 16 (:614-645).
- **Task E:**
  - the thesis chapters 07-11 in `../Explain XAI code/reference-materials/Paper_nr1/thesis-md/`;
  - the thesis code, read with `git show origin/original_PhD_code:"scripts/<file>"`;
  - `NEWS.md`, `2 implementation_summary.txt` §4.2.1, §4.8 and §6;
  - `documentation/*.md`;
  - `.claude/reference/*plan*.md`.

---

## Task A: restore the caller's graphics device

**Change.** In both helpers, replace the bare `png()` + `on.exit(dev.off())` with this, keeping each
helper's own size (hull 600x600, bounds 500x500):

```r
  prev_dev <- grDevices::dev.cur()
  tmp_file <- tempfile(fileext = ".png")
  grDevices::png(tmp_file, width = 600L, height = 600L)
  hull_dev <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(hull_dev)
    if (prev_dev > 1L) grDevices::dev.set(prev_dev)
    unlink(tmp_file)
  }, add = TRUE)
```

- **roxygen:** update to `#' @importFrom grDevices png dev.off dev.cur dev.set`, on one line, in both
  files.
- **NAMESPACE:** regenerated by `devtools::document()`.
- **Geometry:** nothing else changes, so the hull and bounds geometry are byte-identical.

**Tests** go in the new `tests/testthat/test-graphics_device.R`. Use `pdf(NULL)` only and
`on.exit(graphics.off())` in every test.
1. `.build_hull_polygon()` with two devices open: `dev.cur()` and `dev.list()` are unchanged.
2. `.bl_get_plot_bounds()`, the same check. Its biplot object comes from an iris
   `bl_build_projection()`; confirm the field name when implementing.
3. Both helpers with no device open: `expect_null(dev.list())` afterwards.
4. End to end: `bl_prepare_data()` and `bl_build_result()` on the iris fixture with two devices open
   keep `dev.cur()`.
5. The hull polygon coordinates are identical whether or not other devices are open.

**Docs:**
- **NEWS:** add the prompt's bullet under "Bug fixes and other improvements".
- **Walkthroughs:** these follow the rules in "Walkthrough and harness-script rules" below.
  - **`review_section4_to_6.md`, Step 4 (`bl_prepare_data()`):** explain how the outlier hull is
    built. `.build_hull_polygon()` peels `aplpack::plothulls()` at `hull_fraction`; at
    `hull_fraction = 1` the hull is still built but removes nothing.
  - **`review_section4_to_6.md`, Sub-step 6b, items 3 and 6:**
    - Explain why the bounds come from an off-screen render: biplotEZ fixes the axis extents only
      when it draws, and `par("usr")` gives the square window that sets the grid's range.
    - Explain the device handling: the helper opens a temporary `png()`, closes only that device,
      and re-activates the caller's device, so the user's open plots are untouched.
  - **`review_section7_to_8.md`, Stage 2 (clip to the convex hull):** explain where the polygon
    comes from (`bl_grid`, built by `.build_hull_polygon()` in final Z-space) and the same device
    handling.

## Task B: the target's true class

**Step 1, the worked example. Run it before the fix, so the bug is shown genuinely.**
- **Location:** `documentation/select-target-class/` holds `select_target_class_note.md`,
  `select_target_class_example.R` and its figure PNG. This follows the `biplot-rotation/` bundle
  convention.
- **Set-up:** the script is self-contained, using the Task C Pima set-up (filter, prepare, GAM,
  wrap, CVA).
- **Mismatch count:** count the class mismatch between `train_data[i]` and `test_data[i]` for every
  training row. Pima has about 313 training rows and about 79 test rows, so the count has two parts:
  - for `i <= n_test`: how many classes differ;
  - for `i > n_test`: those targets fall back to prediction colour, because the lookup is
    out of range.

  Report both parts. If Pima shows no mismatch, say so, and cite the loan numbers in the prompt.
- **Figure:** one target with two panels.
  - **Left:** coloured by `test_data[i]`'s class. This is the pre-fix behaviour, confirmed against
    a genuine pre-fix run.
  - **Right:** coloured by its own class.

  After the fix, the script reproduces the left panel by setting `class_obs` on a copy of the
  target, and the note says so. This keeps every number and figure reproducible after the fix.

**Step 2, the fix (`R/local_cf.R`):**
- **`bl_select_target()`.** Build the one-row source: `src <- target` for a data frame, or
  `data[row_id, , drop = FALSE]` for an index. Then:
  ```r
  resp <- setdiff(names(src), var_names)
  if ("class" %in% resp) resp <- "class"          # package convention; see note
  class_obs <- if (length(resp) >= 1L)
    as.integer(as.character(src[[resp[1L]]])) else NA_integer_
  ```
  Add `class_obs` to the returned list. Update `@return` to list `class_obs`, and also `scaling`,
  which it omits today.
  - **One small refinement on the prompt, for your approval:** if a `"class"` column is present, use
    it. Otherwise use the first non-feature column, as the plot does. `bl_prepare_data()` and
    `bl_wrap_data()` both name the column `class`. An external frame with an ID column ahead of its
    class column would otherwise pick up the ID.
- **`plot.bl_local_result()`:** replace :884-894 with
  `target_true_class <- bl_target$class_obs; if (is.null(target_true_class)) target_true_class <- NA_integer_`.
- **Unchanged:** the target label, `print.bl_target()` and `shapley.R`.

**Tests** go in the new `tests/testthat/test-local_cf.R`, the package's first Phase 3 test file.
- **Fixture:** iris versicolor, GLM, PCA via `bl_build_result()`, as in `build_shap_pipeline()`.
  Choose the seed so that `train_data[1]` and `test_data[1]` differ in class, and assert that as a
  precondition.
- **Tests:**
  1. An index into `train_data` gives the class of `train_data[1, ]`.
  2. The one-row data-frame target `train_data[1, ]` gives the same `class_obs`.
  3. The default (`test_data`) is unchanged: `class_obs` equals `test_data` row 1's class.
  4. A one-row data frame without a class column gives `NA_integer_`.
  5. `plot()` of a local result whose target comes from `train_data` runs without error on
     `pdf(NULL)`.

**Docs:**
- **NEWS:** add the prompt's bullet.
- **`2 implementation_summary.txt`:** in §4.8 (:326) note the `class_obs` field and where it is
  resolved; in §3 (:79) list it.
- **`review_section9_to_15.md`:**
  - **9b:** add `class_obs` and `scaling` to the field table (:86-94).
    - Explain the method: the true class belongs to the observation, so it is read from the same row
      (or the same one-row frame) that supplies `x_obs`.
    - Explain that `data =` may be any frame.
    - Explain that a target with no class column gets `NA`.
  - **9c (:133) and the Step 12 layer table (:410):** explain the colouring rule. With a known class,
    the colour shows the model's correctness for the target (TP / TN / FP / FN, from `class_obs`
    against `pred_class` at the cutoff). Without one, it shows the predicted class only.
  - **Object flow (:684-693):** add `class_obs`.
  - **Line references:** refresh the 6 `R/*.R:<line>` references in this file wherever the edits to
    `local_cf.R` shift them.

### Walkthrough and harness-script rules (apply to A, B, and every Gate 2 fix)

- **Updating the walkthroughs:** `development/reference/review_section*.md` are updated for every
  change made, covering both the code and the methodology the code captures.
  - Each passage says what the code does, with the function and file, and why: the methodological
    reason.
  - They describe the current code and method only. They carry no history, no "previously", and no
    "new" / "legacy" framing; NEWS.md is the change log.
  - Sections touched by this work are rewritten in that present-tense form. Untouched sections are
    left as they are.
- **Updating the harness:** `development/03_loan_status_Boundary_Logic.R` is updated for the changes
  before Task C starts.
  - **Step 10 block** (`tgt <- bl_select_target(...)`, :485): comments explain that the target's
    true class is resolved at selection time (`class_obs`) and drives the TP / TN / FP / FN colour of
    the local biplot.
  - **Step 10 block, example:** add a commented example of selecting a target from `train_data` via
    `data =`.
  - **Step 16 (external applicant, :634):** a comment explains that a frame with no class column
    gives `class_obs = NA`, so the plot colours by prediction.
  - **Steps 4-6:** a one-line comment that building the hull and the bounds no longer changes the
    active graphics device.
  - Existing object names stay as they are; only lines touched by this work are renamed to the
    CLAUDE.md conventions.
- **Running the harness end to end** is part of Gate 1. Use a scratchpad driver that:
  - sets the working directory to the package root, which the script's `devtools::load_all()` needs;
  - opens `pdf("<scratchpad>/harness_plots.pdf")` first, so no `Rplots.pdf` lands in the repo root;
  - `source()`s the script.

  Run it in the background, because it fits XGB and runs SHAP. Report any error and the step it
  occurs in.

### Gate 1
- Run `devtools::document()` and then `devtools::test()` from a scratchpad `.R` driver.
- Report the PASS / FAIL / WARN counts against the 136 / 0 / 4 baseline.
- Confirm that the walkthroughs and the harness script are updated, and that the harness ran end to
  end.
- Report the files changed. Stop and wait for the user. **Task C does not start until the user
  approves this gate.**

---

## Gate 1 follow-up: fix the pre-existing problems (P1-P3)

The user asked for this on 2026-10-29, after Gate 1 and before Task C. Their choices: fix both the
script and the package helpers (P2), and replace line numbers with code anchors (P3). Everything
follows the walkthrough and harness-script rules above. Task C still waits for the user's approval
afterwards.

### P1: the target-colour claims

- **Already done at Gate 1:** `review_section9_to_15.md` 9c now says `plot_biplotEZ()` draws
  `target_point` yellow (`R/plot_biplot.R:379-384`).
- **Still to fix,** missed in Task B: `2 implementation_summary.txt:396` (§4.8, the description of
  `plot.bl_local_result()`) still says the target's class "can be determined from
  `bl_result$test_data`".
  - Rewrite it as: the colour comes from `bl_target$class_obs` against `pred_class`, and is
    red/blue by prediction when `class_obs` is `NA`.
  - In the same sentence, correct `no_points` to `plot_points`.
- **No other occurrences:** a grep of `R/`, `vignettes/`, `README.md`, the walkthroughs and the
  summary found no other wrong claim.

### P2: the offset and tick pairing

**Root cause** (`R/plot_biplot.R:9-43`): `.make_label_line_vec()` and `.make_ticks_vec()` match
character names to `var_names` and drop the unmatched ones with a warning. They then pair the
values with the surviving variables *by position* (`rep_len(values, length(valid_idx))`). An
unknown name's value is not dropped, so every later value shifts one variable along.

**Symptoms in the harness** (`person_age` was pruned in Step 4b):
- **Steps 8b, 9 and 15 (both plots):** `person_age` is still named. Each produces the warning, and
  every value is shifted. For example, Step 15 gives `loan_amnt` 2 ticks instead of 4, and
  `credit_score` 200 instead of 20.
- **Steps 12 and 14:** `person_age` is commented out of the name vectors, but its `0` is still the
  first entry of `label_offset_dist = c(0,0,0,1,0)`. The offset meant for `loan_percent_income`
  silently lands on `credit_score`.

**P2a, the package fix (display-only).** In both helpers:
```r
n_var <- length(label_offset_var)                 # before any names are dropped
if (length(label_offset_dist) != 1L && length(label_offset_dist) != n_var)
  warning(sprintf(paste0("label_offset_dist has %d values for %d entries in ",
                         "label_offset_var; values are recycled."),
                  length(label_offset_dist), n_var), call. = FALSE)
if (is.character(label_offset_var)) {
  idx <- match(label_offset_var, var_names)
  ...                                             # existing "not found" warning
  if (length(label_offset_dist) == n_var)         # keep each value with its own variable
    label_offset_dist <- label_offset_dist[!is.na(idx)]
  label_offset_var <- idx[!is.na(idx)]
}
```
- `.make_ticks_vec()` gets the same change for `ticks_var` / `ticks_n`.
- A single value still applies to every listed variable.
- The defaults (`0L`, `1.5`) produce no warning.
- Integer index vectors are unchanged.
- **Vignette checked:** its reduced model keeps `person_age` (`vars_to_drop` at :335), and every
  call pairs equal-length vectors, so its output is unchanged.

**P2a, the rest of the package work:**
- **roxygen:** the `@param label_offset_dist` and `@param ticks_n` entries in `R/plot_biplot.R`,
  `R/local_cf.R` and `R/surrogate.R` gain one sentence each. It says: a vector pairs position by
  position with its variable vector; a single value applies to all; an unknown name is dropped,
  with a warning, together with its value.
- **NEWS** (bug fixes): "`label_offset_dist` and `ticks_n` now stay paired with their own
  variables when `label_offset_var` or `ticks_var` names a variable the model does not have.
  Previously the unknown name was dropped but its value was not, so every later value shifted one
  variable along. A warning now also flags a value vector whose length is neither 1 nor the length
  of its variable vector."
- **Tests** in a new `tests/testthat/test-plot_helpers.R`. Call the helpers via
  `boundarylogic:::`, with `var_names = c("a", "b", "c")`.
  1. All names valid: the values land on their positions, and the rest are 0.1 (labels) or
     `ticks_v` (ticks).
  2. An unknown first name with an equal-length value vector: `c("zz","b","c")` and `c(9,1,2)`
     give `c(0.1, 1, 2)`, with the "not found" warning and no other warning.
  3. A scalar value applies to every listed variable.
  4. A length mismatch, such as 3 values for 2 names, warns "values are recycled".
  5. The defaults (`0L`, `1.5`) give no warning and all 0.1.
  6. Tests 2-4 repeated for `.make_ticks_vec()`.

**P2b, the harness script** (`development/03_loan_status_Boundary_Logic.R`):
- **Each fixed vector:**
  - remove `person_age` from every name vector and remove its value from the paired vector;
  - drop the commented-out `#"person_age",` and `#2,` lines in Steps 12 and 14 as well.
- **The intent restored:** this restores the pairing each vector was written for, with
  `person_age` first.
- **The new vectors:**
  - Steps 8b and 9: `label_offset_dist = c(0, 0.5, 0, 0)`, so `loan_int_rate` gets 0.5.
  - Steps 12 and 14: `label_offset_dist = c(0, 0, 1, 0)`, so `loan_percent_income` gets 1;
    `ticks_n` stays `c(4,5,5,5)`.
  - Step 15, both plots: `label_offset_dist = c(0, 0, 0, 1)` and `ticks_n = c(4, 200, 200, 20)`.
- **Assumption:** these vectors follow the original intent. Today's rendering of Steps 12 and 14
  (the offset on `credit_score`) is taken to be the bug, not a deliberate choice.
- **A comment:** add one comment where Phase 2 (second pass) starts. It says `person_age` was
  pruned in Step 4b, and that each `*_dist` / `ticks_n` vector pairs position by position with its
  variable vector.

**P2c, the walkthroughs:**
- **`review_section4_to_6.md` Step 7, "Parameter notes":** add the pairing rule (position by
  position, a scalar applies to all, unknown names dropped with their values, a length mismatch
  warns). Rewrite that touched list and its intro line in present tense, removing "(replaced
  `no_points`)", "is now valid" and "now works".
- **`review_section9_to_15.md` 9c:** update the code block (about :131-137) to match the fixed
  Step 9 call.

### P3: stale references in `review_section9_to_15.md`

Replace the line numbers with code anchors, in present tense:

| where | stale | new text |
|---|---|---|
| :51 | "scripts 00, 01, 03, 05, 06 and both vignettes" | "the development harness (`development/03_loan_status_Boundary_Logic.R`, Step 9) and the loan-default vignette" |
| :51 | the `scripts/02_contour_inspection.R:60` sentence (the file no longer exists) | removed |
| :53 | `R/project_points.R:229` (still correct) | "inside `bl_predict()` (`R/project_points.R`)" |
| :499 | "Method selection (line 178)" | "the `if (p <= exact_max_vars)` branch of `bl_shapley()`" |
| :536 | "Contribute classification (line 206)" | "the `df$Contribute <- dplyr::case_when(...)` step of `bl_shapley()`" |
| :615 | "Algorithm (lines 364-385)" | "inside `bl_find_sparse_cf()`" |
| :621 | "(lines 377-385)" | "the `\"fixed\"` loop over `bl_local$set_filters` in `bl_find_sparse_cf()`" |

`review_section4_to_6.md` has no stale line references; its one match is prose.

### Gate 1b

- Run `devtools::document()` and `devtools::test()`. Expect 158 plus the new helper tests passing,
  0 failures, and the same 4 warnings.
- Re-run the harness through `run_harness.R` in the scratchpad. Expect no "not found" and no
  "values are recycled" warnings.
- Report the PASS / FAIL counts and the files changed. Stop and wait. Task C starts only when the
  user says so.

---

## Gate 1 follow-up 2: the stale-documentation sweep (S1-S6)

The user asked for this on 2026-10-29, after Gate 1b: remove all remaining stale issues. Three
read-only audits found about 150 stale passages. The user's choices:
- **Scope:** all five groups (walkthroughs and script, implementation summary, roxygen and
  messages, README / vignette / NEWS / CLAUDE.md).
- **Numbering:** renumber the walkthrough to match the script.
- **DESCRIPTION:** add the base-package Imports.

**Ground rules:**
- **Verify before editing.** Every audit finding is checked against the code before it is edited.
  A finding that does not hold is skipped and listed in the report; one audit, for example,
  misstated p for the loan model.
- **Code behaviour stays out.** Behavioural bugs stay for Gate 2 / Task D, not this sweep: the
  `used_in_sparse` flag after a fixed revert, the sparse cross ignoring `rotate_deg`, `cutoff` not
  passed to the grid, `xgboost::xgboost()` models failing in `.pred_function()`, and no guard for an
  empty outcome group. Where a document describes such a bug, it now states the current behaviour
  plainly.
- **Present tense everywhere.** The walkthroughs and the implementation summary describe current
  code; NEWS remains the change log.
- **Real numbers.** Every count, value and console sample in the walkthroughs comes from a real run:
  the harness log (`harness_log2.txt` in the scratchpad), or a fresh run if the log lacks it. No
  hand-written console output.

### S1: DESCRIPTION
- Add `grDevices`, `graphics`, `stats` and `utils` to Imports, because NAMESPACE imports from all
  four.

### S2: roxygen and user-facing messages (`R/*.R`, then `devtools::document()`)

| file | fix |
|---|---|
| `local_cf.R` ~358-361 | `dist_z` / `dist_mahalanobis` belong to `best_pair` only; the per-pair values are `all_distances*`; `dist_z` is the selector only when `distance = "euclidean"`; `NA` when there is no solution |
| `local_cf.R` ~296-297 | the target lies in the plotted plane along its projection onto the pair, at full length; align with the Rotation section at :303-305 |
| `local_cf.R` ~302 | `scripts/1.3 Optimal Rotation.R` exists only on the `original_PhD_code` branch; say so |
| `local_cf.R` ~757-759 | the `plot.bl_local_result` target colour comes from `bl_target$class_obs` (missed in Task B) |
| `shapley.R` ~453-454 | the sparse cross is `green3` (valid) / `yellow2` (invalid) |
| `shapley.R` ~150-152 | `seed` is used only when p > `exact_max_vars`; a non-NULL seed calls `set.seed()` |
| `shapley.R` ~254-255 | `plot.bl_shapley` returns the ggplot visibly |
| `projection.R` `@return` | add `cva_classes`, `point_col` and `train_ranges` |
| `biplot_grid.R` `@return`, `@param outlie` | add `hull_fraction` and `b_margin`; the grid vectors are hull-trimmed when `calc_hull = TRUE`; the `outlie` default is `1` |
| `result.R` | <ul><li>Value: add `train_ranges`; `test_data` is not filtered; `polygon`/`hull_fraction` is the grid hull (`outlie`).</li><li>`bl_build_result` returns a model-free `bl_result`, not a `bl_projection` (:206-209, :237-239, :265-266).</li><li>`bl_prepare_data` returns `bl_filter_result` (:38-42, :234-236).</li><li>The workflow section (:20-21, :34-35) drops `bl_filter_outliers()` after `bl_prepare_data()`.</li><li>:11: Phase 3, not Phase 2.</li></ul> |
| `boundary_plot.R` `@return` | add `distance` |
| `data_prepare.R`, `outlier_filter.R` | <ul><li>`@return` gains `scaling`.</li><li>`train_fraction` is in (0, 1].</li><li>`bl_filter_outliers(bl_data)` accepts a `bl_data` (from `bl_wrap_data()`).</li><li>The hull uses the PCs of the *standardised* features.</li></ul> |
| `pick_point.R` ~24-27 example | pass `data = bl_results$test_data` to both `bl_find_boundary()` and `bl_pick_point()` |
| `model_fit.R` ~160-163 | XGB means an `xgb.Booster` from `xgboost::xgb.train()`; `xgboost::xgboost()` models are not supported |
| `predict_utils.R` ~50-56 | GBM is wrapped, not fitted, by `bl_fit_model()` |
| `datasets.R` ~50-51 | the vignette keeps applicants with *no* prior default |
| **Messages** | <ul><li>`outlier_filter.R:62`, `result.R:98-99` and `:298-300`: `bl_data` comes from `bl_wrap_data()` / `bl_prepare_data()` gives `bl_filter_result`.</li><li>`plot_biplot.R:273-274`: drop "did bl_build_result return NULL?".</li><li>`print.bl_result` / `summary.bl_result`: label the grid's `outlie` as the grid hull, not the outlier-filter hull fraction.</li><li>Grep the tests for the old strings first.</li></ul> |

### S3: `2 implementation_summary.txt`

Apply all 31 audit findings. Keep the section numbers §4.2.1 and §4.6.1, because the code cites
them.
- **Signatures:**
  - `bl_assemble` (L203);
  - `plot_biplotEZ` (L214);
  - `bl_shapley` `seed=NULL` (L407);
  - the `plot.bl_local_result` shared and local-only parameter list (L396).
- **Returns:** `bl_build_result` always returns a `bl_result` (L204, L208); `bl_predict` returns a
  data.frame (L72).
- **Facts:**
  - the split is seeded simple random, not stratified (L103);
  - RForest (L124, L473);
  - the 7-step `plot_biplotEZ` pipeline, with no boundary overlay (L216-225);
  - the `bl_projection` fields, with no `Z_train` (L139);
  - the polygon is not visual-only (L193);
  - raw-unit printing of feature values (L239);
  - biplotEZ's role (L492);
  - the loan example uses 4-group CVA (L496);
  - the `bl_` prefix exceptions (L529);
  - the mgcv role (L572).
- **Structure:**
  - §2 file tree: the real files, including `development/`, `documentation/` and the new tests, and
    dropping `scripts/00_*` and `docs/implementation_summary.md`;
  - §3 hierarchy: add `bl_counterfactual`.
- **Status and dependencies:**
  - §7: the tests and vignette exist, `bl_pick_point` is documented, the commit column uses
    branches that exist, and rows are added for scaling and `bl_counterfactual`;
  - §8: per-package file usage, Suggests, and the base packages after S1.
- **Header (L4-5):** the branch and the last-updated date.
- **Also:** every row in the audit's "Low severity" table.

### S4: the walkthroughs, the harness script and one test comment

**Renumbering first.** In `review_section4_to_6.md`, map the walkthrough subsections onto the
script's own step headings:
- the integrated outlier filter becomes part of Step 4;
- 6a / 6b / 6c become 6.1 / 6.2 / 6.3, inside Step 6;
- the plot / project / overlay subsections "7 / 8 / 9" become "6.4 / 6.5 / 6.6", because the script
  runs them inside its Steps 4-6 block.

Then fix every cross-reference: `review_section7_to_8.md:61`, `review_section9_to_15.md:140, 142`,
and any others a grep finds.

**`review_section4_to_6.md`**, audit items 1-24:
- **Data sizes:** the dataset is 11 columns, and `loan_filtered` is 22,142 rows; the train/test
  counts after the hull come from the run log; the `feature_cols` example is the real 9 names.
- **Arguments:** `train_fraction` is in (0, 1]; the Step 6 call has `b_margin = 0.001` and the real
  title.
- **The custom model:** `model_type` is `"custom"` in `.pred_function` and in the `bl_results`
  fields; `bl_mod$model` is the list form.
- **Field tables:** the `bl_dat`, `bl_proj` and `bl_results` tables gain their missing fields.
- **The contour band:** `ct` is at 0.499 / 0.501, and `ct_surrogate` stores only those two levels.
  This includes the item-9 text corrected today.
- **The plot call:** `label_dir = "Hor"`, `label_offset_dist = 0.5`, `filter_to_polygon = TRUE`.
- **Stale references:** the removed `scripts/02_*` sentence goes; L590's script and vignette list
  becomes the real one; Step 5a is noted as commented out; RForest.

**`review_section7_to_8.md`**, items 1-12:
- the Step 7 call is `hist(bl_bnd$B_pred)`;
- the 0.499 / 0.501 levels and `b_margin = 0.001`;
- `"custom"`;
- the `B_pred range` format;
- `bl_pick_point(..., data = bl_results$test_data, bl_boundary = bl_bnd)`, describing test
  observations;
- present tense in place of "no longer", "new" and "legacy";
- the real `plot(bl_bnd)` console block, from the log;
- `distance` in the return list;
- the full `bl_bnd` flow fields.

**`review_section9_to_15.md`**, items 1-40:
- **The reduced model:** 4 features and the full model 9; p = 4, so 2^4 = 16 subsets and 6 pairs;
  `tdp = 13` everywhere.
- **Step numbering:** follows the script, so `bl_select_target` sits in Step 10,
  `bl_find_local_cf` in Step 11, and Phase 3 runs Steps 9-16.
- **Steps 9-10:**
  - the 9c call carries `plot_points = FALSE`, `ticks_v = 3` and `test_data[tdp:10, ]`;
  - the `print(tgt)` sample comes from the log, in raw units;
  - `set_filters` uses `loan_int_rate` / `credit_score` "increase" (no `person_age`), in the call,
    the `print(flt)` sample and the object flow.
- **Step 11:**
  - the "fixed" ±0.5 is in model units (±0.5 SD in v2);
  - `distance` joins the parameter table;
  - the Stage E filter order is actionability, then `train_ranges`;
  - the `bl_local` field list gains `dist_mahalanobis`, `distance`, `all_distances_mahalanobis` and
    `bl_counterfactual`;
  - the `print(bl_local)` sample comes from the log;
  - "update best if `dist_mah` improves" (the default selector).
- **Steps 13-14:**
  - the Shapley sign wording and the sort by signed value;
  - "Shapley table (raw units)" and "Variable summary (raw units)";
  - the raw-unit plot labels;
  - `round_to` acts in model units.
- **Wording:** the history phrases ("the former bug", "added in 2026-05-20", "(legacy)", "now
  available", "the new default", "Older objects") become present tense.
- **The global grid:** it has no 10% padding (L294).

**The harness script** (`development/03_loan_status_Boundary_Logic.R`), comments only, no code:
- L15: CVA, not PCA.
- L94-96: a model-free `bl_result`.
- L286 and L418: `bl_pick_point(..., data = ...$test_data, bl_boundary = ...)`.
- L632: pruned in "Prune least-important variables".
- L692-723: "four variables" / "top four". If the code there hard-codes 5, make it
  `length(feature_cols_v2)`, the one allowed code change, and say so.

**`tests/testthat/test-set_scaling.R:55`:** fix the comment, which says `bl_build_result()` returns
a `bl_projection`.

### S5: README, vignette prose, NEWS, CLAUDE.md

- **README.md:**
  - remove `bl_filter_outliers()` after `bl_prepare_data()`;
  - `rounding` becomes `b_margin`, in the examples and in the band note at :276-280;
  - split the model-types table into `bl_fit_model()` (GLM / SVM / NNET / RForest) and
    `bl_wrap_model()`-only (GAM / GBM / LDA / XGB as `list(model, features)`), with RForest;
  - make the Pima GAM example self-contained: read `pima_diabetes.csv`, `bl_prepare_data`, `gam`,
    `bl_wrap_model(model_type = "GAM")`, then Pima `set_filters`;
  - `bl_prepare_data` returns `bl_filter_result`;
  - `bl_robustness` returns a list;
  - add `bl_set_scaling()` and `bl_pick_point()` to the function tables;
  - add `build_vignettes = TRUE` to the install lines;
  - `docs/` is the pkgdown site.
- **Vignette prose** (no code changes, except the `eval=FALSE` chunk):
  - :133-135: a model-free `bl_result`;
  - :270-272: the hull clips the search; `train_ranges` only flags;
  - :473-474 and :489-490: the pair is chosen by squared Mahalanobis distance;
  - :278-284: `bl_pick_point(..., data = test_data)`.
- **NEWS.md:** correct the false claims at :5-8, :10-13, :39-43 and :107-109, and recheck :118-120
  after S2.
- **CLAUDE.md:**
  - §9 GAM note: `model_type = "GAM"` works directly (`"custom"` is optional);
  - §9 deferred tests: say which tests now exist (`test-local_cf.R` covers `bl_select_target`, and
    `test-result.R` covers `bl_assemble`), and keep the rest deferred.

### S6: verification and Gate 1c

- **Tests:** `devtools::document()` then `devtools::test()`. Expect 176 / 0 / 4.
- **Package check:** `devtools::check(args = "--no-manual")`. Expect no "Namespace dependencies
  missing" warning; report everything else against the earlier 0E / 0W / 2N.
- **Harness:** re-run it through `run_harness.R`. It must stay warning-free.
- **README examples:** extract the R chunks into a scratchpad script and run them, so every example
  runs. Display-only and `bl_pick_point` lines are skipped.
- **Greps:** `rounding`, `RForrest`, `stratified`, `bl_projection` (where it means a return value),
  `tdp = 2`, `person_age = "fixed"`, `6 features`, `0.49` / `0.51` and history words across the
  edited docs. Expect none to remain, except where correct.
- **Report:** the PASS / FAIL counts, the check result, the files changed, and any audit item that
  turned out not to be stale. Stop. Task C starts only when the user says so.
- **Not in scope:** the pkgdown site `docs/` will lag until a rebuild, which is a separate step.

---

## Task C: the Paper 1 Pima example

**Superseded (2026-11-02) by `.claude/reference/paper1-pima-example-plan.md`**, the approved Task C plan for
the current code. That plan has one shared `bl_results`, a Mahalanobis target rule, no
CVA-identity test, the restore comment and a run 6 (`"decrease"`). Its implementation log
records the first run. Where this section differs, that file wins.

All files go in `scripts/paper1-pima/`, with outputs in `scripts/paper1-pima/output/`. `scripts/` is
gitignored and buildignored, so this stays local. The set-up is exactly as the prompt fixes it: 392
rows, p = 7 with BloodPressure dropped in `feature_cols`, `train_fraction = 0.8`, `seed = 121L`,
`hull_fraction = 1`, CVA with the model (K = 3), and `max_pairs = 21L` (the default since 2026-10-31;
it was `10L`) with Mahalanobis distance.

**Design decisions within the prompt's spec:**
- **Loading and paths (`00_setup.R`):**
  - the single repo-specific line is `devtools::load_all("<pkg root>")`;
  - the output directory comes from `Sys.getenv("PAPER1_OUT", "output")`, which lets `tests.R` run
    the pipeline twice into separate folders;
  - the data path is `system.file("extdata", "pima_diabetes.csv", package = "boundarylogic")`, which
    works under both `load_all` and an installed package;
  - it records the package version, the git SHA and a working-tree "dirty" flag (or "installed").
- **`run_all.R`:** a small added driver that sources 01, 02, 03 and then the figure files in order.
  It is used by the reproducibility check.
- **GAM:**
  - fit `mgcv::gam(class ~ s(Pregnancies) + s(Glucose) + s(SkinThickness) + s(Insulin) + s(BMI) +
    s(DiabetesPedigreeFunction) + s(Age), family = binomial, data = train_data)`;
  - wrap it with `bl_wrap_model(model_type = "GAM")`;
  - in `01_data-model.R`, assert that its predictions equal `floor(predict.gam(type = "response") *
    1000) / 1000` on `test_data`. If they don't match, switch to `"custom"` with a `predict_fn`, and
    report which path was used.
- **Choosing the targets:**
  - use `bl_find_boundary(..., data = test_data)$dist_x`. This is the X-space distance on the
    `X_sd`-standardised features, and the only per-observation X-space distance the function
    returns; there is no Mahalanobis field;
  - drop NAs, then take the TP and the TN nearest their group's median `dist_x`;
  - state the field in the review page.
- **Reproducibility of each run:** `run_local()` calls `fresh_result()` every time. Numbers come only
  from that first call; nothing is repeated to "check" it.
- **Runs:** runs 1-5 for each target, exactly as in the prompt.
- **Figures:** 3-9 as `fig-03.R` to `fig-09.R`. Each writes a PDF and a PNG and plots the package's
  own objects only. `plot.bl_shapley()` returns a ggplot, so each call is wrapped in `print()`.
- **Figure 4(a):** the package has no non-interactive boundary-arrow drawer; arrows exist only in
  `bl_pick_point()`.
  - The first run draws the training points on the global biplot only, and logs an issue.
  - The proposed change must not add `boundary =` back to `plot_biplotEZ()` (CLAUDE.md §4). The
    proposal is a `type = "biplot"` (arrows) option in `plot.bl_boundary()`, or a non-interactive
    `rows =` mode for `bl_pick_point()`.
- **`tests.R`:** all the prompt's assertions.
  - The reproducibility check runs `Rscript run_all.R` twice through `system2()`, each run in a fresh
    process with its own `PAPER1_OUT`.
  - It then compares the two `numbers.rds` files, ignoring the version fields.
- **`review.qmd` to `review.html`:** shows the figures, Table 1, a pass/fail line for each
  assertion, and the issues table (symptom, where, effect on the paper, proposed fix, decision left
  blank).

**Issues found while planning.** These are pre-seeded into the issues table so they get a fix-now /
defer decision:

| issue | where |
|---|---|
| A reverted "fixed" variable can still show `used_in_sparse = TRUE` | `R/shapley.R:433` |
| The sparse-CF cross ignores `rotate_deg`; its colours are green3/yellow2, not the documented blue/red | `R/shapley.R:490-499` (D items 8 and 9). **Fixed 2026-10-30:** `rotate_deg` removed from the local plots |
| `cutoff` is not passed to `bl_build_grid()` | `R/result.R:341-350` (D item 7). **Fixed 2026-10-30** |
| An empty outcome group has no guard: with 3 levels `plot.bl_projection` silently colours every point blue | `R/projection.R:253-289`. **Fixed 2026-10-30** |
| `bl_shapley(seed =)` has no effect when p <= 14, which includes Pima; it is a doc point | |
| The implementation summary is stale: "stratified" at :103 (the split is simple random), "RForrest" at :124, and `no_points` / `show_arrows` / `arrow_col` at :394 | `2 implementation_summary.txt` (D item 9) |
| Figure 4(a) arrows, as above | |

**Fixed before Task C (2026-10-30).**
- **The three rows marked fixed above**, plus `xgboost::xgboost()` models failing in
  `.pred_function()`, were fixed ahead of Task C's first run.
  - The plan is `.claude/reference/gate2-code-fixes-plan.md`, archived as
    `.claude/reference/gate2-code-fixes-plan.md`.
  - Task D records prompt items 7 and 8 as fixed, with a pointer, the same way as item 9. They are
    not deferred.
- **The `used_in_sparse` row** was fixed in session 2 (2026-10-29).
- **The random nearest-vertex tie-break** (prompt items 1-2) was fixed the same day; see
  `.claude/reference/nearest-vertex-tie-break-plan.md`.
- **Effect on Task C:** the reproducibility check no longer depends on fresh sessions. A repeated
  call in one session must now give identical results too.

### Gate 2
- Report the path of `review.html`, the assertion results and the issues table. Stop.
- For each issue the user marks "fix now":
  - a short plan, the fix, a test and a NEWS entry;
  - the walkthrough and harness-script updates, as described above;
  - `devtools::test()`;
  - then re-run Task C from a fresh session.
- Each issue marked "defer" goes to Task D.

---

## Task D: record the deferred issues (no code changes)

- **The note:** `documentation/known-issues.md`, with one entry per issue: symptom, file:line, impact
  on results, and suggested fix. It covers prompt items 1-9, the planning findings above, and every
  issue deferred at Gate 2. Prompt item 9 (stale documentation) is fixed by the sweep (S1-S6), so
  it is recorded as fixed, with a pointer, rather than as deferred. It also records the two items
  the user parked on 2026-10-29: the README Pima constraint example returning
  `solution_found = FALSE`, and the harness plot title "(prior defaulters)".
- **Prompt items 1 and 2 are fixed (2026-10-30), not deferred.** Both were the random
  nearest-vertex tie-break in `.nearest_idx_block()`, not the rotation: `best_pair` instability
  across identical calls, and dependence on earlier calls. Record each as fixed, with a pointer to
  `documentation/nearest-vertex-tie-break.md`, the same way as item 9. The rotation's own plane
  under-determination (CLAUDE.md §9 and rotation note §4.1) is a separate, still-open issue: cite it
  rather than repeating it.
- **Corrected 2026-11-04 (Task E):** the next bullet is too strong. The random `max.col()` is only in
  the branch's post-thesis optimized script. The thesis text calls the method deterministic (p.92),
  and Appendix A uses `which.min`.
- **Task E must note:** the thesis's searches used the same random tie-break, inherited from the
  PhD code (`1.4 Biplot Boundary Search_optimized.R:45`). That note's Section 5 gives the scale.
- **CLAUDE.md §9:** one new bullet pointing to the note.
- **Done 2026-11-04** (plan `.claude/reference/known-issues-taskD-plan.md`):
  `documentation/known-issues.md` has 9 open entries (prompt items 3-6, the deferred package
  options of Pima issues 12 and 13, the two parked items) and 8 resolved ones (prompt items 1, 2, 7,
  8 and 9, and the planning findings: `xgboost()` models, empty CVA groups, `used_in_sparse`).
  The arrow layer is left out, by your choice. CLAUDE.md §9 has the pointer bullet.

## Task E: methodology changes since the thesis (for the supervisor)

- **Deliverable:** `documentation/methodology-changes/methodology-changes.qmd`, rendered to `.docx`
  next to it. It records the package commit it was measured at.
- **Prose:** written with the `academic-writing-style` skill, because it is thesis-derived prose by
  the author.
- **Baselines:**
  - the thesis text is the final 27 Feb 2026 version, through the `thesis-md` chapters 07-11, not
    the 30 Sept 2025 PDF in `documentation/`;
  - the thesis code is read only with `git show` / `git diff origin/original_PhD_code`, never
    checked out;
  - where the text and code disagree, report both; the text wins.
- **Part 1: the changes made.** Verify the prompt's nine items against the sources, correct them,
  and add any others.
  - **Already verified:** W = pooled within-group scatter / (n - g) with a ridge fallback
    (`R/projection.R:43-52`); the ±0.5 "fixed" tolerance (`R/local_cf.R:89`); the canonical pair
    order with `max_pairs = 10L` (the default until 2026-10-31; now `21L`).
  - **Candidates to check and add:** global CFs outside `train_ranges` are now kept rather than
    discarded; outlier filtering by `hull_fraction` (default 0.9) inside `bl_prepare_data()`; grid
    probabilities clamped to [0, 1]. Exclude anything that only affects speed or the API, such as
    the filter order and the `kappa(V)` check.
  - **The pair (1,7) claim:** check it against the thesis.
- **Measurements:** `scripts/paper1-pima/04_method-change-effects.R`, using the Task C set-up with a
  fresh `bl_results` for each run. It writes `documentation/methodology-changes/method-change-effects.csv`.
  - **Change 1:** Mahalanobis against Euclidean importance, as the per-variable ranking and median
    from `plot.bl_boundary()` / `bl_robustness()`.
  - **The pair selector:** Mahalanobis against Euclidean, for both targets: `best_pair`, the CF,
    √d²_M and `solution_found`.
    - **Flag (2026-10-30, open for decision):** under CVA the two selectors are identical by
      construction. Every rotated basis is W-orthonormal, so d²_M = (n−g)·dist_z²; the loan harness
      measures a ratio of 15992 for every pair (rotation note §3.3).
    - **Decide:** measure the comparison under PCA, or report the identity, depending on Task C's
      projection.
  - **Change 3:** `max_pairs = 21L` against `10L`, for both targets.
  - Every other change is described and marked "not measured".
- **For each change:** the thesis version (with section), the package version (with function), the
  reason, the effect on Pima, and a decision table (approve / query / reject, plus a comment cell).
- **Appendix:** one line on the bug fixes from A, B and C.
- **Part 2: future changes.** The prompt's ten items, each with what it is, why, the expected effect,
  the effort and a recommendation.

### Gate 3
- Report the drafts of D and E and the path to the `.docx`. Stop.

---

## Error cases

| case | action |
|---|---|
| Part 0: the push is rejected (non-fast-forward) | Stop and report. No force push, no rebase without asking. |
| Part 0: something unexpected is staged, or a file is missing | Stop before committing and show `git status`. |
| A line reference has moved | Find the function by name, and say so. |
| The prompt conflicts with CLAUDE.md | CLAUDE.md wins; say where. Known: the branch, Rscript through PowerShell. |
| A device test misbehaves without a screen | Use `pdf(NULL)` only. |
| An outcome group is empty, so CVA fails | Stop with a message; do not change the set-up. |
| A constrained run returns `solution_found = FALSE` | Report it in Table 1; do not re-pick the target. |
| `best_pair` differs between two fresh runs | Log it as an issue; never re-run to hide it. |
| `s(Pregnancies)` fails | Lower its `k` and flag it as a set-up change for approval. Not expected: there are 17 distinct values. |
| The GAM path's predictions mismatch | Use `"custom"` + `predict_fn`, and report it. |
| A Task E change cannot be measured | Describe it and say "not measured". |
| The thesis text and code disagree | Report both; the text defines the thesis version. |

## Verification

- **Part 0:**
  - `git log -1 --stat` shows the documentation-only commit;
  - `git status -sb` shows `method_developments` level with origin;
  - `git branch -vv` shows `paper_nr1_changes` tracking `origin/paper_nr1_changes` at the same SHA.
- **Gate 1:**
  - `devtools::document()` gives a clean NAMESPACE diff (adds `dev.set`);
  - `devtools::test()` passes, with the new `test-graphics_device.R` and `test-local_cf.R` included;
  - `03_loan_status_Boundary_Logic.R` runs end to end from the scratchpad driver;
  - the three walkthroughs are updated for A and B, in present tense.
- **Gate 2:**
  - `tests.R` passes or reports each failure;
  - the reproducibility check passes across two fresh processes;
  - `review.html` renders.
- **Gate 3:** `methodology-changes.docx` renders, and every "measured" claim traces to a CSV row.
- **Finally:**
  - run `devtools::test()` again after any Gate 2 fixes;
  - update `progress.md` with the files changed, decisions, test status and next step;
  - copy this plan to `.claude/reference/` with an IMPLEMENTED banner;
  - wait for instructions on commits, a merge into `method_developments` / `main`, and a release
    tag for `install_github("AI-Rowan/boundary_logic@<tag>")`.
