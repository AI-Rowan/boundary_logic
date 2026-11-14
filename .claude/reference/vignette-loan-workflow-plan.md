> **IMPLEMENTED 2026-11-06** (committed in `17d3c18`). The loan-default vignette rewritten for
> package users and precomputed (`vignettes/precompile.R`).

# Plan: rewrite the loan-default vignette for package users, in line with the harness

Date created: 2026-11-06
Date last edited: 2026-11-06

> **Housekeeping (done 2026-11-06):** this plan was written in plan mode into
> `before-updating-the-paper-wiggly-wave.md`, which holds the algorithms-document plan (waiting at its
> review gate). Step 0 moved it here and restored that file.

## Implementation notes (2026-11-06)

Steps 0–9 are done:
- `devtools::test()`: 527 PASS, 0 FAIL, 4 WARN;
- `devtools::check()`: 0 errors, 0 warnings, 0 notes;
- pkgdown rebuilt. pkgdown copies `vignettes/figures/` to `docs/articles/figures/` (24 PNGs, none missing).
  It also brought 7 stale reference pages up to date and added `dot-grid_in_hull.html`. Its site report
  flags `url` missing in `_pkgdown.yml`, which was already the case.

The harness ran in 60 s, and its key results match the knitted vignette line for line (step 5). The knit takes about 63 s and produces 24 figures (about 0.6 MB). The rendered HTML is 1 MB, with
every image embedded.

**Deviations from the harness found during implementation** (display only; no result changes), in addition
to the four planned:

5. **Biplot figures are square:**
   - 7 × 7 in for the four-feature plots;
   - 10 × 7 in (`biplot_wide`) for the nine-feature plots, whose long horizontal labels were clipped.
6. **Label offsets on the nine-feature CVA biplots** (Step 6 and its test overlay, Step 8):
   - `loan_amnt` and `credit_score` are lifted by one line (they overlapped `loan_int_rate` and
     `person_income`);
   - on the Paral plots, `cb_person_cred_hist_length` and `person_age` are moved out by 1 and 2 lines;
   - Step 6's `bl_build_result()` title ends in `"\n"`, so it clears the lifted labels (commented in the
     code).
7. **Step 16:**
   - `plot(bl_sparse_ext)` is left out;
   - the text explains that the new applicant (credit score 420, below the training minimum) lies outside
     the local plot's window, which is set by the training data (`algorithms.md` F4), so neither the
     applicant nor its arrow shows;
   - the text says the sparse counterfactual keeps the out-of-range score of 420, because
     `bl_find_sparse_cf()` checks only the crossing.
8. **Extra prints for learners:**
   - `print(bl_bnd)` in Step 7;
   - `print(tgt)` in Step 10;
   - `head(pred_summary)` and `pred_summary[tdp, ]` in Step 9;
   - Step 7b's repeat of the reduced biplot is not drawn again.
9. **Hidden `*-claims` chunks** (`stopifnot`) guard each qualitative statement in the text, as the
   methodology note's `claims` chunk does.

## Context

`vignettes/Boundary_Logic_loan_default_workflow.Rmd` has drifted from the development harness
`development/03_loan_status_Boundary_Logic.R`. The main differences:
- it draws a 4,000-row subsample;
- it uses the built-in XGB wrapper for both models;
- it keeps `person_age` and drops four variables, not five;
- it does not standardise the reduced model or call `bl_set_scaling()`;
- it picks the first false negative as its target, not test row 13;
- it fixes `person_age`, where the harness fixes `loan_amnt` and uses `c("increase", "integer")`;
- its external applicant skips `bl_apply_scaling()`;
- it calls `bl_robustness()` as a step, which CLAUDE.md §4 and §8 forbid.

The user wants the vignette to follow the harness's code and data, and to be written for an external
audience who need to learn how to apply each function and its parameters. It must still end with the
comparison with SHAP.

**User decisions (2026-11-06):**

| Topic | Decision |
|---|---|
| Data | The full data, as in the harness, **precomputed**: the source is `vignettes/Boundary_Logic_loan_default_workflow.Rmd.orig`, knitted locally into the shipped `.Rmd`, with static output and figures |
| SHAP feature values | Raw loan units, in the vignette only: shapviz's `X` (display only) gets the values converted back from z-scores. The harness is not changed. |
| pkgdown | Rebuild `docs/` afterwards, so the published article matches |

## Files to read first

- `development/03_loan_status_Boundary_Logic.R`: the source of every call, parameter, target and constraint.
- `vignettes/Boundary_Logic_loan_default_workflow.Rmd`: the current text, especially the SHAP section, which is
  kept and adapted.
- The roxygen headers of every function the vignette introduces, for the defaults in the parameter tables:
  `bl_prepare_data()`, `bl_wrap_model()`, `bl_build_result()`, `plot_biplotEZ()`/`plot.bl_result()`,
  `bl_project_points()`, `bl_find_boundary()`, `plot.bl_boundary()`, `bl_set_scaling()`, `bl_surrogate()`,
  `bl_predict()`, `bl_select_target()`, `set_filters()`, `bl_find_local_cf()`, `plot.bl_local_result()`,
  `bl_shapley()`, `bl_find_sparse_cf()`, `bl_apply_scaling()`, `bl_pick_point()`.
- `development/reference/review_section4_to_6.md`, `review_section7_to_8.md` and `review_section9_to_15.md`:
  the walkthroughs. Use them for the reasons behind the harness's choices (for example, why `person_income`
  is dropped), not for their numbers.
- `documentation/scaling_axis_relabel_note.md` ("Display conventions"), CLAUDE.md §4, §5 and §8 (what a
  vignette must not show).

## Design

### What follows the harness exactly

- The data (all 45,000 rows, then the no-prior-default filter), the seed 121 and the 0.8 split.
- `hull_fraction = 0.9` in the model fits and 1 in the exploratory biplot.
- The full XGBoost model (9 features), wrapped with `model_type = "custom"` and a `predict_fn`.
- The five dropped variables and the reasons the harness gives.
- The reduced model:
  - standardised on the training split only, with `bl_set_scaling()`;
  - wrapped with `model_type = "XGB"` and `list(model, features)`;
  - built with `b_margin = 0.01`.
- The target: test row 13, a true negative with a predicted probability of 0.111. Its constraints:
  - `loan_amnt = "fixed"`;
  - `loan_int_rate = "increase"`;
  - `credit_score = c("increase", "integer")`;
  - `hull = TRUE` and `train_ranges = TRUE`, set explicitly.
- The unconstrained alternative, and the external applicant converted with `bl_apply_scaling()`.
- The SHAP section: the full ranking set beside the distance ranking, global SHAP on the reduced model
  (bar and beeswarm), and local SHAP for row 13 (waterfall and force).
- **Step numbers follow the harness** (1–3, 4–6, 7, the pruning step, 4b–6b, 7b, 8, 8b, 9–16), so the two can be
  compared line for line.

### Deliberate differences from the harness

1. **No `bl_robustness()` step** (CLAUDE.md §4, §8). Step 7 keeps the invisible return value,
   `bnd_imp <- plot(bl_bnd)`, and the pruning step reads `bnd_imp$sum_of_distance`. This returns the same
   numbers without a redundant call.
2. **SHAP feature values in raw units:** convert with `bl_results_v2$scaling$center` and `$scale`
   (`x * scale + center`) and pass the result to `shapviz(X = ...)`. The SHAP values are unchanged.
3. **The harness's working lines are left out:**
   - `hist(bl_bnd$B_pred)`;
   - the stray `bl_results$test_data`;
   - the `tdp:10` point subset drawn with `plot_points = FALSE` (Step 9 draws only the target).
4. **Interactive calls** (`bl_pick_point()`) are shown with `eval = FALSE`, as CLAUDE.md §8 requires.

### Writing for package users

- **Each step:**
  - what it does and why, in one or two sentences;
  - the code;
  - a **Parameters** table (`Parameter | Default | This call | Effect`) for each function the step
    introduces, listing only the arguments set here and the defaults that change the result, with a pointer
    to `?fn` for the rest;
  - how to read the output.
- **Notes users need:**
  - the `c(min, max)` constraint is in model units (z-scores for the reduced model);
  - `bl_set_scaling()` changes only the display;
  - new raw-unit data goes in through `bl_apply_scaling()`;
  - `cutoff` stays at 0.5;
  - the two ways to wrap a model (`"custom"` with `predict_fn`, or a built-in type), and `bl_wrap_data()`
    for data that is already split (prose only).
- **A short data dictionary in Step 1:** each kept column, its meaning and its units.
- **Quoted numbers:** written as inline `r` expressions in the `.Rmd.orig`, so they become static text from
  the same run and cannot drift. The only typed numbers are settings.
- **Things not to show** (CLAUDE.md §8): `outlie`, `calc_hull`, `boundary =` in `plot_biplotEZ()`, and
  `bl_filter_outliers()` as a step (it may be mentioned in prose).
- **Style:** plain technical English, present tense. Keep the ending "Summary" table (step | function |
  purpose), updated.

### The precomputed build

| File | Role |
|---|---|
| `vignettes/Boundary_Logic_loan_default_workflow.Rmd.orig` | The source that is edited: real chunks, with the SHAP chunks still guarded by `has_shap` |
| `vignettes/precompile.R` | Loads the development package with `pkgload::load_all()` and stops unless `pkgload::is_dev_package("boundarylogic")` is true, so it never knits against a stale installed copy. It runs `knitr::knit()` from `vignettes/` with `fig.path = "figures/loan-"`, then confirms that every referenced figure exists |
| `vignettes/Boundary_Logic_loan_default_workflow.Rmd` | Generated, and shipped. It keeps the same YAML and `VignetteIndexEntry`, so the article URL stays the same. A comment at the top says it is generated and must not be edited |
| `vignettes/figures/loan-*.png` | Generated figures. `html_vignette` embeds them when the vignette is built |
| `.Rbuildignore` | Add `^vignettes/.*\.Rmd\.orig$` and `^vignettes/precompile\.R$` |

## Steps

0. **Housekeeping:**
   - restore the algorithms plan to `.claude/plans/before-updating-the-paper-wiggly-wave.md`, word for word;
   - write this plan to the plan-mode working copy (removed 2026-11-07).
1. **Get the reference numbers:** run the harness end to end with a scratchpad driver, as CLAUDE.md
   describes. Its log is the reference for the vignette's numbers.
2. **Check every default:** read the signature and roxygen of each function in the list above.
3. **Write the `.Rmd.orig`:**
   - Overview;
   - Setup;
   - Phase 1 (Steps 1–6);
   - Phase 2 (Step 7 and pruning);
   - the refit (Steps 4b–6b);
   - Phase 2 again (Steps 7b, 8, 8b);
   - Phase 3 (Steps 9–16);
   - Comparison with SHAP;
   - Summary.
4. **Write `vignettes/precompile.R`, update `.Rbuildignore` and run the precompile.** Then read the knitted
   `.Rmd`:
   - every quoted number comes from the run;
   - every figure is legible. Use `ticks_v`, `ticks_var` and `ticks_n` where an axis lacks ticks
     (CLAUDE.md, "Biplot ticks in figures"), and check each PNG by eye.
5. **Compare with the harness log from step 1.** These must be identical:
   - the predicted probability for row 13;
   - the variable ranking;
   - the local and sparse counterfactual tables, including `credit_score` 634.22 rounded to 635;
   - the surrogate agreement figures;
   - the external applicant's counterfactual.
6. **Render the shipped `.Rmd`** with `rmarkdown::render()` into `tempdir()` (the Windows path limit) and
   read the HTML.
7. **Update the related files:**
   - **NEWS.md:** one bullet under 0.2.0, "Bug fixes and other improvements", saying the vignette now follows
     the harness on the full data, is written for package users, and is precomputed;
   - **README:** refresh the one-line vignette description;
   - **CLAUDE.md §5:** one line saying the vignette is precomputed (edit the `.Rmd.orig`, run
     `vignettes/precompile.R`, never edit the knitted `.Rmd`);
   - **`2 implementation_summary.txt`:** check its vignette mentions still hold.
8. **Run the checks:**
   - `devtools::test()`, a sanity check, since no R code changes;
   - `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")` with `RSTUDIO_PANDOC`
     set. The target is 0 errors, 0 warnings and 0 notes. If a `.docx` is open in Word, use the CLAUDE.md
     work-around.
9. **Rebuild pkgdown**, following the CLAUDE.md procedure:
   - a temporary `.R` script with `RSTUDIO_PANDOC`, `clean_site(force = TRUE)` and `build_site()`;
   - `unlink()` the generated `CLAUDE`/`progress` pages;
   - delete the helper script;
   - check that `docs/articles/` holds the new article and its figures.
10. **Stop and report:**
    - the test and check counts;
    - the changed files;
    - the deliberate differences from the harness.

    Do not commit. Remind the user about `progress.md` and `/revise-claude-md`.

## Error cases and risks

- **Knitting against a stale installed package:** `precompile.R` stops unless the development copy is loaded.
- **The harness and the vignette disagree on a number:** step 5 compares them. Any difference is a bug in the
  `.Rmd.orig`; the harness is the reference.
- **Figures missing from the pkgdown article:** step 9 checks `docs/articles/`. If pkgdown does not copy
  `vignettes/figures/`, set `fig.path` so the figures sit next to the `.Rmd` and re-knit.
- **The Windows path limit:** precompile inside `vignettes/`, which has short relative paths, and render
  checks into `tempdir()`.
- **The shipped `.Rmd` goes stale:** it no longer runs code at build time, so R code changes do not reach it.
  The CLAUDE.md line (step 7) and the comment at the top of the generated file say to re-run `precompile.R`
  after code changes.
- **Long precompile:** SHAP on about 16,000 training rows with `nsim = 50` may take minutes. This is run
  only by the maintainer, so it is acceptable.
- **`fastshap` or `shapviz` missing when precompiling:** both are installed here. `has_shap` stays as the
  guard for other maintainers.

## Verification

- The precompile runs without errors. The knitted `.Rmd` contains no live chunks and refers only to
  figures that exist.
- The key numbers match the harness log exactly (step 5).
- The rendered HTML reads correctly and every figure is legible.
- `devtools::test()` passes. `devtools::check()` gives 0 errors, 0 warnings and 0 notes.
- The pkgdown article shows the new vignette with its figures, and no `CLAUDE`/`progress` page is published.
- `git status` shows only these files:
  - the vignette files;
  - `.Rbuildignore`, `NEWS.md`, `README.md` and `CLAUDE.md`;
  - `docs/`;
  - the earlier uncommitted files from this session.
