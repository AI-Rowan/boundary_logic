# Progress

## Session summary (2026-05-13)

### Completed this session

1. **Created three code walkthrough documents** in `.claude/reference/` — deep function-by-function reference for the loan-default script (`scripts/03_loan_status_Boundary_Logic.R`).

   | Document | Steps covered | Functions documented |
   |---|---|---|
   | `review_section4_to_6.md` | Steps 3–6 (Phase 1) | `bl_prepare_data`, `bl_filter_outliers`, `bl_fit_model`, `bl_build_result` (+ `bl_build_projection`, `bl_build_grid`, `bl_assemble`), `plot_biplotEZ`, `bl_project_points` |
   | `review_section7_to_8.md` | Steps 7–8 (Phase 2) | `bl_find_boundary`, `print.bl_boundary`, `plot.bl_boundary` (jitter + boxplot), `bl_robustness`, `plot_biplotEZ` boundary overlay |
   | `review_section11_to_17.md` | Steps 11–17 (Phase 3) | `bl_predict`, `bl_select_target`, `set_filters`, `bl_find_local_cf` (+ `.bl_rotate` SVD rotation), `plot.bl_local_result`, `bl_shapley` (exact + approximate), `plot.bl_shapley`, `bl_find_sparse_cf`, `print/plot.bl_sparse_result` |

2. **Set up reference document infrastructure:**
   - Created `.claude/reference/` folder and moved review docs there
   - Updated `CLAUDE.md` Section 2 with a two-tier reference structure (generic + script-specific) and per-document update triggers
   - Added `reference_docs.md` to memory and updated `MEMORY.md` index
   - Updated `.gitignore` to track `.claude/reference/` while keeping `.claude/plans/` and `.claude/settings` ignored

3. **Committed all outstanding work** — commit `bf025e7` on `method_developments`:
   - Phase 3 source files (`local_cf.R`, `shapley.R`)
   - Source fixes (`model_utils.R`, `plot_biplot.R`, `surrogate.R`)
   - Test suite fixes (72 PASS, 0 FAIL)
   - Full `man/` regeneration (14 updated, 14 new `.Rd` files)
   - Loan dataset (`inst/extdata/loan_data.csv`, `R/datasets.R`)
   - Project docs (`CLAUDE.md`, `progress.md`, reference docs, loan workflow script, design reference txt files)

4. **Fixed script bug** — restored `points = test_pts_v2` in Step 11 of `scripts/03_loan_status_Boundary_Logic.R` (had been changed to `NULL`, causing implicit `bl_project_points` call that triggered a debug browser).

---

### Dead ends this session

- **Parallel agent spawning failed** — first attempt to launch three Explore agents for the Step 11–17 review hit API connection errors (`ConnectionRefused`, `FailedToOpenSocket`) on all three simultaneously. Fell back to reading source files directly with the Read tool.
- **`.gitignore` exception pattern failed** — initial attempt to use `!.claude/reference/` to un-ignore the reference folder didn't work because git cannot un-ignore files inside an ignored parent directory. Fixed by restructuring `.gitignore` to ignore specific subdirectories (`.claude/plans/`, `.claude/settings*.json`) rather than the whole `.claude/` folder.

---

### Current state

- Branch: `method_developments` — 1 commit ahead of `origin/method_developments`, not yet pushed
- Working tree: clean
- Test suite: **72 PASS, 0 FAIL**

---

### Next steps

**Immediate fixes (pick one to start):**

1. **Fix stale CLAUDE.md Section 6 note** — correct "`bl_build_result()` returns a `bl_projection`" to reflect actual behaviour (always returns `bl_result`).

2. **Accuracy fixes (recommended next implementation):** implement `.claude/plans/accuracy-correctness-fixes.md` — 5 purely defensive fixes, no method changes, highest scientific value:
   - Condition number check before `solve(V)`
   - Grid probability clamping to [0, 1]
   - Boundary pruning fallback warning
   - Surrogate hull coverage context in accuracy reporting
   - Shapley seed independence (change default from `1L` to `NULL`)

3. **Usability bug fixes:** implement `.claude/plans/usability-bug-fixes.md` — 7 fixes, several outright crashes:
   - NULL model crash in `bl_assemble`
   - `bl_build_result` silent NULL return
   - `print.bl_points` crash when no model
   - RForrest→RForest typo + backward-compat shim
   - Silent ignored hyperparameters
   - Non-0.5 cutoff warning
   - Rounding doc mismatch (4 d.p. → 3 d.p.)

4. **Documentation gaps:** implement `.claude/plans/documentation-gaps.md` — `devtools::document()` plus stale roxygen fixes.

**Deferred (do not implement without instruction):** see CLAUDE.md Section 9 for the full list (Phase 3 unit tests, zero-length arrow warning, GAM fix, vignette, CRAN prep, Mahalanobis alternative).

### Verification command
```r
"/c/Program Files/R/R-4.6.0/bin/Rscript" -e "devtools::test()"
```

---

## Session summary (2026-05-08)

### Completed this session

1. **Fixed 4 pre-existing test failures** — all test-side only, no source changes.

   | File | Fix applied |
   |---|---|
   | `tests/testthat/test-predict_utils.R` | Removed stale `2L` rounding arg from all `.pred_function()` calls; updated `* 100` integer check to `* 1000` (predictions are 3 d.p., not 2); fixed error test to use a non-workflow model object |
   | `tests/testthat/test-model_fit.R` | Removed `result$rounding` from `.pred_function()` call; removed `rounding = 2L` from `bl_fit_model()` call; updated `* 100` → `* 1000` assertion |
   | `tests/testthat/test-projection.R` | Updated `expect_error()` regexp to `"method = 'CVA' requires"` |
   | `tests/testthat/test-result.R` | Removed extra `NULL` positional arg from all `bl_assemble()` calls; updated `hull_fraction` test to check it's a numeric in [0,1] from `bl_grid` (not NULL) |

   **Final test result: 72 PASS, 0 FAIL, 4 WARN** (warnings are pre-existing CVA/iris 2-class notices from biplotEZ — not failures).

2. **Comprehensive code review** — full audit of all R/ files across three dimensions: accuracy/correctness, usability, and documentation.

3. **Created 3 improvement plan files** in `.claude/plans/`:

   | Plan file | Contents |
   |---|---|
   | `accuracy-correctness-fixes.md` | 5 fixes: condition number check for `solve(V)`, grid prob clamping [0,1], boundary pruning warning, surrogate hull coverage reporting, Shapley seed independence |
   | `usability-bug-fixes.md` | 7 fixes: NULL model crash in `bl_assemble`, `bl_build_result` silent NULL, `print.bl_points` crash, RForrest→RForest typo+shim, silent hyperparams, cutoff warning, rounding doc mismatch |
   | `documentation-gaps.md` | Run `devtools::document()` to generate 11 missing `.Rd` files, fix stale roxygen (`rounding` param, 4 d.p. → 3 d.p.), remove orphaned rlang imports |

4. **Created `.claude/settings.json`** — set `"plansDirectory": ".claude/plans"` so future plans save to the project directory.

---

### Current state

- Branch: `method_developments`
- Uncommitted changes: `R/plot_biplot.R`, `R/local_cf.R`, `R/surrogate.R` (the `par(new=TRUE)` fix from last session), plus 4 test files updated this session
- Test suite: **72 PASS, 0 FAIL**

---

### Known issues

- **Zero-length arrow warning** — `plot_biplotEZ()` with `boundary =` generates "zero-length arrow is of indeterminate angle and so skipped" warnings from R's graphics device. Three suppression approaches were tried and all failed (pre-filter, `suppressWarnings()`, `withCallingHandlers()`). The warning fires at deferred device-flush time, outside the R call stack. All attempted fix code has been reverted. Recommended next step: change `show_arrows` default to `FALSE` in `plot_biplotEZ()` — at large observation counts the boundary crosses suffice without the arrows. See CLAUDE.md Section 9 for full details.

---

### Next steps (pick any plan to implement)

1. **Accuracy first (recommended):** implement `.claude/plans/accuracy-correctness-fixes.md` — purely defensive, no method changes, highest scientific value.
2. **Quick wins:** implement `.claude/plans/usability-bug-fixes.md` — several are outright crashes.
3. **Documentation:** implement `.claude/plans/documentation-gaps.md` — mostly running `devtools::document()` plus small source fixes.

After completing any plan, commit the changes on `method_developments` and update this file.

### Verification command
```r
"/c/Program Files/R/R-4.6.0/bin/Rscript" -e "devtools::test()"
```
