# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Standing items / carry-overs

- **Paper 1 work is on branch `paper_nr1_changes`** (branched from `method_developments` at
  `b7aee65` on 2026-10-29, pushed with upstream). Start Paper 1 sessions on that branch
  (CLAUDE.md §10). The task list is `.claude/plans/boundary-logic-fix-prompt.md` (Tasks A–E,
  three review gates); the working plan, with an implementation log, is
  `.claude/plans/commit-and-push-the-mutable-hearth.md`. Both are gitignored, local only.
- **`main` was last merged and pushed on 2026-08-03.** `method_developments` is at `b7aee65`
  (pushed 2026-10-29). Merging `paper_nr1_changes` back, and the release tag the paper will pin,
  wait for the user's instruction after Gate 3.
- **`paper_nr1_changes` state:** `f2a076b` (session 1) and `7aa1344` (session 2) are committed
  but **not pushed**. Session 3's Gate 2 code fixes (below) are **uncommitted** in the working
  tree. Commit only when the user says so.
- **`scripts/` is gitignored** — the user's personal scratch playground, never published. Task C
  will write `scripts/paper1-pima/` there (local only, by design).
- **`.bl_rotate()` plane under-determination** — documented, deliberately not actioned.
  See `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1 and §6
  (a proposal, not implemented), and CLAUDE.md §9. Do not "fix" without an explicit instruction.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 278 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are pre-existing
  biplotEZ CVA 2-class notices). **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-10-30,
  session 3).
- **Verification:** run R as `.R` files via the PowerShell tool,
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <driver>.R`, with the driver calling
  `devtools::document(<pkg root>)` then `devtools::test(<pkg root>)`. The dev harness runs end to
  end through a scratchpad driver (CLAUDE.md §5).
- **`.claude/plans/` is gitignored** — plan files are local only. `.claude/reference/` *is*
  tracked (generic docs + archived plans).

---

## Session summary (2026-10-30, session 3, latest) — Gate 2 code issues fixed before Task C

Goal: the user asked for the background and fixes of the four code issues pre-seeded in Task C's
Gate 2 issues table, then approved fixing all four now. The plan is
`.claude/plans/can-you-summarise-the-replicated-magpie.md`, archived as
`.claude/reference/gate2-code-fixes-plan.md`. Session 2's work was committed first as `7aa1344`
(user's instruction; not pushed).

### What was done

1. **`rotate_deg` removed from the local plots** (`R/local_cf.R`, `R/shapley.R`).
   - The user's rule: `rotate_deg` is only for viewing the global biplot. Whenever `.bl_rotate()`
     is used, it is in effect 0, and it never turns a rotated individual point.
   - `plot.bl_local_result()` lost the argument and the rotation block. Passing `rotate_deg` to it
     or to `plot.bl_sparse_result()` now stops with an error before any device opens.
   - This removes the sparse-cross bug. `plot_biplotEZ()` and `plot.bl_surrogate()` keep
     `rotate_deg`.
2. **Grid cutoff** (`R/biplot_grid.R`, `R/result.R`).
   - `bl_build_grid()` defaults to `cutoff = bl_model$cutoff` and stores `bl_grid$cutoff`.
   - `bl_assemble()` stops if the grid's and the model's cutoffs differ. A grid with no `cutoff`
     field (an older object) skips the check.
3. **`xgboost::xgboost()` models** (`R/predict_utils.R`, `R/model_fit.R`).
   - The XGB branch scores `"xgboost"`-class models from the data frame. Their predict method
     refuses an `xgb.DMatrix`.
   - `bl_wrap_model()` stops unless such a model is `binary:logistic` on levels 0/1. Fitted with a
     numeric response, `xgboost()` gives a regression.
4. **Empty CVA group** (`R/projection.R`).
   - A message names an empty TP/TN/FP/FN group, and empty levels are dropped before CVA and `W`.
   - `point_col` uses the confusion palette for any subset of the four labels.
5. **Documents:**
   - NEWS (one breaking change, three fixes) and roxygen;
   - `2 implementation_summary.txt` (§4.3 `bl_grid$cutoff`, §4.4 and §4.8 the `rotate_deg` note,
     the XGB passages, CVA empty groups);
   - walkthroughs `review_section4_to_6.md` and `review_section9_to_15.md` (the `rotate_deg` note
     after the `plot(bl_local)` Parameters table);
   - README, and the rotation note §4.6;
   - CLAUDE.md §4 ("never add `rotate_deg` back to the local plots") and §6 (CVA level order);
   - the Task C plan's issues table marks the four as fixed.

### Modified files (uncommitted)

- **Package:** `R/biplot_grid.R`, `R/local_cf.R`, `R/model_fit.R`, `R/plot_biplot.R`,
  `R/predict_utils.R`, `R/projection.R`, `R/result.R`, `R/shapley.R`, `R/surrogate.R`, 10 files
  in `man/`, and the tests `test-biplot_grid.R`, `test-local_cf.R`, `test-predict_utils.R`,
  `test-projection.R`, `test-result.R`.
- **Documents:** `2 implementation_summary.txt`, `CLAUDE.md`, `NEWS.md`, `README.md`,
  `development/reference/review_section4_to_6.md`, `development/reference/review_section9_to_15.md`,
  `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`, `progress.md`.
- **New:** `.claude/reference/gate2-code-fixes-plan.md`.

### Key decisions (all the user's)

- All four issues are fixed now, not deferred.
- `rotate_deg` is removed from the local plots, with an error if passed, rather than made to
  rotate the sparse cross.
- The xgboost fix includes the wrap-time objective and level checks. These read xgboost object
  attributes that are not documented API.

### Finding (recorded in CLAUDE.md §6)

- The first version built the confusion factor with levels `TP, TN, FP, FN`. That changed every
  Phase 3 number in the harness:
  - target 13's constrained CF went from 0.547 to 0.509;
  - Step 15 went from 0.5554 / 0.542 to 0.6094 / 0.517.
- The level order feeds `biplotEZ::CVA()`, and the change in `V` moves `.bl_rotate()`'s
  under-determined plane.
- The factor now keeps the alphabetical order of `as.factor()`, and a test pins it.

### Test status

- 278 PASS / 0 FAIL / 4 WARN / 0 SKIP, up from 250.
- `R CMD check`: Status OK, 0 / 0 / 0.
- Harness: the full console log matches a run of the committed `7aa1344` code line for line, apart
  from the timestamps, so every walkthrough number stands.

### Open points / blockers

- **No blockers.**
- **Parked by the user (2026-10-29), not to be raised again unless the user brings them up.** They
  go into Task D's known-issues note:
  - The README's Pima constraint example finds no counterfactual for its first predicted-diabetic
    patient (`solution_found = FALSE`, no error).
  - The harness plot title still says "Loan default (prior defaulters)"; the data are applicants
    with no prior default.
- 32 pre-existing non-ASCII characters remain in `R/` comments (allowed by `R CMD check`, but
  against CLAUDE.md §4).
- The pkgdown site (`docs/`) has not been rebuilt for sessions 2-3's man-page and vignette changes.
  Rebuild before a release (CLAUDE.md §5).

### Exact next step

1. The user reviews the session 3 fixes and says whether to commit them (nothing is committed or
   pushed without an explicit instruction).
2. Then start **Task C**: the Pima example in `scripts/paper1-pima/` (gitignored, local only).
   - Read `.claude/plans/commit-and-push-the-mutable-hearth.md` (Task C section; the four issues
     are now marked fixed) and `.claude/plans/boundary-logic-fix-prompt.md`.
   - Plan in `.claude/plans/` and wait for approval.
   - Build the example, ending at Gate 2.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-10-29 (session 2)** — walkthrough review completed:
  - Parameters tables in the three walkthroughs, and a variable-pruning section;
  - `"fixed"` restore after the search (`.restore_fixed()`, `fixed_check`);
  - one shared observed / counterfactual comparison table for the three prints;
  - Shapley "No effect" class and the greedy sparse CF (`sparse_path`).
  - Committed as `7aa1344` on 2026-10-30 (not pushed). 250 PASS / 0 FAIL / 4 WARN.

- **2026-10-29 (session 1)** — Paper 1 set-up:
  - created `paper_nr1_changes` (from `b7aee65`);
  - Task A (graphics-device restore) and Task B (`bl_select_target()` stores `class_obs`; Pima
    worked example in `documentation/select-target-class/`);
  - follow-ups P1–P3 (axis offset/tick pairing) and S1–S6 (stale-documentation sweep, about 150
    passages, walkthroughs renumbered to the script's steps);
  - Gate 1 reached. Committed as `f2a076b` (not pushed). 176 PASS / 0 FAIL / 4 WARN.

- **2026-09-10** — integrated the `documentation/biplot-rotation/` note-and-script bundle (copied
  from the `Explain XAI code` project), checked that the deleted rotation-plane note's content was
  carried over, and repointed every stale reference. Documentation only; committed 2026-10-29 as
  `b7aee65`.
- **2026-09-01** — wrote the biplot-rotation and pair-selection technical write-up from the PhD
  thesis excerpt (later superseded by the 2026-09-10 bundle). Dead end recorded in CLAUDE.md §5:
  a legacy binary `.doc` saved as `.docx` defeats pandoc, and Word COM automation is unreliable.
- **2026-08-03** — split `scripts/` into a gitignored personal playground and a new tracked
  `development/` folder (main dev/testing harness + its three code walkthroughs); repaired path
  references; committed and fast-forward-merged `method_developments` into `main` (Mac-tester gate
  consciously overridden), then pushed both branches. 136 PASS / 0 FAIL / 4 WARN unchanged.
- **2026-07-15** — script 17 (`scripts/17_loan-gam-pca-3d-4var-surface.R`) target-alignment
  rotation experiments surfaced the `.bl_rotate()` plane under-determination finding: bases
  differing by 4.6e-15 give rotation planes 65.9 deg apart, and only 46% of the plane lies in the
  `proj_pair` eigenplane it is named after. Documented in CLAUDE.md §9; no package code changed;
  folded into commit `c6ca185`.
- **2026-06-17** — `355e3e3` round the Shapley plot counterfactual label to 2 d.p. (observed
  stays 3 d.p.).
- **2026-06-16** — `0f47602` fix a stale import and non-ASCII chars in `biplot_grid.R`; doc updates.
- **2026-06-16** — raw-unit scaling support (`00b43a5`): `bl_set_scaling()` +
  `.bl_rescale_biplot_axes()` (raw-unit biplot axes), `.scale_to_raw()` (raw-unit feature values in
  Shapley/sparse/target prints; `kind="level"` vs `"delta"` is the correctness rule), and
  `bl_local$bl_counterfactual` mirroring `bl_target`. All display-only. Script 03 fixed to
  standardise on the training split only (was leaking test data). 77 -> 136 PASS.
- **2026-06-16** — `new_title = NA` extended to every biplot `plot()` method (`dc0c3c4`);
  established as a CLAUDE.md §5 convention.
- **2026-06-08** — editor-diagnostics fixes in `shapley.R`/`local_cf.R` (non-ASCII chars,
  multi-line `@importFrom`, ggplot2 NSE globals).
- **2026-06-08** — script cleanup + pkgdown site rebuild. Dead-end: pkgdown 2.2.0
  `build_home()` publishes root `*.md` (incl. `CLAUDE.md`/`progress.md`) as public pages —
  fixed via `llm-docs: false` + post-build `unlink()` (gotcha now in CLAUDE.md S5). Also
  accidentally deleted 4 untracked user files — the lesson is now CLAUDE.md S5 bulk-`rm` rule.
- **2026-06-08** — loan-default vignette with SHAP comparison (replaced Pima vignette;
  `fastshap`/`shapviz` added to Suggests).
- **2026-06-08** — biplotEZ visual-features note added to `review_section4_to_6.md`.
- **2026-06-08** — per-variable axis tick-mark counts (`ticks_var`/`ticks_n`) +
  `.make_ticks_vec()` helper.
- **2026-06-08** — cleanup + unified plot interface comments + reference-doc sync.
- **2026-06-05** — unified biplot `plot()` interface; helpers `.make_label_line_vec()`,
  `.apply_biplot_rotation()`; `no_points` -> `plot_points`; `label_dir` default `"Paral"`.
- **2026-06-05** — Phase 2 step merge (7+8 -> 7); Phase 3 renumber 10-17 -> 9-16; loan
  end-to-end smoke test passed.
- **2026-05-29** — implemented the three improvement plans: accuracy-correctness (5 fixes),
  usability-bug (6 fixes), documentation-gaps. All now done, not outstanding.
- **2026-05-25** — committed `bl_filter_outliers` merge (`d509201`); xgboost API fix
  (`xgboost()` -> `xgb.train()` with `params`/`learning_rate`).
- **2026-05-23** — synced `reference_data_prep_functions.md`; established memory<->reference
  sync rule (CLAUDE.md S2).
- **2026-05-23** — `bl_filter_outliers()` merged into `bl_prepare_data(hull_fraction=)`;
  all scripts/vignettes/tests updated (72 -> 77 PASS).
- **2026-05-23** — `.claude/reference/` made git-tracked; `bl_fit_model()` slimmed to
  GLM/SVM/NNET/RForest; commit-gate rule established (CLAUDE.md S5).
- **2026-05-22** — `rounding` -> `b_margin` migration; filter-order swap in
  `bl_find_local_cf()` (set_filters before train_ranges); usage maps added to review docs.
- **2026-05-20** — Mahalanobis distance migration (Plans A+B): metric `W`/`metric_inv` fields,
  `chol2inv` inversion, `distance=` param; `.bl_rotate()` bug fix; technical note added.
- **2026-05-19** — boundary-arrow migration to `bl_pick_point()`; XQuartz graphics-device
  guard; R CMD CHECK to 0E/0W/2N; new loan scripts; committed `bcfc3fb`.
- **2026-05-13** — created 3 code-walkthrough reference docs + reference infrastructure;
  committed `bf025e7`.
- **2026-05-08** — fixed 4 test failures (72 PASS); full code review; wrote 3 improvement
  plans (all later implemented on 2026-05-29).
