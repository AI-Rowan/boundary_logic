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
- **`scripts/` is gitignored** — the user's personal scratch playground, never published. Task C
  will write `scripts/paper1-pima/` there (local only, by design).
- **`.bl_rotate()` plane under-determination** — documented, deliberately not actioned.
  See `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1 and §6
  (a proposal, not implemented), and CLAUDE.md §9. Do not "fix" without an explicit instruction.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 176 PASS / 0 FAIL / 4 WARN** (the 4 WARN are pre-existing biplotEZ
  CVA 2-class notices). **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-10-29).
- **Verification:** run R as `.R` files via the PowerShell tool,
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <driver>.R`, with the driver calling
  `devtools::document(<pkg root>)` then `devtools::test(<pkg root>)`. The dev harness runs end to
  end through a scratchpad driver (CLAUDE.md §5).
- **`.claude/plans/` is gitignored** — plan files are local only. `.claude/reference/` *is*
  tracked (generic docs + archived plans).

---

## Session summary (2026-10-29, latest) — Paper 1: branch set-up, Tasks A and B, two follow-ups, documentation sweep

Goal: start the Paper 1 (JCGS) package work from `.claude/plans/boundary-logic-fix-prompt.md`
on a new branch, through Gate 1, then clear every pre-existing and stale issue found on the way
before Task C.

### What was done

1. **Commit, push, branch.** The 2026-09-10 documentation work was committed as `b7aee65` on
   `method_developments` and pushed. `paper_nr1_changes` was created there and pushed. The fix
   prompt moved from the repo root into `.claude/plans/` (it must not be published).
2. **Task A — graphics device.** `.build_hull_polygon()` (`R/hull_utils.R`) and
   `.bl_get_plot_bounds()` (`R/biplot_grid.R`) now close only their own temporary `png()` and
   re-activate the caller's device. New `tests/testthat/test-graphics_device.R`.
3. **Task B — target class.** `bl_select_target()` stores `class_obs`, read from the same row as
   the features; `plot.bl_local_result()` colours the target from it, not from
   `test_data[row_id]`. New `tests/testthat/test-local_cf.R`. Worked example bundle in
   `documentation/select-target-class/`: on Pima, 36 of 79 compared rows had the wrong class
   and 234 lost their confusion colour; the bug was confirmed on the unfixed code first.
4. **Follow-up P1–P3.** `.make_label_line_vec()` / `.make_ticks_vec()` (`R/plot_biplot.R`) now
   drop an unknown variable's value together with its name (values used to shift one axis
   along) and warn on a length mismatch; new `tests/testthat/test-plot_helpers.R`. The harness's
   reduced-model plot vectors were fixed (`person_age` was pruned). Stale line references in the
   walkthroughs were replaced with code anchors.
5. **Follow-up S1–S6, the stale-documentation sweep** (three audits found about 150 passages):
   - `DESCRIPTION` imports `grDevices`, `graphics`, `stats`, `utils`.
   - Roxygen corrected across `R/` (e.g. `bl_build_result(bl_model = NULL)` returns a
     model-free `bl_result`; `dist_z` / `dist_mahalanobis` are the winning pair's; sparse cross is
     green/yellow; XGB needs `xgb.train()`); `print`/`summary` of `bl_result` label the grid hull;
     `plot.bl_boundary()` "25%%" printing bug fixed.
   - `2 implementation_summary.txt`: all 31 findings (signatures, file tree, status, dependencies).
   - Walkthroughs renumbered to the script's steps (`review_section4_to_6.md` 6.1–6.6;
     `bl_select_target()` now in Step 10 of `review_section9_to_15.md`), every number and console
     sample taken from a real harness run, written in present tense.
   - README examples fixed (they errored: `rounding =`, `bl_filter_outliers()` after
     `bl_prepare_data()`, a Pima example on iris objects); vignette prose, NEWS and CLAUDE.md §9
     corrected; harness comments fixed.
6. **CLAUDE.md** gained the branch (§10) and five working notes (walkthrough style, device
   restore, running R and the harness, `bl_pick_point` data rows).

### Key decisions

- Work on `paper_nr1_changes`, not `method_developments` (user's instruction; overrides the prompt).
- `class_obs` prefers a `class` column, otherwise the first non-feature column (a refinement on
  the prompt's rule, to avoid picking up an ID column).
- Walkthroughs describe current code and methodology only; NEWS is the change log.
- Code-behaviour bugs found during the sweep are **not** fixed yet; they go to Task C's Gate 2
  issues table: `used_in_sparse` stays `TRUE` for a reverted "fixed" variable; the sparse-CF cross
  ignores `rotate_deg`; `cutoff` is not passed to `bl_build_grid()`; `xgboost::xgboost()` models
  fail in `.pred_function()`; no guard for an empty outcome group.

### Test status

176 PASS / 0 FAIL / 4 WARN; `R CMD check` 0 / 0 / 0; every README chunk runs; the harness runs
clean in about 1 minute. All of this is committed on `paper_nr1_changes` (not pushed).

### Open points / blockers

- **The user is reviewing the three `development/reference/review_section*.md` walkthroughs.**
  Task C does not start until they approve.
- The README's Pima constraint example finds no counterfactual for its first predicted-diabetic
  patient (`solution_found = FALSE`, no error). Ask whether to change the target or constraints.
- The harness plot title still says "Loan default (prior defaulters)"; the data are applicants
  with no prior default. Left as is (a code string); ask.
- 32 pre-existing non-ASCII characters remain in `R/` comments (allowed by `R CMD check`, but
  against CLAUDE.md §4).

### Exact next step

Wait for the user's review of the walkthroughs. Then, on the user's "approved", start **Task C**
(the Pima example in `scripts/paper1-pima/`) from the plan file, ending at Gate 2.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

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
