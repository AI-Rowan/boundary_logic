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
- **`paper_nr1_changes` state:**
  - **Committed, not pushed:** `f2a076b` (session 1), `7aa1344` (session 2) and `f4f3647`
    (session 3, Gate 2 fixes).
  - **Uncommitted** in the working tree:
    - session 3's Part B (nearest-vertex fix) and Part C (rotation folder);
    - session 4's merge of the two notes;
    - session 5's package changes (`max_pairs = 21`, the Mahalanobis check, the 2-dp print
      format, `bl_apply_scaling()` and the Step 16 fix).
  - **Staged:** the rename `documentation/mahalanobis_technical_note.md` →
    `documentation/boundary-distance-and-robustness.md`, from `git mv`.
  - Commit only when the user says so.
- **`scripts/` is gitignored** — the user's personal scratch playground, never published. Task C
  will write `scripts/paper1-pima/` there (local only, by design).
- **The boundary searches are deterministic** since 2026-10-30 (exact nearest vertex;
  `documentation/nearest-vertex-tie-break.md`).
- **The rotation and pair selection are documented in one note:**
  `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`.
  - **Its parts:**
    - §2: the original PhD method;
    - §3: the implementation;
    - §4: the problems;
    - §5: the proposal (not implemented);
    - §6: future developments;
    - Appendix A: the measurements.
  - **Phase 2 distances** are in `documentation/boundary-distance-and-robustness.md`.
  - **The plane under-determination** (§4.1) stays documented and deliberately not actioned. Do
    not "fix" it without an explicit instruction (CLAUDE.md §9).
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 348 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are pre-existing biplotEZ CVA
  2-class notices). **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-10-31, session 5).
- **Verification:** run R as `.R` files via the PowerShell tool,
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <driver>.R`, with the driver calling
  `devtools::document(<pkg root>)` then `devtools::test(<pkg root>)`. The dev harness runs end to
  end through a scratchpad driver (CLAUDE.md §5).
  - **The installed `boundarylogic` is stale** (installed 2026-09-10, before the tie-break fix).
    Drivers must `devtools::load_all()` the repo; `library(boundarylogic)` alone runs old code.
- **`.claude/plans/` is gitignored** — plan files are local only. `.claude/reference/` *is*
  tracked (generic docs + archived plans).

---

## Session summary (2026-10-31, session 5, latest) — Step 16 raw units, `max_pairs = 21`, the Phase 2 check, 2-dp tables

**Plan:** `.claude/plans/step16-pairs-check-decimals.md` (approved, implemented). Archive it with
session 4's plan (`.claude/plans/can-you-compare-the-nifty-flask.md`) in `.claude/reference/` when
the user commits.

### What was done (the user's five requests)

1. **PCA selector default: deferred** (documentation only).
   - **The decision:** Paper 1 compares by Mahalanobis distance only, under CVA, where both
     selectors are identical. The PCA default (`W_binary`) is revisited only once the rotation
     problem is fully addressed.
   - **Recorded in:** rotation note §4.4 and §6.2, and CLAUDE.md §6 and §9.
2. **`max_pairs` default 10 → 21** (`bl_find_local_cf()`).
   - **The user's choice:** 21, after being told that 15 never reaches component 7.
   - **What it covers:** every pair for p ≤ 7, so Paper 1's Pima search (p = 7) tries all 21
     pairs, including the thesis's (1,7). From p = 8 the cap still leaves pairs out.
   - **Effect:** the loan harness (p = 4) is unchanged. The README Pima example (p = 6) now tries
     all 15 pairs, and its result is unchanged.
   - **Also updated:** the local plan files `boundary-logic-fix-prompt.md` (Task C set-up,
     figure 6, a Task E note) and `commit-and-push-the-mutable-hearth.md`.
3. **The Phase 2 Mahalanobis check now describes the plot** (`plot.bl_boundary()`).
   - **The fix:** its "diagonal" is `sum(vec_to_boundary_sd^2)` (1/W_jj, as plotted), instead of
     `diag(W^-1)`.
   - **New output:** a renamed print, and a new `mahalanobis_check` return field.
   - **Measured on the harness:**
     - Step 7: +17.7% (the old diagnostic printed −42.6%, mostly the mismatch);
     - Step 7b: +15.9%;
     - neither reaches the 25% note.
   - **Tests:** new `tests/testthat/test-boundary_plot.R`.
4. **Fixed 2-dp print tables.**
   - **Covered:** the comparison tables of `print(bl_local)`, `print(bl_shapley)` and
     `print(bl_sparse)`, the Target/Counterfactual blocks, the Shapley column and the sparse path.
   - **The format:**
     - fixed notation with 2 decimals;
     - widened per feature row (up to 4) where a change would print as 0.00;
     - the sparse path's predictions keep 4 decimals within 0.005 of the cutoff.
   - **Helpers:** `.fmt_digits()`, `.fmt_fixed()`, `.format_cf_table()`, `.print_values()`.
   - `print.bl_target()` and `print.bl_counterfactual()` gain `digits = NULL`.
   - **The cause of the scientific notation:** the old harness log had it only in Step 16, from
     the raw-unit bug (a loan of 72 million).
5. **Step 16 raw units fixed.**
   - **New export:** `bl_apply_scaling(data, bl_result)` converts raw-unit data to the model's
     units.
   - **New warning:** `bl_select_target()` warns when a data-frame target lies outside
     `train_ranges`, with a `bl_apply_scaling()` hint when the target looks like raw units.
   - **The harness:** Step 16 converts the applicant.
   - **The new result:**
     - Z 0.0135 (was 115.76);
     - the applicant is predicted 0.296, class 0 (was 0.787, class 1);
     - pair (1,2) wins, at √d²_M 1.7026;
     - the counterfactual's credit score is 489.44, at the training minimum (489).
   - **The applicant still warns:** credit score 420 is below the training range (489–754).

**Documentation:**
- NEWS (one Breaking-change bullet, one New-feature bullet, two Bug-fix bullets);
- `2 implementation_summary.txt`;
- `1 Foundation intro documents.txt` §10.7;
- CLAUDE.md (§2, §3, §5, §6, §9);
- both notes;
- walkthroughs `review_section4_to_6.md`, `review_section7_to_8.md` and `review_section9_to_15.md`
  (every quoted block refreshed from the harness).

`review_section9_to_15.md` also had a stale sentence saying repeated calls can differ; it was
replaced.

### Test status

- **Tests:** 348 PASS / 0 FAIL / 4 WARN / 0 SKIP.
- **`17_3`:** every assertion passes; reference values unchanged (1.147879 / 1.312277).
- **Harness:** runs end to end; no scientific notation in the log.
- **`R CMD check`:** 0 errors / 0 warnings / 0 notes.

### Open points / blockers

- **No blockers.**
- **Parked by the user (2026-10-29), not to be raised again unless the user brings them up.** They
  go into Task D's known-issues note:
  - the README's Pima constraint example finds no counterfactual for its first predicted-diabetic
    patient;
  - the harness plot title says "Loan default (prior defaulters)".
- **The vignette's external applicant** (raw units, no scaling, so correct) now gets the
  target-range warning too: its credit score of 420 is below the training range.
- 32 pre-existing non-ASCII characters remain in `R/` comments.
- **The pkgdown site** (`docs/`) has not been rebuilt for sessions 2–5. Rebuild before a release
  (CLAUDE.md §5).

### Next session (exact next steps)

**0. Commit decision.** Sessions 3 (Parts B and C), 4 and 5 are uncommitted; the session-4 rename
is staged. Ask the user. A natural split is five commits:

- the nearest-vertex fix (code, test, note);
- the rotation folder brought to the current state;
- the note merge and Phase 2 note;
- `max_pairs = 21`, together with the PCA deferral;
- the Mahalanobis check, the 2-dp format, and `bl_apply_scaling()` with the Step 16 fix.

Some files carry changes from several of these (NEWS, CLAUDE.md, the implementation summary, the
walkthroughs), so splitting them needs `git add -p`, or fewer, larger commits.

**1. The remaining rotation concerns.** The ranked list is in the rotation note §6.2. The first is
cross-platform reproducibility (§4.1), which matters for Paper 1. The PCA selector waits for the
rotation work.

**2. Then Task C** (the Pima example in `scripts/paper1-pima/`, gitignored).

- Read `.claude/plans/commit-and-push-the-mutable-hearth.md` and
  `.claude/plans/boundary-logic-fix-prompt.md`; the Task C set-up now uses `max_pairs = 21L`.
- Plan, wait for approval, and build to Gate 2.
- The Task E pair-selector measurement: under CVA the two selectors are identical, so report the
  identity. Paper 1 uses Mahalanobis comparisons only.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-10-31 (session 4)** — the rotation and Mahalanobis notes were reconciled and merged
  (`biplot_rotation_and_pair_selection.md`, rewritten; plan `can-you-compare-the-nifty-flask.md`).
  - The Phase 2 note was renamed to `boundary-distance-and-robustness.md` (staged `git mv`).
  - Under CVA the two pair selectors are identical: d²_M = (n−g)·dist_z²; `17_3` §11 asserts it.
  - Documentation only; 314 PASS. Uncommitted.

- **2026-10-30 (session 3)** — three pieces of work:
  - Gate 2 code fixes, committed as `f4f3647`: `rotate_deg` removed from the local plots; the
    grid cutoff taken from the model; `xgboost()` models; empty CVA groups.
  - The deterministic nearest boundary vertex (`ties.method = "first"`;
    `documentation/nearest-vertex-tie-break.md`, new `test-boundary.R`). Uncommitted.
  - The rotation folder brought to the current state (re-measured; `17_3` checks asserted).
    Uncommitted.
  - 314 PASS / 0 FAIL / 4 WARN; `R CMD check` 0 / 0 / 0.

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
