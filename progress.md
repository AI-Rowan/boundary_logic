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
- **`paper_nr1_changes` state:** `f2a076b` (session 1) is committed but **not pushed**. All of
  session 2's work (below) is **uncommitted** in the working tree. Commit only when the user
  says so.
- **`scripts/` is gitignored** — the user's personal scratch playground, never published. Task C
  will write `scripts/paper1-pima/` there (local only, by design).
- **`.bl_rotate()` plane under-determination** — documented, deliberately not actioned.
  See `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1 and §6
  (a proposal, not implemented), and CLAUDE.md §9. Do not "fix" without an explicit instruction.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 250 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are pre-existing
  biplotEZ CVA 2-class notices). **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-10-29,
  session 2).
- **Verification:** run R as `.R` files via the PowerShell tool,
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <driver>.R`, with the driver calling
  `devtools::document(<pkg root>)` then `devtools::test(<pkg root>)`. The dev harness runs end to
  end through a scratchpad driver (CLAUDE.md §5).
- **`.claude/plans/` is gitignored** — plan files are local only. `.claude/reference/` *is*
  tracked (generic docs + archived plans).

---

## Session summary (2026-10-29, session 2, latest) — Walkthrough review completed; fixed restore, comparison tables, greedy sparse CF

Goal: the user's review of the three `development/reference/` walkthroughs. Each review comment
was resolved in the harness, the walkthroughs and, where the method needed it, the package.
**The review of `development/` is complete (user, 2026-10-29).**

### What was done

1. **Parameters tables.** Every function the three walkthroughs introduce has a
   `Parameter | Default | This call | Effect` table after its call. The table lists only the
   arguments the script sets and the defaults that change the result. All 70 Default cells were
   checked against `formals()`. The convention is in CLAUDE.md §2.
2. **Harness, Phase 2.**
   - Steps 8 and 8b (surrogates) now sit together after Step 7b, with a comment that nothing
     downstream uses them.
   - The pruning block prints the overall robustness score, labelled "not used here".
   - A comment gives why `person_income` is dropped and `credit_score` kept.
   - `review_section7_to_8.md` gained a "Variable pruning" section: `bl_robustness()` →
     `feature_cols_v2` → refit → `bl_results_v2`.
3. **Harness, Phase 3 (user edits).** The Step 10 constraints are now `loan_amnt = "fixed"`,
   `loan_int_rate` / `credit_score` `"increase"`. Step 15 tick counts are `c(4, 10, 10, 10)`.
   Step 15 prints `bl_sparse_free`. `round_to = NULL` was removed from the
   `bl_find_sparse_cf()` calls.
4. **Package: `"fixed"` restore** (`R/local_cf.R`).
   - The search keeps the ±0.5 model-unit tolerance, now the constant `.fixed_tol`.
   - After the search, `.restore_fixed()` sets each moved fixed feature back to its observed
     value. It keeps that point when the point still crosses the boundary and is within
     `train_ranges`.
   - New field `bl_local$fixed_check`.
   - Technical note `documentation/fixed-constraint-restore.md`.
5. **Package: comparison tables.**
   - `print(bl_local)` ends with an observed / counterfactual / change / constraint table, plus a
     note explaining any moved fixed feature.
   - `print(bl_shapley)` and `print(bl_sparse)` use the same layout.
   - All three tables share one builder, `.cf_compare_table()`.
6. **Package: sparse CF** (`R/shapley.R`).
   - A Shapley value of 0 is now `"No effect"`. The PhD rule counted 0 as "Supports".
   - `bl_find_sparse_cf()` is greedy: it changes the "Supports" features largest |Shapley| first,
     re-scores the model after each, and stops at the first crossing. New field `sparse_path`,
     printed.
   - `used_in_sparse` is `FALSE` for fixed features (this closes that Gate 2 item).
   - `round_to` is unchanged; its issues are recorded in `documentation/sparse-counterfactual.md`
     §3.
7. **Documents:**
   - `2 implementation_summary.txt` §4.8–4.9; NEWS; `1 Foundation intro documents.txt`; README;
     the loan vignette's prose;
   - the three walkthroughs (every number from a harness run);
   - plan archives in `.claude/reference/` (`fixed-constraint-restore-plan.md`,
     `sparse-counterfactual-plan.md`);
   - CLAUDE.md §2, §5, §6 and §9.

### Modified files (all uncommitted)

- **Package:** `R/local_cf.R`, `R/shapley.R`, `tests/testthat/test-local_cf.R`, and
  `man/` (`bl_find_local_cf`, `bl_find_sparse_cf`, `bl_shapley`, `plot.bl_shapley`,
  `print.bl_local_result`, `print.bl_shapley`, `print.bl_sparse_result`, `set_filters`).
- **Harness and walkthroughs:** `development/03_loan_status_Boundary_Logic.R`, and
  `development/reference/review_section4_to_6.md`, `review_section7_to_8.md`,
  `review_section9_to_15.md`.
- **Documents:** `1 Foundation intro documents.txt`, `2 implementation_summary.txt`, `NEWS.md`,
  `README.md`, `vignettes/Boundary_Logic_loan_default_workflow.Rmd`, `CLAUDE.md`, `progress.md`.
- **New:** `documentation/fixed-constraint-restore.md`, `documentation/sparse-counterfactual.md`,
  `.claude/reference/fixed-constraint-restore-plan.md`,
  `.claude/reference/sparse-counterfactual-plan.md`.

### Key decisions (all the user's)

- **Parameters tables:** selective, in one four-column format.
- **Pruning:**
  - It reads the per-variable totals; the overall robustness score is printed but not used.
  - `person_income` is dropped as redundant: `loan_percent_income` = `loan_amnt / person_income`
    rounded to 2 d.p. in 44,992 of 45,000 rows.
  - `credit_score` is kept as actionable.
- **`"fixed"`:** the search tolerance stays at 0.5 model units, then the post-search restore.
  `dist_*` and the per-pair tables stay the search's (they describe the pair choice).
- **"No effect":** `abs(shapley) < 1e-8`, a departure from the PhD `>= 0` rule.
- **Greedy sparse CF:** it re-scores the model rather than summing Shapley values (a Shapley value
  is an average over coalitions). It is not guaranteed minimal.
- **`round_to`:** unchanged (model units, nearest). Investigation is deferred (CLAUDE.md §9).
- **Tables:** rows follow `var_names`. `print(bl_shapley)` wraps its `contribute` column at an
  80-column console; this is cosmetic and left as is.

### Test status

- 250 PASS / 0 FAIL / 4 WARN / 0 SKIP, up from 176 (48 + 26 new expectations in
  `test-local_cf.R`).
- `R CMD check`: 0 / 0 / 0.
- The harness runs clean (about 1 minute). Target 13:
  - the restored counterfactual keeps `loan_amnt` at 12,000 and predicts 0.547;
  - the sparse path is `loan_int_rate` (0.360), then `loan_percent_income` (0.547, crosses);
  - Step 15 (no constraints) gives 0.5554 / 0.542.

### Open points / blockers

- **No blockers.**
- **Gate 2 issues table (Task C) still holds four code issues:**
  - the sparse-CF cross ignores `rotate_deg`;
  - `cutoff` is not passed to `bl_build_grid()`;
  - `xgboost::xgboost()` models fail in `.pred_function()`;
  - no guard for an empty outcome group.
- **Parked by the user (2026-10-29), not to be raised again unless the user brings them up.** They
  go into Task D's known-issues note:
  - The README's Pima constraint example finds no counterfactual for its first predicted-diabetic
    patient (`solution_found = FALSE`, no error).
  - The harness plot title still says "Loan default (prior defaulters)"; the data are applicants
    with no prior default.
- 32 pre-existing non-ASCII characters remain in `R/` comments (allowed by `R CMD check`, but
  against CLAUDE.md §4).
- The pkgdown site (`docs/`) has not been rebuilt for this session's man-page and vignette
  changes. Rebuild before a release (CLAUDE.md §5).

### Exact next step

Start **Task C**: the Pima example in `scripts/paper1-pima/` (gitignored, local only).
1. Read `.claude/plans/commit-and-push-the-mutable-hearth.md` (Task C section) and
   `.claude/plans/boundary-logic-fix-prompt.md`.
2. Plan in `.claude/plans/` and wait for approval.
3. Build the example, ending at Gate 2 with the issues table above.

Before starting, ask the user whether to commit session 2's work first (nothing is committed or
pushed without an explicit instruction).

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

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
