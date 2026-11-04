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
- **`paper_nr1_changes` is committed and pushed** (2026-11-04): `f2a076b` (session 1), `7aa1344`
  (session 2), `f4f3647` and `f0e3dd4` (session 3), `a89ff1f` (session 4), `49c9067` (session 5),
  `67a1424` (session 6: the Gate 2 fixes and the documentation audit), and the session-7 commit
  (Task D, CLAUDE.md, this file). Commit and push only when the user says so.
- **Tasks A–D are done; Gate 2 is closed.** Task E and Gate 3 remain.
- **`scripts/` is gitignored** — the user's personal scratch playground, never published. The
  Paper 1 Pima example (Task C) lives in `scripts/paper1-pima/` (local only, by design): its plan
  is `.claude/plans/paper1-pima-example.md`, its review page `scripts/paper1-pima/review.html`,
  the issues `issues.csv` (the page shows the open ones). How to run it: CLAUDE.md §5.
  - **Phases 1–2 use the 7-feature GAM. Phase 3 uses the GAM refitted on Pregnancies, Glucose,
    Insulin and Age** (`02b_prune-refit.R`), after the global search ranked DiabetesPedigreeFunction,
    SkinThickness and BMI lowest.
- **Known issues** (Task D): `documentation/known-issues.md` (9 open, 8 resolved). The user does
  not want any of them implemented for Paper 1.
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
- **Deferred / future work** — see `CLAUDE.md` Section 9 and `documentation/known-issues.md` (not
  duplicated here).
- **Baseline test status: 434 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are pre-existing biplotEZ CVA
  2-class notices). **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-11-04).
- **Verification:** run R as `.R` files via the PowerShell tool,
  `& "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" <driver>.R`, with the driver calling
  `devtools::document(<pkg root>)` then `devtools::test(<pkg root>)`. The dev harness runs end to
  end through a scratchpad driver (CLAUDE.md §5).
  - **The installed `boundarylogic` is stale** (installed 2026-09-10, before the tie-break fix).
    Drivers must `devtools::load_all()` the repo; `library(boundarylogic)` alone runs old code.
- **`.claude/plans/` is gitignored** — plan files are local only. `.claude/reference/` *is*
  tracked (generic docs + archived plans).

---

## Session summary (2026-11-04, session 7, latest) — the Pima review, the pruning and Task D

**Plan:** `.claude/plans/i-have-reviewed-the-robust-rossum.md`, archived as
`.claude/reference/known-issues-taskD-plan.md` with an implementation log. The user's decisions are
in `scripts/paper1-pima/pruned-model-issues-review.md`.

### What was done

1. **The user's review of `review.html`** (scripts only, `scripts/paper1-pima/`):
   - a paragraph before Figure 3: the biplot points are the 313 training rows; the targets come
     from the 79 test rows;
   - **pruning:** `02b_prune-refit.R` checks that DPF, SkinThickness and BMI are the three weakest
     in the global search (48.68, 44.11, 40.05), re-splits the same rows (same split, checked) and
     refits the GAM on four features. Phase 3 runs on it. Run 6 is Glucose and Insulin
     `"decrease"`;
   - **figures:** 3 = the 7-feature biplot; 4(a) = the importance boxplot ("global search"
     subtitle); 4(b) = Figure 3 on four features. The old 4(a) and every mention of the arrow issue
     are gone;
   - every biplot axis shows at least 2 ticks (`ticks_v = 4L`, more per axis where needed);
   - the Shapley subtitles give the target's and the counterfactual's predictions;
   - the review page follows the workflow (Phase 3 after the pruning) and shows only open issues.
2. **The user's decisions on the pruned model's issues:**
   - 4(a): Figure 7 is TP run 4, TN run 1 and TN run 6;
   - 5(b): run 2 is out of Table 1 (it equals run 1 for both targets);
   - 9(b): the review page shows the Mahalanobis check for both models;
   - 11, 12, 13: paper text; their package options go to Task D.
3. **Task D:** `documentation/known-issues.md` — 9 open issues (the approximation gap, only the
   winning pair's plane stored, `plot.bl_projection()`'s colours, the hard-coded `"fixed"`
   tolerance, the follow-up searches for `"fixed"` and `"integer"`, the README constraint example,
   the harness title) and 8 resolved ones with their commits. The arrow layer is left out (the
   user's choice). CLAUDE.md §9 points to the note.
4. **CLAUDE.md:** the note in §2; how to run the Pima example and how to get enough biplot ticks
   (§5); `bl_prepare_data()`'s split depends only on rows and seed (§6).
   `fixed-constraint-restore.md` §4 says its Pima example is the 7-feature run.

**Files changed (package repo):** `documentation/known-issues.md` (new),
`.claude/reference/known-issues-taskD-plan.md` (new), `CLAUDE.md`,
`documentation/fixed-constraint-restore.md` (in `67a1424`), `progress.md`.

**Key results (4-feature model):**
- **Targets:** TP = test row 23 (prediction 0.963), TN = test row 25 (0.374).
- **Runs without a counterfactual:** TP run 3 (no boundary point), TN runs 3 and 4 (the strict
  restore fails: 0.475, 0.487).
- **The TP counterfactual ages the patient** (Age 51 -> 59.58). The GAM's risk drops after 55, on
  19 of 313 training patients over 50.
- **Integer rounding:** for the TN target, Pregnancies 10.21 -> 11 triples √d²_M (0.14 -> 0.42).
- **The Mahalanobis check:** -35.2% on four features, against -27.4% on seven.
- **Figure 7(a) shows the same point as Figure 6(b)** (TP runs 4 and 5 give the same counterfactual
  on pair (1,2)).

### Test status

- **Tests:** 434 PASS / 0 FAIL / 4 WARN / 0 SKIP (re-run before the commit; no `R/` changes this
  session).
- **`R CMD check`:** 0 / 0 / 0 (2026-11-04).
- **Pima `tests.R`:** 13 of 13 PASS (the new check 1b covers the pruning; 8a–8c are identical).

### Open points / blockers

- **No blockers.**
- **For the user:** Figure 6(b) and Figure 7(a) show the same counterfactual — drop one, or point
  to it in the caption.
- **For the Task E plan:** on the 4-feature model only 6 pairs exist, so `max_pairs = 21` against
  `10` measures nothing there; the max_pairs effect must be measured on the 7-feature model.
- **Paper text for Task E:** open issues 5, 9, 11, 12, 13 (`issues.csv`), and Gate 2 issues 3
  and 6. Gate 2 issue 1 is removed and issue 10 is closed.
- **The pkgdown site** (`docs/`) has not been rebuilt for sessions 2–7. Rebuild before a release
  (CLAUDE.md §5).
- **The Claude Doc copy of the Task C plan** (claude.ai artifact) is out of date; left alone.
- 32 pre-existing non-ASCII characters remain in `R/` comments.

### Next session (exact next step)

1. **Start planning Task E** (in plan mode): the methodology-changes note for the supervisor,
   `documentation/methodology-changes/methodology-changes.qmd` rendered to `.docx`. The full
   specification is in `.claude/plans/boundary-logic-fix-prompt.md` (Task E, Parts 1 and 2), with
   the working plan's Task E section. The plan must settle:
   - which model the measurements run on (see the max_pairs point above);
   - the measurement script in `scripts/paper1-pima/` and its `method-change-effects.csv`;
   - the paper-text items above;
   - the thesis's random tie-break (inherited from the PhD code) and the CVA distance identity.
2. **Then Gate 3**, and after it the merge, the release tag, and the Pima scripts' move to the
   paper repo (reinstall the package first).

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-11-02 to 2026-11-04 (session 6)** — Task C (the Paper 1 Pima example, Gate 2: 10 issues),
  then the Gate 2 fixes: `dist_mahalanobis` in `bl_find_boundary()`, strict `"fixed"`, the
  `"integer"` constraint (`round_to` deprecated), `"observed -> counterfactual"` Shapley labels;
  then a documentation audit. Committed as `67a1424` on 2026-11-04. 434 PASS.

- **2026-10-31 (session 5)** — `bl_apply_scaling()` (Step 16 raw units), `max_pairs` default 21,
  the Phase 2 Mahalanobis check describes the plot, fixed 2-dp print tables; plan
  `step16-pairs-check-decimals-plan.md`. Committed as `f0e3dd4`, `a89ff1f`, `49c9067`. 348 PASS.

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
