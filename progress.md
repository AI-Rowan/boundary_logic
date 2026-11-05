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
  wait for the user's instruction (part of the hand-over, below).
- **`paper_nr1_changes` is pushed up to `540f9ad`** (2026-11-04): `f2a076b` (session 1), `7aa1344`
  (session 2), `f4f3647` and `f0e3dd4` (session 3), `a89ff1f` (session 4), `49c9067` (session 5),
  `67a1424` (session 6: the Gate 2 fixes and the documentation audit), `4ecdc8c` and `540f9ad`
  (session 7: Task D, run 5). Session 8 (Task E and the PCA distance rule) is committed and pushed
  in the commit after `540f9ad` (2026-11-05). Commit and push only when the user says so.
- **Tasks A–E are done (session 8), committed and pushed.** The next step is the hand-over back
  to the paper writer (below).
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

## Session summary (2026-11-04, session 8, latest): Task E, the methodology-changes note

**Plan:** `.claude/plans/please-start-planning-task-frolicking-rossum.md`, archived as
`.claude/reference/methodology-changes-taskE-plan.md` with implementation notes.

### What was done

1. **Planning.** Three explorations were run: the thesis text, the PhD code against the package, and
   the package history. They gave five findings:
   - **Three baselines.** The thesis text, its Appendix A (pp.158-218), and the branch
     `original_PhD_code`, committed **2026-03-09, after submission**. A and B differ: A uses `which.min`,
     B the random `max.col()`; A builds the local hull from target + training data, B from target +
     test data.
   - **Cite the final PDF** (printed page = PDF page - 16). The rotation note and CLAUDE.md cite the
     30 Sept 2025 draft, about 20 pages higher.
   - **The (n-g) identity is not in the thesis**, whose W is the SSP (factor 1).
   - **The thesis Pima example did use pair (1,7)** (p.118).
   - **The package departs from the thesis text** on the local hull (p.117), on the global
     training-range filter (p.93), and on the pair selector under PCA (p.86).
2. **The user's structure (at plan review):**
   - a history of changes;
   - then the **core methodology** (11 changes + 2 unstated properties, one decision box each);
   - then the **code improvements** summarised (one box);
   - then future changes and interpretation points.
3. **Decisions:**
   - `max_pairs` is measured on the 4-feature model only (10 = 21);
   - the pair selector is the CVA identity, confirmed on Pima;
   - Pima issues 9, 12 and 13 go into their core entries; issues 5 and 11 go into Section 6.
4. **The measurement script `scripts/paper1-pima/04_method-change-effects.R`** (local) writes
   `documentation/methodology-changes/method-change-effects.csv` (511 rows, commit `540f9ad`). It has
   five self-checks, all passing:
   - the baseline reproduces `numbers.rds`;
   - the derived CVA labels reproduce V;
   - `max_pairs` 10 = 21;
   - the ratio d2_M/dist_z^2 = n-g = 309;
   - no empty values.
5. **The note** is `documentation/methodology-changes/methodology-changes.qmd`, rendered to `.docx`
   (about 6,300 words). Every measured number is an inline `r` call on the CSV, and a `claims` chunk
   asserts the qualitative statements. The prose follows the academic-writing-style skill; the style
   check is within range except the short-sentence share, inflated by "p." citations, and "This"
   chains at about 8%.
6. **§3.3 rewritten (the user's question, 2026-11-05).** It now gives the recorded reason for keeping
   global CFs outside the training ranges (`ea6a12f`, prompted by PCA). The CF is a rank-2
   back-projection (`B_z %*% tVr` sets the non-plotted components to their mean), so its features can
   leave the ranges even inside the hull, and a pre-filter picked a farther point. The section also
   adds:
   - that the package first followed the text (`e42fa03`);
   - that the commit's "CVA was unaffected" does not hold on Pima;
   - the distances-against-feasibility choice for the supervisor;
   - whether the targets' global CFs are within range (new CSV rows).

   The same reason is now in the walkthrough (`review_section7_to_8.md` Stage 8) and the
   implementation summary §6.11. Section 5 adds a per-row feasibility flag.
7. **The thesis distance rule under PCA (2026-11-05, the user's request).** Plan:
   `.claude/reference/pca-euclidean-distance-plan.md`.
   - **The change:** `distance = "auto"` (the new default) applies the thesis's rule by projection, Euclidean in
     X-space (each feature divided by its total SD) under PCA and Mahalanobis under CVA. It applies in:
     - `bl_find_local_cf()` (the pair comparison and the `"integer"` rounding);
     - `plot.bl_boundary()` / `bl_robustness()` (the total SD under PCA).
   - **Helpers and fields:** new private `.resolve_distance()` (`R/utils.R`) and `.cf_sq_distance()`
     (`R/local_cf.R`); new fields `dist_x` and `all_distances_x`, which the print and plot show.
   - **An explicit `"euclidean"`** now means the X-space distance; it was the plane distance.
   - **Paper 1 (CVA) is unchanged:**
     - `numbers.rds` is identical;
     - `tests.R` 13/13;
     - the loan harness's numbers are identical to `HEAD`.
   - **The rotation scripts** `17_1`-`17_3` pin `distance = "mahalanobis"` and reproduce the note's numbers.
   - **All notes updated:**
     - CLAUDE.md;
     - the implementation summary;
     - the Foundation document;
     - NEWS;
     - the rotation, boundary-distance and fixed-constraint notes, and the rotation folder's README;
     - `known-issues.md`;
     - the methodology note (§3.1, §3.5, the summary, Tables 1, 2 and 5; new M5-PCA / M1-PCA measurements:
       57 of 79 Pima test rows choose a different pair under PCA);
     - both walkthroughs (samples refreshed from the harness log).
8. **Corrections:**
   - CLAUDE.md §6: the thesis pages and the (n-g) factor;
   - NEWS.md: the "No effect" entry (class-0 targets only);
   - an erratum line in `.claude/reference/mahalanobis_implementation_plan.md`;
   - the working plan's tie-break bullet;
   - CLAUDE.md §2: pointers to the note and the archived plan.

**Key results (Pima):**
- **M1:** the thesis's total-SD rule gives the same ranks and prunes the same three features.
- **M3:** on the 4-feature model, 10/313 training and 5/79 test global CFs lie outside the training
  ranges (kept).
- **M4:** 10 = 21; all 6 pairs against (1,2): search distance TP 1.92 -> 1.50, TN 0.17 -> 0.14.
- **M5:** the same pair and CF under both selectors.
- **M10:** greedy is sparser by one variable for both targets.
- **K1** (`b_margin` 0.01): small changes; the rule picks TP row 1.
- **K2** (`outlie` 0.9): run 1 unchanged; the rule picks TN row 69.
- **K3** (thesis group order): V changes by 5e-17, yet the TP pair moves (2,3) -> (1,3) and its
  distance 1.42 -> 1.70. This is evidence of the plane's under-determination (§3.12).
- **2.12(a):** on the 4-feature model p - K = 1, so component 4 is fixed up to sign. The TN
  counterfactual uses pair (2,4).

**Files changed:**
- package repo:
  - `documentation/methodology-changes/` (new: qmd, docx, csv);
  - `development/reference/review_section7_to_8.md` (Stage 8);
  - `2 implementation_summary.txt` (§6.11);
  - `CLAUDE.md`;
  - `NEWS.md`;
  - `.claude/reference/mahalanobis_implementation_plan.md`;
  - `.claude/reference/methodology-changes-taskE-plan.md` (new);
  - `progress.md`;
  - **the PCA distance rule (2026-11-05):**
    - `R/utils.R`, `R/local_cf.R`, `R/boundary_plot.R`, `R/projection.R`, `R/result.R` and `man/`;
    - `tests/testthat/test-distance_by_method.R` (new) and `test-local_cf.R`;
    - `1 Foundation intro documents.txt`;
    - `documentation/biplot-rotation/` (the note, its README, and `17_1`/`17_2`/`17_3` pinned);
    - `documentation/boundary-distance-and-robustness.md`, `fixed-constraint-restore.md` and
      `known-issues.md`;
    - `development/reference/review_section9_to_15.md`;
    - `.claude/reference/pca-euclidean-distance-plan.md` (new);
- local only:
  - `scripts/paper1-pima/04_method-change-effects.R`;
  - the working plan `.claude/plans/commit-and-push-the-mutable-hearth.md`.

### Test status

- **Tests:** 465 PASS / 0 FAIL / 4 WARN / 0 SKIP, after the PCA distance rule (was 434).
- **Pima:** `run_all.R` gives `numbers.rds` identical to before (meta excluded); `tests.R` 13/13;
  `04_method-change-effects.R` passes all 6 self-checks.
- **The loan harness and the rotation scripts:** numbers identical to the `HEAD` baseline.

### Open points / blockers

- **No blockers.**
- **Committed and pushed** (2026-11-05): everything listed above, in one commit on
  `paper_nr1_changes`.
- **Local only, by design:** `scripts/paper1-pima/` (including `04_method-change-effects.R`) and
  `.claude/plans/`.
- **Not done, out of scope:** the rotation note's thesis citations stay on the 30 Sept 2025 draft (it
  says so). The pkgdown site is still not rebuilt.

### Next session (exact next step)

1. **Prepare the hand-over back to the paper writer** (the Paper 1 project,
   `../Explain XAI code/reference-materials/Paper_nr1/Exploring_boundaries/`). Plan it first. It is
   the report the task prompt asks for (`.claude/plans/boundary-logic-fix-prompt.md`, "A
   recommendation, then your report"):
   - the `devtools::test()` counts: 465 PASS / 0 FAIL / 4 WARN;
   - the files changed;
   - the Task C issues with the user's decision on each (`scripts/paper1-pima/issues.csv`);
   - the paths of `scripts/paper1-pima/review.html` and
     `documentation/methodology-changes/methodology-changes.docx`;
   - the commit SHA;
   - the reinstall command.

   **Also for the paper writer:**
   - the paper-text items: the methodology note's proposed sentences (issues 5, 9, 11, 12 and 13);
   - the changes since the prompt that affect the paper: the `"integer"` constraint, strict
     `"fixed"`, the pruned 4-feature Phase 3 model, and the 2026-11-05 PCA distance rule (no effect
     on Paper 1, which uses CVA);
   - how the Pima scripts move to the paper repo: the one `library(boundarylogic)` line in
     `00_setup.R`; `04_method-change-effects.R` writes into this repo's `documentation/`, so it stays
     here or takes `METHOD_CHANGES_OUT`.
2. **On the user's instruction, as part of the hand-over:**
   - merge `paper_nr1_changes` into `method_developments` (and `main`);
   - tag the release the paper pins, `install_github("AI-Rowan/boundary_logic@<tag>")`;
   - reinstall the package locally (the installed copy is stale);
   - copy the Pima scripts to the paper repo and re-run them against the installed package.
3. **The supervisor's decisions** on the methodology note (Gate 3) may come back with edits: apply
   them to the `.qmd`, re-run `04_method-change-effects.R` if a measurement changes, and re-render.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-11-04 (session 7)** -- the user's review of the Pima page; pruning to four features (`02b_prune-refit.R`, Phase 3 on the 4-feature GAM); Task D (`documentation/known-issues.md`, 9 open / 8 resolved); figure 6(b) removed and run 5 added to Table 1. Committed as `4ecdc8c` and `540f9ad`, pushed. 434 PASS.

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
