# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Next session: start here

**Update the Paper 1 plan**, then write its algorithm boxes. Start on branch `paper_nr1_changes`
(CLAUDE.md §10). The plan to update is
`../Explain XAI code/reference-materials/Paper_nr1/Exploring_boundaries/.claude/plans/paper-1-plan.md`
(last edited 2026-10-29). Plan the update first; what it must take in is under "Outstanding items",
Paper 1, below. The algorithm boxes have their own plan, `.claude/plans/before-updating-the-paper-wiggly-wave.md`
(stage 2 of `documentation/algorithms.md`).

---

## Outstanding items (consolidated 2026-11-07)

Every plan in `.claude/plans/` was audited against the code and git on 2026-11-07; 33 were finished
and deleted (the ones without an archive copy were archived to `.claude/reference/` first). These are
the items still open. Nothing here is a blocker.

### Paper 1 (next session; mostly in the paper repo)

1. **Update `paper-1-plan.md`.** It must take in:
   - **parameters and wording:** `max_pairs` is 21, not 10 (paper plan lines 11, 38, 76, 242);
     remove "no hull after rotation" and "call-history dependence" (both resolved);
   - **the method:** the convex hull of the training data on each rotated plane, with the empty-region
     rule (B10; `documentation/local-convex-hull.md`); the switches in `set_filters()` (the
     unconstrained counterfactual is one call); the variable importance measures the in-plane gap,
     while the global distances are reported only (`algorithms.md` F3);
   - **figures and Table 1:** figure 6 is the dot chart only (no 6(b)); figure 7(a) is the TP target's
     run 4 searched without the hull; figures 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) changed; Table 1 has
     runs 1, 3, 4, 5, 6 and 7;
   - **the TP narrative:** plausibility against actionability (the issue 11 sentence chosen in
     `future-improvements-review.md`), replacing the "older patient" story;
   - **§5 (local interpretation):** the hull as part of the method, and the A2, B3 and B7 sentences;
   - **§8 (properties and limitations):** the A1 result (upper-bound sentence, `05_full-space-gap.R`)
     and the B2 platform sentence;
   - **§9 (further work):** the four themed paragraphs built from the agreed sentences (G1 in the form);
   - **the Pima paper sentences** for issues 5 (why run 2 is not in Table 1), 9 (the Mahalanobis
     check), 11, 12 (fixed runs often have no counterfactual) and 13 (the integer rounding cost):
     decided in `scripts/paper1-pima/issues.csv` and the form, G1;
   - **the final thesis's page numbers** (CLAUDE.md §5 anchors);
   - **tasks T3** (paper repo infrastructure, TinyTeX) **and T4** (port `scripts/paper1-pima/` to the
     paper repo, switch `00_setup.R`'s `devtools::load_all("../..")` to `library(boundarylogic)`, check
     `numbers.rds` matches); figures 1–2 (schematics, made in the paper repo).
   Source: the archived Pima and task plans (`.claude/reference/paper1-pima-example-plan.md`,
   `paper1-task-spec.md`), `future-improvements-review.md`.
2. **The paper's algorithm boxes** (stage 2 of `documentation/algorithms.md`): four LaTeX boxes for
   CVA in the paper's notation, with a traceability table. Plan:
   `.claude/plans/before-updating-the-paper-wiggly-wave.md`.
3. **Hand-over to the paper writer:** the test counts (561 PASS; `R CMD check` 0/0/0), the commit
   (`0689341` on `main`), `scripts/paper1-pima/review.html`,
   `documentation/methodology-changes/methodology-changes.docx`, the Pima issues and decisions
   (`issues.csv`), and the install command. Source: `.claude/reference/paper1-task-spec.md`, "Report back".
4. **Gate 3: the supervisor reviews the methodology note.** Afterwards decide whether the note lists the
   design-document-only items B5 (a test of the CVA identity) and B8 (a `pairs =` argument); B21–B23
   are in CLAUDE.md §9. Source: `.claude/reference/future-improvements-hull-plan.md` (Phase 4).

### Release and branches: done 2026-11-14 (session 13)

5. ~~Tag `v0.2.0`~~ **Done:** annotated tag `v0.2.0` on `0689341`, pushed to origin. The paper
   installs it with `devtools::install_github("AI-Rowan/boundary_logic@v0.2.0")`.
6. ~~`method_developments`~~ **Fast-forwarded** to `main` and pushed (kept as the development branch;
   CLAUDE.md §10).
7. ~~Check the public install~~ **Done:** `remotes::install_github("AI-Rowan/boundary_logic@v0.2.0",
   build_vignettes = TRUE)` without a GitHub token, into an empty library placed first on
   `.libPaths()` (dependencies from the user library), installs 0.2.0 (`RemoteSha` `0689341`), 23
   exports, and the loan-default vignette's HTML (1.1 MB) opens with `vignette()`.
8. ~~The old Pima review Claude Doc~~ **Deleted** ("Task C plan: the Paper 1 Pima example",
   `A3XvG3gydFu8DUm2umqaLc`); `scripts/paper1-pima/review.html` replaces it.

### Housekeeping in this repo

9. ~~Commit the rotation-script loading fix~~ **Committed and pushed** (2026-11-14) with the plan
   archives and this file; `main` fast-forwarded to it.
10. **Re-render the Pima review page** (`scripts/paper1-pima/review.qmd`): it still names `ef014ca`,
    but the outputs are from `17d3c18` (CLAUDE.md §5).
11. **Refresh the root documents' trees and counts** (they drifted after the 2026-11-04 audit):
    - `2 implementation_summary.txt`: "Last updated" date; the test tree (15 files: add
      `test-distance_by_method.R`, `test-hull_utils.R`); the documentation tree (add
      `local-convex-hull.md`, `algorithms.md`, `methodology-changes/`; drop the deleted presentation
      PDF); the expectation count (561 PASS, from a real `devtools::test()` run);
    - `1 Foundation intro documents.txt`: the "last revised" date; §8's file tree (it still lists the
      deleted `scripts/00_pima_*.R`, `docs/` files that do not exist, 13 test files; it omits
      `pick_point.R`, `datasets.R`, `loan_data.csv`, `development/`, `documentation/`,
      `vignettes/*.Rmd.orig` and `precompile.R`); §9's branch list (add `paper_nr1_changes`).
12. **A test gap:** `tests/testthat/test-set_scaling.R` has no raw-unit relabel round-trip case for a
    PCA with `standardise = FALSE`, and checks only that the label midpoint lies near the centre;
    assert exact labels (`std * scale + center`) at matching ticks (source: the 2026-06-16 raw-scale
    axes plan, `.claude/reference/scaling_axis_relabel_plan.md`).

---

## Standing items / carry-overs

- **Paper 1 work is on branch `paper_nr1_changes`** (branched from `method_developments` at
  `b7aee65` on 2026-10-29, pushed with upstream). Start Paper 1 sessions on that branch
  (CLAUDE.md §10).
- **`paper_nr1_changes`, `main` and `method_developments` (local and origin) are all at the
  2026-11-14 housekeeping commit** (rotation-script loading fix, plan archives), one documentation-only
  commit after `0689341`. **Tag `v0.2.0` = `0689341`** (the package code is the same at both). Commit
  and push only when the user says so.
- **The installed `boundarylogic` is 0.2.0** (reinstalled 2026-11-07 from the working tree, equal to
  `main`). Reinstall after each merge to `main` (CLAUDE.md §5).
- **`scripts/` is gitignored** (the user's playground). The Paper 1 Pima example lives in
  `scripts/paper1-pima/` (local only): review page `review.html`, issues `issues.csv`, the
  approximation gap `05_full-space-gap.R`, the findings effects `06_findings-effects.R`. How to run it:
  CLAUDE.md §5. Phases 1–2 use the 7-feature GAM; Phase 3 the GAM refitted on Pregnancies, Glucose,
  Insulin and Age. Table 1 has runs 1, 3, 4, 5, 6 and 7.
- **Every open item is decided** (2026-11-05): `documentation/methodology-changes/future-improvements-review.md`
  records each decision (fix now, future improvement in the paper or the design document only, no
  action) and the agreed paper sentences.
- **The methodology note** (`documentation/methodology-changes/methodology-changes.qmd` / `.docx`)
  names commit `17d3c18`. It is organised as "Changes made" and "Future work", with no approval
  boxes. The supervisor has not seen it yet (item 4).
- **The local search's plausibility rule** (`set_filters(hull = TRUE)`): the convex hull of the
  training data on each rotated plane (the target is not part of it since 2026-11-06), and closed
  contours only if they enclose training data, the same empty-region rule as `bl_find_boundary()`.
  The plane under-determination (rotation note §4.1) stays documented and deliberately not actioned
  (CLAUDE.md §9).
- **Test status: 561 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are biplotEZ 2-class CVA notices).
  **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-11-06, session 12).
- **Verification:** run R as `.R` files (`Rscript.exe <driver>.R`), with drivers that
  `devtools::load_all()` the repo. The rotation scripts load the working tree themselves.
- **`.claude/plans/` is gitignored** (local only); `.claude/reference/` is tracked. Since the
  2026-11-07 consolidation it holds only the algorithm-boxes plan.

---

## Session summary (2026-11-14, session 13, latest)

- Release and branch housekeeping (items 5–9): the session-12 changes committed and pushed; tag
  `v0.2.0` on `0689341` pushed; `main` and `method_developments` fast-forwarded and pushed; the
  public `install_github(...@v0.2.0, build_vignettes = TRUE)` checked; the old Pima Claude Doc
  deleted. No package code changed: 561 PASS stands.
- **Next:** the Paper 1 plan update (start here, above).

---

## Session summary (2026-11-06 to 2026-11-07, session 12)

- **The vignette** (archived plan `.claude/reference/vignette-loan-workflow-plan.md`): rewritten for
  package users, following the harness on the full data; precomputed
  (`vignettes/Boundary_Logic_loan_default_workflow.Rmd.orig`, `vignettes/precompile.R`, hidden
  `*-claims` chunks); SHAP plots in raw units; no `bl_robustness()` step.
- **`documentation/algorithms.md`**: every stepwise algorithm (A1–A4, B1–B3, C1–C8, supporting S1–S9),
  their order and interaction, and findings F1–F8 with the user's decisions (§7).
- **The findings fixes** (archived plan `.claude/reference/algorithms-findings-fixes-plan.md`):
  - the empty-region rule judged on the training data in `bl_find_boundary()` (before: the searched
    rows) and added to the local search with the hull on (`.encloses_data()`);
  - each rotated plane's hull from the training data alone;
  - the local grid clamped to [0, 1]; the local plot widened to show a target outside the training
    data; `grid_col` split at the cutoff; `inside_data` removed; F1 noted in the code;
  - F3 (B1's gap against B2's) measured: the off-plane part is 59–80% of B1's gap, so switching B2 to
    it was tried and reverted at the gate; B1's distances are reported only;
  - Paper 1 unchanged (Table 1, every run); on the loan harness test row 13's constrained
    counterfactual (an empty closed contour held its nearest point) and the new applicant's move.
  - New `scripts/paper1-pima/06_findings-effects.R` (`findings-effects.csv`, read by the note).
- **The W note:** `boundary-distance-and-robustness.md` §1.2 says why the package computes W itself.
- **Committed and pushed** (on the user's instruction): `17d3c18`, `0689341`; `main` fast-forwarded.
- **The rotation scripts** (archived plan `.claude/reference/rotation-scripts-package-loading-plan.md`):
  they failed because `library(boundarylogic)` loaded the stale installed 0.1.0 (no
  `set_filters(hull =)`). They now load the working tree inside the repo and check the version
  elsewhere; the package was reinstalled (0.2.0). Output identical to before. Uncommitted (item 9).
- **The plans consolidated** (`/summarise-plans`, 2026-11-07): 34 plans audited; 33 deleted. Five had
  no archive copy and were archived first: `paper1-task-spec.md` and `paper1-working-plan.md` (the
  Paper 1 task specification and working plan), `paper1-pima-example-plan.md`,
  `vignette-loan-workflow-plan.md`, `rotation-scripts-package-loading-plan.md`. Every reference to a
  deleted plan (CLAUDE.md §2 and §10, `future-improvements-review.md`, this file, the archive banners)
  now points to its archive. Merged into the items above: session 9's "Next session" list (item 1) and
  the stale standing items (main, the tag, the reinstall, the test count). The algorithm-boxes plan is
  kept, trimmed to stage 2. Unsure, left open: Phase 4's B5/B8 (item 4), the public install check
  (item 7), the old Claude Doc (item 8).
- **Status:** 561 PASS / 0 FAIL / 4 WARN; `R CMD check` 0 / 0 / 0.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-11-06 (session 11)** -- the plausibility (`hull`) and feasibility (`train_ranges`) switches in `set_filters()`; the local plot colours its grid only inside the plane's hull; Pima figure 7(a) from run 4 without the hull; the final thesis PDF in `documentation/` and every citation checked. Committed and pushed as `9bc8e1c`. 527 PASS.

- **2026-11-05 (session 10)** -- `.poly_clip()` defined once (`R/hull_utils.R`), new `test-hull_utils.R`. Committed with session 11. 473 PASS.

- **2026-11-05 (session 9)** -- every open item decided (`future-improvements-review.md`, A1-A9, B1-B24); the convex hull on each rotated plane (B10); the approximation gap (`05_full-space-gap.R`, median 1.19); run 7; version 0.2.0; `hull_fraction` default 1; the methodology note reorganised; `known-issues.md` removed. Committed as `c13fe44`, `ee64c3a`. 465 PASS.

- **2026-11-04/05 (session 8)** -- Task E: the methodology-changes note (`methodology-changes.qmd`, `04_method-change-effects.R`, 5 self-checks); the §3.3 reason for keeping global CFs outside the ranges; the thesis distance rule by projection (`distance = "auto"`: Euclidean under PCA, Mahalanobis under CVA). Committed and pushed as `12b38e4`. 465 PASS.

- **2026-11-04 (session 7)** -- the user's review of the Pima page; pruning to four features (`02b_prune-refit.R`, Phase 3 on the 4-feature GAM); Task D (`documentation/known-issues.md`, 9 open / 8 resolved); figure 6(b) removed and run 5 added to Table 1. Committed as `4ecdc8c` and `540f9ad`, pushed. 434 PASS.

- **2026-11-02 to 2026-11-04 (session 6)** — Task C (the Paper 1 Pima example, Gate 2: 10 issues),
  then the Gate 2 fixes: `dist_mahalanobis` in `bl_find_boundary()`, strict `"fixed"`, the
  `"integer"` constraint (`round_to` deprecated), `"observed -> counterfactual"` Shapley labels;
  then a documentation audit. Committed as `67a1424` on 2026-11-04. 434 PASS.

- **2026-10-31 (session 5)** — `bl_apply_scaling()` (Step 16 raw units), `max_pairs` default 21,
  the Phase 2 Mahalanobis check describes the plot, fixed 2-dp print tables; plan
  `step16-pairs-check-decimals-plan.md`. Committed as `f0e3dd4`, `a89ff1f`, `49c9067`. 348 PASS.

- **2026-10-31 (session 4)** — the rotation and Mahalanobis notes were reconciled and merged
  (`biplot_rotation_and_pair_selection.md`, rewritten; plan `notes-merge-rotation-phase2-plan.md`).
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
