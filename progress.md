# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Standing items / carry-overs

- **Paper 1 work is on branch `paper_nr1_changes`** (branched from `method_developments` at
  `b7aee65` on 2026-10-29, pushed with upstream). Start Paper 1 sessions on that branch
  (CLAUDE.md §10).
- **`paper_nr1_changes` is pushed up to the sessions 10–11 commit** (2026-11-06: `.poly_clip()` once,
  the `set_filters()` switches, the local plot's hull trim, the final thesis in `documentation/`),
  after `c13fe44` (2026-11-05: the hull on each rotated plane, version 0.2.0) and its clean-up commit
  `ee64c3a`. Commit and push only when the user says so.
- **`main` was last merged and pushed on 2026-08-03; `method_developments` is at `b7aee65`.** Merging
  `paper_nr1_changes` back, the `v0.2.0` tag the paper will pin, and reinstalling the package
  locally wait for the user's instruction.
- **`scripts/` is gitignored** (the user's playground). The Paper 1 Pima example lives in
  `scripts/paper1-pima/` (local only): review page `review.html`, issues `issues.csv`, the
  approximation gap `05_full-space-gap.R`. How to run it: CLAUDE.md §5. Phases 1–2 use the 7-feature
  GAM; Phase 3 the GAM refitted on Pregnancies, Glucose, Insulin and Age. Table 1 has runs 1, 3, 4, 5,
  6 and 7.
- **Every open item is decided** (2026-11-05): `documentation/methodology-changes/future-improvements-review.md`
  records each decision (fix now, future improvement in the paper or the design document only, no
  action) and the agreed paper sentences. `known-issues.md` is removed (in git at `c13fe44`).
- **The methodology note** (`documentation/methodology-changes/methodology-changes.qmd` / `.docx`) is
  organised as "Changes made" (history, core methodology §3.1–3.12, code improvements) and "Future
  work", with no approval boxes. The supervisor has not seen it yet.
- **The local search applies the convex hull on each rotated plane** (since 2026-11-05;
  `documentation/local-convex-hull.md`). Since 2026-11-06 its switch, and a new one for the training
  ranges, are in `set_filters()` (`hull = TRUE`, `train_ranges = TRUE`); `set_filters = NULL` keeps
  both on. The plane under-determination (rotation note §4.1) stays documented and deliberately not
  actioned (CLAUDE.md §9).
- **Test status: 527 PASS / 0 FAIL / 4 WARN / 0 SKIP** (the 4 WARN are biplotEZ 2-class CVA notices).
  **`R CMD check`: 0 errors / 0 warnings / 0 notes** (2026-11-06, session 11).
- **Verification:** run R as `.R` files (`Rscript.exe <driver>.R`), with drivers that
  `devtools::load_all()` the repo. **The installed `boundarylogic` is stale** (2026-09-10).
- **`.claude/plans/` is gitignored** (local only); `.claude/reference/` is tracked.

---

## Session summary (2026-11-06, session 11, latest): the switches in `set_filters()`, the hull on the local plot, the final thesis

**Plan:** `.claude/plans/poly-clip-is-defined-twice-dynamic-garden.md` (the plan-mode file name),
archived as `.claude/reference/filter-switches-in-set-filters-plan.md` (IMPLEMENTED banner).

- **Done:** `set_filters(bl_target, ..., hull = TRUE, train_ranges = TRUE)` holds all three of the
  thesis's rules: actionability, plausibility and feasibility (p.116). `bl_find_local_cf(hull = )`
  is removed. The training ranges gain a switch (before, the only way off was
  `bl_result$train_ranges <- NULL`). `train_ranges = FALSE` also skips the range checks of
  `.restore_fixed()` and `.round_integers()` (new `train_ranges` argument). Both off with no
  constraint gives the unconstrained counterfactual (p.115).
- **Key decisions (user):**
  - the switches are named `hull` and `train_ranges`, documented as plausibility and feasibility;
  - the old `hull` argument is removed (0.2.0 untagged);
  - `set_filters = NULL` keeps both rules on, so no existing result changes
    (`local-convex-hull.md` §6 gives the reasons).
- **New:** private `.filter_switches()`, `.check_switch()`; result field `train_ranges`; `Not applied`
  print line when a rule is off; `print(<bl_filters>)` shows both switches; reworded no-solution message
  and target-range warning.
- **Callers updated:** harness (Step 10 sets both switches; the Step 15 title fixed), vignette,
  README, rotation scripts `17_1`-`17_3`, Pima `04` and `05`.
- **Verified unchanged against `HEAD`:**
  - the loan harness;
  - the rotation scripts;
  - the Pima outputs (`table1.csv`, `numbers.rds` apart from the message wording, `full_space_gap`,
    `method-change-effects.csv` apart from the meta rows).
- **Verification runs:** Pima `tests.R` 13/13, `04` 6/6, `05` both checks. `review.html` and the
  methodology note `.docx` were re-rendered.
- **New behaviour (harness):** for the external applicant (credit score below the training range),
  √d²_M is 1.7026 by default, 0.8521 with `train_ranges = FALSE` and 0.4750 with both off.
- **Docs updated:** NEWS, the implementation summary, the Foundation document, CLAUDE.md (§2, §4 new
  rule, §6), `local-convex-hull.md` (§6 new), `fixed-constraint-restore.md`, the rotation note, the
  methodology note, the walkthrough `review_section9_to_15.md`, `review.qmd`.
- **Then, the local plot shows the hull** (plan archived as
  `.claude/reference/local-plot-hull-trim-plan.md`):
  - `plot.bl_local_result(trim_to_hull = NULL)` colours the grid only inside the plane's hull when the
    search applied it, as Figure 5(a) does globally (`.grid_in_hull()`, `R/hull_utils.R`);
    `bl_find_local_cf()` stores `hull_polygon` also with the hull off;
  - Pima Figure 7(a) is now a plot: the TP target's run 4 searched without the hull (pair (1,2), Age
    51 → 59.5), its cross in the white area; the caption says why there is no counterfactual
    (`03_local.R` `run4_no_hull`, `numbers.rds` `run4_no_hull`);
  - the user chose: look like Figure 5(a), no legend, captions only;
  - figures 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) changed: the Paper 1 plan must take them in.
- **Last, the final thesis in the repository** (plan archived as
  `.claude/reference/final-thesis-citations-plan.md`):
  - **The user's change:** `documentation/FINAL_ADRIAAN ROWAN_27022026.pdf` added (untracked;
    byte-identical to `../Explain XAI code/reference-materials/1 FINAL_ADRIAAN ROWAN_27022026.pdf`, which
    the notes already cited, so no page moved); the 30 September 2025 draft and
    `XAI PhD presentation - March 2026.pdf` deleted (both tracked; nothing cites the presentation).
  - **Audit:** two read-only passes checked every thesis citation against the final PDF. All 27 in
    `R/`, NEWS, README, the vignette and tests matched, and about 60 in the documents did. Fixed:
    - draft and old-location mentions (CLAUDE.md §6, the rotation note, the boundary-distance note);
    - the draft's §3.9.3 for text now in §3.10 (rotation note, 5 places; CLAUDE.md §2's section list);
    - exact quotes (rotation note pp.86, 92, 118; boundary-distance p.103; local-convex-hull p.115,
      now a paraphrase since the thesis writes "unconstraint");
    - `17_3`'s printed label "thesis (§2.3)" → "thesis (§3.5, p.55)";
    - stale links to the rotation note's old §5.3/§5.9, now Appendix A.2 (CLAUDE.md §6,
      `nearest-vertex-tie-break.md`).
  - **Added:** §3.9.2 (pp.88–90), where the thesis defines the three rules, next to p.116 in
    `set_filters()`'s roxygen, `local-convex-hull.md` §6 and the walkthrough.
  - **Kept, by the user's decision:** draft mentions in the archived plans and in
    `future-improvements-review.md` (history records). The documents name the thesis "the final
    version of 27 February 2026".
- **CLAUDE.md** (`/revise-claude-md`): the thesis PDF path and page anchors; `Rscript -e` with `|`
  fails from Git Bash; Figure 7(a)'s source; the local plot's hull trim.
- **Status:** 527 PASS / 0 FAIL / 4 WARN; `R CMD check` 0 / 0 / 0. Harness, rotation scripts and
  Pima numbers unchanged.
- **Committed and pushed** (2026-11-06, on the user's instruction): sessions 10 and 11 in one commit
  on `paper_nr1_changes`, with the new thesis PDF, the two PDF deletions, `man/dot-grid_in_hull.Rd`,
  `tests/testthat/test-hull_utils.R` and the four new `.claude/reference/` plans. The pkgdown site
  (`docs/`) is stale until the next `build_site()`.
- **Next step:**
  1. The methodology note prints the commit it was measured at (`ee64c3a` "with uncommitted
     changes"): re-run `run_all.R` and `04_method-change-effects.R`, then re-render the note, so it
     names the new commit (CLAUDE.md §5).
  2. Then session 9's "Next session" below (update the Paper 1 plan). The paper plan should also take
     in the switches (the unconstrained counterfactual is now one call), the changed local figures
     (5(b), 7(a), 7(b), 7(c), 9(a), 9(b)) and the final thesis's page numbers.

---

## Session summary (2026-11-05, session 10): `.poly_clip()` defined once

**Plan:** `.claude/plans/poly-clip-is-defined-twice-dynamic-garden.md`, archived as
`.claude/reference/poly-clip-single-definition-plan.md` (IMPLEMENTED banner).

- **Done:** removed the copy of `.poly_clip()` in `R/boundary.R`. It was dead code: R loads
  `boundary.R` before `hull_utils.R`, so the `hull_utils.R` copy always replaced it. The live function
  is identical before and after (checked with `deparse()`).
- **New test:** `tests/testthat/test-hull_utils.R` (5 blocks, 8 expectations).
- **Docs updated:** `man/dot-poly_clip.Rd` (regenerated), the implementation summary, CLAUDE.md
  (§2 archive list, §6), `review_section7_to_8.md` and NEWS.
- **Left stale by choice:** the pkgdown page `docs/reference/dot-poly_clip.html` and
  `docs/search.json`, until the next full `build_site()`.
- **Status:** 473 PASS / 0 FAIL / 4 WARN; `R CMD check` 0 / 0 / 0. Committed with session 11.

---

## Session summary (2026-11-05, session 9): every open item decided; the hull restored

**Plan:** `.claude/plans/cryptic-snuggling-bee.md`, archived as
`.claude/reference/future-improvements-hull-plan.md` (IMPLEMENTED banner).

### What was done

1. **The reply form** `documentation/methodology-changes/future-improvements-review.md`: the 9 known
   issues (A1–A9) in full and 24 other items (B1–B24) in a table, decided over five rounds with the
   user. Outcome:
   - **fix now:** A1, A3, A8, A9, B3, B9, B10, B24;
   - **future, in the paper:** A2, A4–A7, B1, B2, B6, B7, B12–B18 (four themed paragraphs in §9: the
     local plane, constraints, attribution, scope);
   - **future, design document only:** B5, B8, B11, B21–B23 (B21 and B5: no tests, the package need
     only work on the user's computer);
   - **no action:** B19, B20 (removed from CLAUDE.md §9).
2. **B10, the convex hull on each rotated plane** (the user's decision): `bl_find_local_cf(hull = TRUE)`
   builds the hull of the training data and the target per pair (`.build_hull_polygon(..., outlie =
   1)`) and filters contour vertices with `.poly_clip()` before the other filters; new fields `hull`,
   `hull_polygon`. `hull = FALSE` reproduces the earlier search exactly; the rotation scripts pin it.
   Design choices: training data + target (Appendix A), switchable, no enclosure rule.
3. **A1, the approximation gap:** `scripts/paper1-pima/05_full-space-gap.R` (nloptr SLSQP from the
   package's points and 50 random directions; geometry for the 4D hull check).
4. **Small fixes:** one colour map (`.confusion_cols`, A3); README constraint example with output
   (A8); harness title (A9); "component pair" wording (B3); centroid warning (B9).
5. **Outlier removal:** `bl_prepare_data()`'s `hull_fraction` default is now 1 (no removal); the loan
   harness and vignette keep 0.9, Pima 1 (the user's instruction).
6. **B24:** version 0.2.0; thesis citations in the rotation and boundary-distance notes moved to the
   final thesis pages and section numbers (each quote located with `pdftotext`); pkgdown rebuilt.
7. **The methodology note:** rewritten for the hull and run 7, the A1 result, the chosen paper
   sentences, the future-work decisions; then reorganised as "Changes made" and "Future work" without
   approval boxes (the user's instruction). Corrected the earlier false statement that the
   supervisor had seen a version.
8. **Documents:** `documentation/local-convex-hull.md` (new); the implementation summary; the
   Foundation document (its "no hull, intentionally" rule); CLAUDE.md; NEWS; both walkthroughs; the
   rotation note §4.8 and its hull notes; the Pima `issues.csv` (issues 5, 11, 12) and `review.qmd`.
9. **Clean-up:** `known-issues.md` removed and every reference repointed.

**Key results (Pima, Paper 1):**
- **TN target:** unchanged in every run.
- **TP target, run 1:** pair (2,3) -> (3,4); search distance 1.50 -> 1.84; counterfactual 1.42 -> 1.82.
  It now lowers every variable (Age 51 -> 41.2; without the hull it raised Age to 59.6, outside the
  training data). Runs 2, 3, 4, 5 and 7 have no counterfactual; figure 7(a) shows a "no
  counterfactual" panel.
- **Run 7 (new):** Pregnancies fixed, Glucose "decrease": the TN target has a counterfactual (0.29).
- **Over the 79 test rows:** the hull moves 11, none loses its counterfactual.
- **A1:** the package's distance is a median 1.19 times the full-space distance (quartiles 1.09 and
  1.31; range 1.00 to 2.75; 15% within 5%); TP 1.52, TN 1.05.
- **K3 (thesis group order):** the TP target now loses its counterfactual (was (2,3) -> (1,3)).
- **The loan harness:** unchanged (the hull on and off give identical results).

### Test status

- **Tests:** 465 PASS / 0 FAIL / 4 WARN / 0 SKIP (no new tests, by the user's decision).
- **`R CMD check`:** 0 / 0 / 0.
- **Pima:** `tests.R` 13/13; `04_method-change-effects.R` 6/6 self-checks (new M6 = the hull; M7
  includes run 7); `05_full-space-gap.R` both checks pass.
- **Rotation script `17_3`:** output identical to `HEAD` (timings only).

### Open points / blockers

- **No blockers.**
- **On the user's instruction:** merge into `method_developments` / `main`, tag `v0.2.0`, reinstall.
- **The Pima scripts are not under version control:** they move to the paper repo with the
  supplement (task T4); `00_setup.R` has the one `library(boundarylogic)` line to switch.
- **A fragility, not a defect:** the TP target's winning pair is narrow ((3,4) 1.8443 against (2,4)
  1.8449); another platform may pick another plane (B2's paper sentence covers it).

### Next session (exact next step)

1. **Update the paper plan**: `../Explain XAI code/reference-materials/Paper_nr1/Exploring_boundaries/.claude/plans/paper-1-plan.md`.
   Plan it first. It must take in:
   - **figures:** figure 6 is the dot chart only (no 6(b)); figure 7(a) is now a "no counterfactual"
     panel for the TP target's run 4; Table 1 gains run 7;
   - **the TP narrative:** plausibility against actionability (the issue 11 sentence chosen in the
     form), replacing the "older patient" story;
   - **§5 (local interpretation):** the hull as part of the method (B10), and the A2, B3 and B7
     sentences;
   - **§8 (properties and limitations):** remove "no hull after rotation" and "call-history
     dependence" (both resolved); add the A1 result (upper-bound sentence) and the B2 platform
     sentence;
   - **§9 (further work):** the four themed paragraphs built from the agreed sentences (G1 in the
     form);
   - **the outline's tasks T3 (paper repo infrastructure) and T4 (port the Pima scripts).**
2. **Then the hand-over to the paper writer:** the test counts, the files changed, the Pima issues
   with the user's decision on each, the review page and the methodology note, the commit SHA, the
   reinstall command, and how the Pima scripts move to the paper repo.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

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
