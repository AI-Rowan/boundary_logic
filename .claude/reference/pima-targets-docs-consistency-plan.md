> **IMPLEMENTED 2026-11-16.** The plan that added the Pima example's FP target, the stable target
> rule, run 8 and patients C and D; corrected the documentation against the code (about 100
> findings); merged every future development into `documentation/future-improvements-review.md`;
> swept the notes to the current-code rule; and closed with the site, progress.md and memory.
> Archived from `.claude/plans/before-the-next-steps-temporal-storm.md`.

# Plan: Pima example targets and documentation consistency

Date created: 2026-11-15
Date last edited: 2026-11-16

## Context

Before the Paper 1 plan is updated (progress.md item 1), two questions had to be answered:

1. Do the Pima example's two targets illustrate the design issues well, or are there better rows?
2. Are the documentation, the code comments and CLAUDE.md consistent with the code and with each
   other, and are `gate2-issues-review.md` and `pruned-model-issues-review.md` still relevant?

The review is done (three read-only audits and a read-only scan of all 79 test rows under every
constraint set). This plan records what was found and how to correct it.

**Status (2026-11-15):** Step 0, WP1 and WP2 are done. The session-13 work is committed as `c866b27`
(not pushed). The Pima example has its third target and the stable rule; `tests.R` passes 15 of 15;
every TP and TN number is unchanged; the two review files are deleted. **WP1b is done too:** run 8,
patients C and D, `fig-rule.R` and `07_rule-effects.R` are in; `tests.R` passes 18 of 18; runs 1–7 of
the three targets, `run4_no_hull` and patients A and B are identical to before; every result matches
Appendix B. Patient C lies outside the convex hull of the training data in four dimensions, so on
pair (3,4) it lies outside the plane's hull: a property of the method, written up in
`documentation/local-convex-hull.md` §3 (user, 2026-11-15). The empty-region switch is listed in
`future-developments.md` §2.2. `04` and `06` re-run after the pin: every value identical, all seven
self-checks pass (only the commit stamp changed, so the two CSVs were restored to the committed copy).

**Step 2 done (2026-11-16): WP3, WP4, WP4b, WP6, WP7.** Five agents corrected the Foundation document,
the implementation summary, NEWS/README/vignette/walkthroughs, seven notes, and built the single
`documentation/future-improvements-review.md` (35 items; the reply form as its appendix;
`future-developments.md` removed); I did `algorithms.md`, `local-convex-hull.md`, CLAUDE.md, the R/
text (incl. 30 non-ASCII characters in comments), the pointer pass and the memory sync.
`devtools::document()` rewrote 9 Rd files; tests 578 PASS / 0 FAIL; vignette precompile passed.
`documentation/Boundary_Logic_Pima_diabetes_workflow.html` removed (assumption 1). Found for the user:
the Phase 2 full Shapley decomposition has a closed form, φ_j = v_j (W⁻¹v)_j, so the documents' cost
argument for not using it does not hold (decision unchanged; noted in item 1.6). WP5 (the history sweep, four agents) and WP8 (closing) done on 2026-11-16.

## Decisions taken (user, 2026-11-15)

| Question | Decision |
|---|---|
| Pima targets | Keep TP row 23 and TN row 25; **add FP row 26** as a third target |
| The two review files | **Delete both** (after moving the one missing item into `issues.md`) |
| Documentation scope | Correct every wrong or stale fact, fill the omissions, **and do the history sweep now** (closes progress.md item 14) |
| CLAUDE.md | Correct it and **merge duplicated rules**; dates stay |
| Future developments | **One file holds them all, with the reasons**: `documentation/future-improvements-review.md` (the list first, the dated replies as an appendix). `future-developments.md` and the copy in `methodology-changes/` are removed. CLAUDE.md only points to it |
| The three behaviours no target showed | **Range constraint:** run 8 (Glucose in 70 to 140) on the three targets. **Empty-region rule:** made-up patient C (4/95/400/35). **Training ranges:** made-up patient D (4/130/30/35) |
| Reach of that rule | CLAUDE.md and the notes: CLAUDE.md keeps one "never change `.bl_rotate()`" line (§4) and one pointer (§9); each note's future-work section becomes a pointer. The rotation note's §5 (the full derivation of the proposal) stays, and the single file points to it |

## Review results

### 1. Pima example: what the targets show

Scan: every test row (79) under runs 1–7 with the hull, and runs 1, 2 and 4 without it (790 searches,
printed only; nothing on disk). "Runs with a CF" counts Table 1's six runs (1, 3, 4, 5, 6, 7).

| Row | Group | Features (Preg/Gluc/Ins/Age) | Runs with a CF | What it shows | What it lacks |
|---|---|---|---|---|---|
| **23** (current TP) | TP | 7/181/192/51 | 2 of 6 (run 6 = run 1) | Hull moves the CF; runs 2 and 4 exist only without the hull (plausibility against actionability); the only row whose run 5 has no whole-number point | Fixed success; binding direction; "No effect"; a sparser sparse CF (it needs all 4) |
| **25** (current TN) | TN | 10/129/122/39 | 4 of 6 | Run 6 binds; run 7 fixed restore succeeds; runs 3 and 4 restore fails; rounding not to the nearest whole number (10.21 to 11); pair gain; sparse 3 of 4 | Hull effect (no TN or FN row has one); "No effect" in a plotted figure |
| **26** (new FP) | FP | 7/136/135/51, pred 0.688 | **6 of 6** | Runs 3, 4 and 7 fixed restores succeed; "No effect" in run 1 (two in runs 3 and 4); sparse 2 of 4; rounding not to the nearest (7.93 to 7); pair gain (run 5 differs); a misclassified target | Hull effect; binding direction (runs 2 and 6 equal run 1) |
| 63 (not chosen) | TP | 8/124/600/52 | 4 of 6 | Richest single row, but Insulin 600 is atypical and every approved TP number would change | |
| 72 (not chosen) | TP | 8/167/231/43 | 4 of 6 | Run 2 moves the distance 0.93 to 2.46; no hull effect | |

Across the test set: the hull moves 11 rows (8 TP, 3 FP); run 3 (both fixed) succeeds for 4 of 79 rows;
run 5 fails only for row 23.

**The target rule is fragile.** `pick_target()` (`03_local.R`) takes `which.min(|d - median|)`. In a
group of even size the two middle rows are equally far from the median; row 23 beats row 10 by
4.4e-16, so the choice rests on rounding and the `# first on ties` comment is not what happens.
The rule "the row at rank `floor(n/2) + 1` of the group's sorted distances" gives the same rows (TP 23,
TN 25, FP 26) with a gap between the two middle rows of 0.0099 (TP) and 0.0975 (FP), so it is
platform-stable.

**Still not illustrated by any target** (to be stated on the review page, not fixed): a range
constraint `c(min, max)`; the `train_ranges` filter removing a point; the empty-region rule (test row
16 only).

### 2. The two review files

Both are decision records that `issues.md` has replaced. Every decision in them is implemented or
carried in `issues.md` / `documentation/future-developments.md`; their targets (rows 68 and 73), runs
and numbers are out of date; both point to the removed `documentation/known-issues.md`. Only archived
plans mention them. One item is in neither `issues.md` nor elsewhere: Gate 2 issue 1 (no arrows in the
global figure; closed by the `bl_pick_point()` example in the harness).

### 3. Documentation

About 100 findings (Appendix A). The ones that matter most:

- **The local hull is described as "training data and the target"** in the rotation note (3 places),
  `1 Foundation intro documents.txt` §4D, `2 implementation_summary.txt` §4.8 (2 places), NEWS.md
  (2 places), the Phase 3 walkthrough's object flow, `plot.bl_local_result()`'s roxygen and the Pima
  review page. The code uses the training data alone.
- **Loan-harness numbers from before the 2026-11-06 changes** (0.547 / 634.22 / 12,398.5) in
  `fixed-constraint-restore.md`, `sparse-counterfactual.md`, `algorithms.md` C6, CLAUDE.md §6,
  `nearest-vertex-tie-break.md`. Current: 0.693 / 634.05 / 12,523.4.
- **`local-convex-hull.md` §3** says the sparse counterfactual is checked against the constraints and
  ranges; `bl_find_sparse_cf()` only re-scores.
- **The Foundation document** (which takes precedence on method questions) still describes removed
  API: eight `bl_fit_model()` types, `rounding`, `boundary =`, old step numbers, old field names.
- **NEWS.md 0.2.0** contradicts itself in six places and never announces `bl_set_scaling()`, the
  `bl_counterfactual` object or `new_title`.
- **CLAUDE.md**: a dangling "§7", 12 unlisted archive files, an incomplete R/ file map, "Filter 4"
  for filter 5, a README loan example that does not exist, rules stated up to four times.
- **Memory sync rule broken**: the two copies of `reference_data_prep_functions.md` differ.
- **The public pkgdown site** (`docs/search.json`, `docs/sitemap.xml`) still indexes CLAUDE.md and
  progress.md.

## Files to read first (at execution)

- `scripts/paper1-pima/`: `00_setup.R`, `03_local.R`, `tests.R`, `review.qmd`, `issues.md`, `fig-07.R`,
  `fig-08.R`, `fig-09.R`, `fig-ext.R` (the pattern for review-only figures), `run_all.R`.
- Each note before it is edited, plus `documentation/methodology-changes/methodology-changes.md`
  (to check that a piece of history is recorded there before it is removed from a note).
- `R/local_cf.R` (`bl_find_local_cf()`, `.restore_fixed()`, `.round_integers()`), `R/shapley.R`,
  `R/hull_utils.R`, `R/utils.R` for every statement that is corrected.
- CLAUDE.md §2 (current-code rule) and §5 (how to run the harness, the Pima example, the patches).

## Work packages

### Step 0. Starting state

- **Commit the session-13 work first** (the user's word, 2026-11-15: "commit before starting this
  plan"), so this plan's changes form their own diff:
  1. `git status` must show exactly the session-13 set: the two staged deletions
     (`methodology-changes.docx`, `.qmd`), the 12 modified files, and the three new files
     (`.claude/reference/outstanding-items-closeout-plan.md`, `documentation/future-developments.md`,
     `documentation/methodology-changes/methodology-changes.md`). Anything else: stop and report.
  2. `devtools::test()`: 578 PASS / 0 FAIL expected. A failure: stop, no commit.
  3. `git add` those files by name and commit on `paper_nr1_changes` (message: the session-13
     close-out; ends with the Co-Authored-By line).
  4. **Commit only.** No push and no fast-forward of `main` / `method_developments` (progress.md
     item 15's other half) until the user says "push".
- Then take a `git status` snapshot and copy `scripts/paper1-pima/output/numbers.rds` and `table1.csv` to
  the scratchpad (the before-state for the equality check in WP1).

### WP1. Pima example: third target, stable rule (scripts/paper1-pima, local only)

Design:

- `00_setup.R`: add `target_groups <- c("TP", "TN", "FP")`; no change to `runs_spec`.
- `03_local.R`, `pick_target(group)`:
  - add `FP = true_cls == 0L & pred_cls == 1L`;
  - replace `which.min(abs(d - med))` by: order the group's rows by `d_M_test`, take rank
    `floor(n/2) + 1`; store also `mid_gap` (the distance gap between the two middle rows; `NA` for odd
    n);
  - `targets <- lapply(setNames(nm = target_groups), pick_target)`.
  - Fix the header ("Six runs" is seven) and the run 2 comment (it equals run 1 for the TN target only).
- Table 1 (`table1.csv`), `numbers.rds`, `local.rds`: the loops already run over `names(targets)`;
  the FP block is added, nothing else changes.
- New `fig-fp.R` (review page only, as `fig-ext.R`; added to `run_all.R`): the FP target on its
  rotated biplot for run 1 and run 3 (both fixed, restore applied), its Shapley values (run 1, with
  the "No effect" bar) and its sparse counterfactual (run 1, 2 features). `ticks_v = 4L`, each PNG
  checked by eye. Whether these replace figures 8(a)/9(a) in the paper is decided in the paper plan.
- `tests.R`: new check 11 (targets are rows 23, 25, 26; each even group's `mid_gap` > 1e-6); check 8b
  runs over all targets. Checks 3–7 and 9 already loop over `names(nums$runs)`.
- `review.qmd`:
  - hull wording: "training data and the target" becomes the training data alone;
  - a section for the FP target (what it adds: fixed restores that succeed, "No effect", sparse of 2);
  - a coverage table: each design behaviour, and where it is shown (Pima run or figure, loan harness
    step, or not shown), from the matrix in Appendix B;
  - the plausibility paragraph names runs 2 and 4 only (run 3 has no counterfactual without the hull
    either).
- `issues.md`: #11 (runs 2 and 4, not 3); #12 (TP run 7's failed restore in the symptom; the FP target
  as the case where the fixed runs succeed); closed #4 and #14 marked as describing the state before
  the hull; #6 ("integer wherever it is not fixed"); #13 unchanged; new closed entry 1 (arrows).
- Stale comments: `fig-07.R` ("figure 6(b)"), `fig-09.R` ("7(a)" is 7(b)).
- `05_full-space-gap.R`, `pkg_point()`: the `in_plane_hull` diagnostic uses `rbind(Z_train_rot,
  Z_target)`; change it to `L$hull_polygon` (the package's hull). Ratios do not depend on it.
- `04`, `05` (checks) and `06` keep `tg_rows` as TP and TN: the methodology note's claims are about
  those two targets.

Error cases:

- A TP or TN number differs from the before-state: stop; the rule change must not move rows 23 and 25.
- An FP run that the scan showed as found is not found (or the reverse): stop and compare with
  Appendix B before continuing.
- `tests.R` check fails: fix the script, not the check.

### WP2. Delete the two review files

1. Add Gate 2 issue 1 to `issues.md` (closed list) (done in WP1).
2. Delete `scripts/paper1-pima/gate2-issues-review.md` and `pruned-model-issues-review.md`
   (not in git: permanent; the user chose this).
3. The archived plans that name them are verbatim archives and stay as they are; progress.md records
   the deletion.

### WP1b. Points for the three behaviours not yet shown (scripts/paper1-pima, local only)

Found by two read-only scans (the 79 test rows and a grid of 193 made-up patients, each searched with
the training ranges off and with the empty-region rule off) and a trial of eight range constraints on
the three targets. Expected results are in Appendix B.

**1. Range constraint: run 8.**

- `00_setup.R`, `runs_spec$run8`: label "Glucose in the normal range", filters
  `list(Pregnancies = "integer", Glucose = c(70, 140))`, `max_pairs = NULL`. Glucose is the 2-hour
  plasma glucose of the tolerance test, and 140 is the usual upper limit of a normal value.
- The run loop, `numbers.rds` and `local.rds` take it up as they are. `03_local.R`: `table1_cols` and
  `loc_out` gain `run8`.
- `04_method-change-effects.R` and `06_findings-effects.R` loop over `names(runs_spec)`: one line
  after `source("00_setup.R")` in each pins them to runs 1–7 (`runs_spec <- runs_spec[paste0("run",
  1:7)]`), so the methodology note's CSVs keep their rows.

**2. Empty-region rule: patient C. 3. Training ranges: patient D.**

- `00_setup.R`: a new list `rule_patients`, apart from `external_patients` (A and B stay as they are,
  with check 10): `C = 4/95/400/35` (rule `"empty_region"`) and `D = 4/130/30/35` (rule
  `"train_ranges"`), both inside every training range, each with a label and one line on what it
  shows.
- `00_setup.R`, new helper `without_empty_region_rule(expr)`: the rule has no switch in
  `set_filters()`, so the helper replaces `.encloses_data()` in the package namespace by a function
  that returns `TRUE`, evaluates `expr`, and restores the original in `on.exit()` (the method of
  CLAUDE.md §5 "to measure one rule alone"; it works with `load_all()` and with the installed package).
- `03_local.R`, new block after the external patients: each patient is searched with run 1's filters
  twice, the default and the comparison (C: inside `without_empty_region_rule()`; D:
  `train_ranges = FALSE`). Recorded per search with `extract_numbers()` (reused): pair, both distances,
  counterfactual, prediction, whether the target and the counterfactual lie in the plane's hull
  (`boundarylogic:::.grid_in_hull()`, as for A and B), and the features of the counterfactual outside
  `train_ranges`. Saved as `numbers.rds$rule_patients` and `local.rds$rule_patients`.
- New `fig-rule.R` (review page only; added to `run_all.R`): C with the rule and without it, the
  second with the training rows drawn (`plot_points = TRUE`) so the closed line with no data inside is
  visible; D with the ranges and without them. `ticks_v = 4L`, at least 2 ticks per axis, each PNG
  checked by eye.

**Tests (`tests.R`).**

- Check 5 also tests a range: the counterfactual's value lies within `[min, max]`.
- Check 12 (run 8): the TP target's run 8 has a counterfactual with Glucose in 70 to 140 that differs
  from run 1's; the TN and FP targets' run 8 equals run 1.
- Check 13 (patient C): inside the training ranges; both searches find a counterfactual; they differ;
  the one without the rule is nearer; `.encloses_data()` is the package's own function afterwards.
- Check 14 (patient D): inside the training ranges; without the ranges the counterfactual has a
  feature outside them; with the ranges it has none and is further.

**Review page and issues.**

- `review.qmd`: run 8 in `run_labels`, Table 1 and one sentence; a section "Patients that show one
  rule" (a table per patient and the four panels); the coverage table's three "Not shown" rows become
  run 8, patient C and patient D, each with a claim the render checks; the text says that the ranges
  change 6 of the 79 test rows and the empty-region rule 1 (test row 16), both read from
  `output/rule_effects.rds`, which a new script `07_rule-effects.R` writes (every test row searched
  with the ranges off and with the rule off; about 4 minutes; run after `run_all.R`, like
  `05_full-space-gap.R`, so `run_all.R` and check 8c stay at about 1 minute). The page does not quote
  a scan.
- `issues.md`: closed entry 17 (the three behaviours, and what now shows each); issue 15 updated (the
  TP target has a counterfactual in 3 of 7 Table 1 runs).

Error cases:

- Any existing number changes (runs 1–7 of the three targets, `run4_no_hull`, patients A and B): stop.
- A result differs from Appendix B, or C's two searches do not differ: stop and report; do not pick
  another point silently.
- The namespace is not restored after an error inside the helper: `on.exit()` restores it, and check
  13 tests it.

### WP3. Wrong and stale facts (tracked documents)

Source of every loan number: one fresh harness run (scratchpad driver, CLAUDE.md §5); its console log
is kept for the session. Corrections by file are in Appendix A; the patterns:

- **Hull = training data alone, plus the empty-region rule**: rotation note (intro, §3.3 stage 4, §3.6
  table), Foundation §4D, implementation summary §4.8 (step text and constraint table), NEWS.md,
  walkthrough 9–15 object flow.
- **Harness numbers**: `fixed-constraint-restore.md` §5, §6, §9; `sparse-counterfactual.md` §1, §2,
  §3.3 (the constrained sparse path is one step; the two-step path is Step 15's);
  `algorithms.md` C6; `nearest-vertex-tie-break.md` §5.1 and §2; rotation note §4.1;
  `boundary-distance-and-robustness.md` §4.2 (+15.9%: re-read from the harness log's Step 7b check).
- **`local-convex-hull.md`** §3 (sparse CF is re-scored only), §4 (run 2 added), §6 table row label
  (Step 15, not Step 10), "Appendix A7".
- **Rotation note** §3.2 and §4.7 (Paper 1's Phase 3 has 4 features and 6 pairs), §4.4 (the
  superseded "Decision" paragraph removed in the sweep, WP5).
- **Foundation document**: `bl_prepare_data()` signature and return class; four `bl_fit_model()` types
  (GAM, XGB wrapped); `rounding` replaced by `b_margin` throughout; `plot_biplotEZ()` signature;
  `bl_result` fields and the class list; workflow step numbers as the harness (Steps 1–16); §7.6 with
  all five "not to show" rules; §10.6, §11.6, §12 (point to `future-developments.md`).
- **Implementation summary**: §4.8 (hull, filter letters (d)/(e), "five-stage", m = 200), the
  "Instruction manual §3A–§3G" cross-references, §5 (`calc_gini()` use, missing helpers), §4.5 fields,
  §4.6.1 row count, §6.1, §6.10.
- **NEWS.md 0.2.0**: reconciled to one consistent statement of the released behaviour (hull, loan
  harness, `hull_fraction`, `distance`, XGB, the `hull` argument), a date on the heading.
- **Walkthrough 9–15**: Step 15 text (Step 14 is one change), Step 11's −0.0043, object flow (hull,
  `flt` fields), "unconstrained" wording, the `counterfactual` column name. Walkthrough 4–6: the
  pointer to CLAUDE.md §9.
- **Vignette `.Rmd.orig`**: "random forest" is a decision tree; "Step 5b"; the `"fixed"` range case.
  Then `Rscript vignettes/precompile.R`.
- **The reply form** (`methodology-changes/future-improvements-review.md`): its stale rows are handled
  in WP4b, where it becomes the appendix of the single future-developments file.

### WP4. Omissions

- `algorithms.md`: C4 signature (`verbose`, the `distance` choices) and missing output fields (C4, B1,
  B3, C7, C8, A4); C2 (`.check_switch()`, reserved names `hull` / `train_ranges`); §2 inventory
  (`.encloses_data()`); S3 argument names; C4 note against the roxygen on a target outside the data;
  C7 (the PhD zero rule by class).
- `scaling_axis_relabel_note.md`: `print.bl_counterfactual()`, `print.bl_local_result()`.
- `fixed-constraint-restore.md`: `all_distances_x`, `dist_x`; the §4 sentence on the TP target.
- `boundary-distance-and-robustness.md` §2: the Pima 4-feature off-plane share.
- `biplot-rotation/README.md`: the scripts also pin `set_filters(tgt, hull = FALSE)`.
- Future developments not yet listed anywhere (added to the single file in WP4b): F1 (the
  `bl_find_boundary()` stop when every contour lies on one side), F5, F7; `new_title` for
  `plot.bl_projection()`; `URL` / `BugReports` in DESCRIPTION; a help page for `pima_diabetes.csv`;
  the uncovered display arguments.
- README.md: `bl_apply_scaling()` in the function reference; the `@v0.2.0` install line and the tag;
  `bl_project_points()` / `bl_predict()` filed correctly.
- NEWS.md 0.2.0: `bl_set_scaling()` / `bl_apply_scaling()`, `bl_counterfactual`, `new_title`.
- "Component pair" for "eigenvector pair" in the three notes and the Foundation document.

### WP4b. One file for every future development and its reason

Target: `documentation/future-improvements-review.md`. Sources: `documentation/future-developments.md`
(162 lines) and `documentation/methodology-changes/future-improvements-review.md` (1,200
lines, tracked).

Layout of the merged file:

- **Header**: title, `Date created` / `Date last edited`, one paragraph: this is the only place for
  future developments; none is implemented or part of Paper 1; each needs an explicit instruction.
- **Part 1, the list** (from `future-developments.md` §1–§4: method; code and interface; tests;
  documentation and release). Each item: what it is; **the reason it is future work and not done now**
  (the decision and its reason, taken from the reply on that item); the code it touches; its effect on
  results; effort; the reply it comes from (A1–A9, B1–B24) or its other source.
  - Item 1.1 takes in CLAUDE.md §9's `.bl_rotate()` paragraph (the rank-1 cross-product, the
    measurements, why it is not covered by the "outputs proven equivalent" rule) and points to the
    rotation note §5 for the derivation.
  - New items from WP4 (F1, F5, F7, `new_title`, DESCRIPTION fields, the Pima help page, the display
    arguments).
  - Whatever a note's future-work section holds that the list lacks.
- **Part 2, appendix: the replies of 2026-11-05**, kept as the dated record of the decisions. Changes
  only: a status line naming what later work superseded (training-only hull, the switch in
  `set_filters()`, the B2 result, TP run 4); the three dangling companion references repointed; the
  13 line-number citations replaced by function names.

Steps:

1. `git mv documentation/methodology-changes/future-improvements-review.md
   documentation/future-improvements-review.md` (the rename stays staged; CLAUDE.md §5 on `git mv`).
2. Build Part 1 above the replies; move the reasons in; add the new items.
3. Check by script that every reply decided "future improvement" is cited by a Part 1 item, and that
   every Part 1 item has a reason.
4. `git rm documentation/future-developments.md` (tracked since the Step 0 commit), after confirming
   that every item and sentence of it is in the merged file.
5. Notes reduced to a pointer: `fixed-constraint-restore.md` §8, `boundary-distance-and-robustness.md`
   §7, the rotation note §6 (its §5 stays; its §4, the current problems, stays), the Foundation
   document §12, and any other future-work remark found by grep in `documentation/`, the root
   documents and the walkthroughs. A fact about the current code stays where it is; only the proposal
   moves.
6. CLAUDE.md (with WP6): §9 is the rule and the pointer, nothing else; §4 gains one line ("never
   change `.bl_rotate()` or implement its proposed replacement without an explicit instruction");
   §2 lists the one file and drops the two old entries; the rotation-note and fixed-restore entries in
   §2 and the Shapley-decomposition remark in §6 point to the file; the current-code rule names the
   file's appendix as the third place where history may stand.
7. Repoint every reference to the two old paths (grep: CLAUDE.md, progress.md, the notes, the root
   documents, the walkthroughs, `scripts/paper1-pima/issues.md` and `review.qmd`, the memory files).
   The archives in `.claude/reference/` are verbatim and keep the old paths.

Error cases: an item in a note's future section that is not in the merged file is moved before the
section is cut; a reply with no clear reason is listed for the user, not given an invented one.

### WP5. History sweep (current-code rule, CLAUDE.md §2)

For each note: describe the current code in the present tense; an alternative in the conditional, with
its numbers measured by removing that one rule (namespace patch in a scratchpad driver) or left out;
before a piece of history is removed, check it is in the methodology note or NEWS.md, and add a line to
NEWS.md if it is in neither (the supervisor-approved methodology note is edited only where a statement
would otherwise be wrong). Thesis citations stay; research-code citations stay only in the rotation
note §2.

| Note | Work |
|---|---|
| `nearest-vertex-tie-break.md` | Rewritten: the rule, why exact, what random ties would do (measured by patching `.nearest_idx_block()`); the diagnosis and before/after tables go |
| `sparse-counterfactual.md` | §1–§2 as current behaviour; §3.1 (`round_to` history) removed |
| `fixed-constraint-restore.md` | Revision header, "What changed", the Gate 2 7-feature example |
| `boundary-distance-and-robustness.md` | The dated markers and the "Changed 2026-10-31" paragraph |
| `select-target-class/select_target_class_note.md` | From a bug record to the current class rule of `bl_select_target()` |
| `biplot_rotation_and_pair_selection.md` §3–§4, `biplot-rotation/README.md` | About 42 dates, "since"/"until"; §4.4's superseded decision; §2 untouched |
| `scaling_axis_relabel_note.md` | Light: two markers |
| `2 implementation_summary.txt` | "since 2026-…", "until", commit hashes |
| Walkthroughs 7–8 and 9–15 | Three history lines |
| R/ code comments | "(since 2026-11-06)", "decided 2026-11-06", "removed on 2026-11-06" |

The rotation note is the largest; it is done last and may need its own session.

### WP6. CLAUDE.md

- **Corrections**: §2 ("three places" is four; the "§7" of `local-convex-hull.md`; the walkthrough
  table: Steps 9–16 and the missing source files); §5 (check date; pkgdown bullet: two pages, and
  strip the entries from `search.json` and `sitemap.xml`; `new_title` scope: every biplot plot method
  except `plot.bl_projection()`); §6 (README has no loan example; 634.05; "Filter 5"; `hull = FALSE`
  reproduces the earlier search under CVA, under PCA only with `distance = "mahalanobis"`); §7 step 5
  (no table); §10 (three branches, the tag, `origin/Github_deployment`, `origin/original_PhD_code`).
- **Omissions**: §2 lists the 12 unlisted archive files, `documentation/select-target-class/`,
  `biplot-rotation/README.md`; §3 file map adds `datasets.R` and the missing helpers
  (`.resolve_distance()`, `.confusion_cols`, `.encloses_data()`, `.grid_in_hull()`,
  `.points_in_polygon()`, `.scale_to_raw()`, `plot.bl_projection()`); §5 raw-unit bullet adds
  `.fixed_text()` / `.integer_text()`; §5 Pima bullet updated for three targets and the stable rule.
- **Future developments**: none described in CLAUDE.md; §9 and one §4 line only (WP4b step 6).
- **Merged duplicates** (one full statement, pointers elsewhere): the "never show" rules (§8 keeps
  them, §4 points); the local hull rule (§4 keeps the rule, §6 the mechanics); `distance = "auto"`
  (§5 keeps it); the display rules (pointer to the scaling note only); the R-invocation and Windows
  path-limit bullets (one group each).

### WP7. R/ source: text only, no behaviour change

- `plot.bl_local_result()` roxygen ("and the target"); `bl_find_sparse_cf()` title (greedy, not
  "zeroing"); `result.R` step numbers; the console note's file name in `plot.bl_boundary()`;
  `loan_data` documentation in `datasets.R` (a CSV in `inst/extdata`, not a data object); the harness
  comment `library(boundary_logic)`.
- ASCII only (CLAUDE.md §4). Then `devtools::document()` and `devtools::test()`.
- The version stays 0.2.0 and the tag stays where it is; HEAD differs from the tag in text only.

### WP8. Closing

- **Memory**: make `reference_data_prep_functions.md` identical in `.claude/reference/` and the memory
  folder (keep the `bl_set_scaling()` section, add `bl_apply_scaling()`, 10 fields, drop the pointer
  to a script that does not exist); update `project_boundarylogic.md` and `reference_docs.md`.
- **pkgdown**: rebuild `docs/` by the CLAUDE.md §5 recipe (NEWS and README changed), with the
  post-build step extended to remove the CLAUDE and progress entries from `search.json` and
  `sitemap.xml`.
- **progress.md**: item 14 closed; item 1 gains the FP target (Table 1 layout, which FP figures go in
  the paper, the target rule's wording); items 11 and 13 corrected; the deletion of the review files.
- **Hand-over Claude Doc**: the third target and the rule (after the commit, with its SHA).
- Archive this plan to `.claude/reference/` with an IMPLEMENTED banner; list it in CLAUDE.md §2.
- Run `/revise-claude-md`.

## Order

1. Step 0, WP1, WP2 (Pima). **Done.**
1b. WP1b (run 8, patients C and D). **Done.**
2. WP3, WP4, WP4b, WP6, WP7 (facts, omissions, the single future-developments file, CLAUDE.md,
   roxygen). Report the test count, wait.
3. WP5 (sweep), the rotation note last.
4. WP8.

Apart from the Step 0 commit, no commit or push at any point without the user's word.

## Assumptions to confirm on approval

1. `documentation/Boundary_Logic_Pima_diabetes_workflow.html` (a May 2026 export of a vignette that no
   longer exists) is removed with `git rm` and its listings dropped. Recoverable from git.
2. `plot.bl_projection()` does not get `new_title` now (a behaviour change after the tag Paper 1
   pins); it is listed as a future development and CLAUDE.md states the rule's scope.
3. In the merged `future-improvements-review.md` the replies of 2026-11-05 stay a dated record
   (appendix) and are not rewritten; the methodology note's §5 "Future work" (supervisor-approved)
   stays, with one sentence pointing to the single file.
4. DESCRIPTION is not edited (its date and the missing `URL` go to the next version).
5. The methodology note is edited only where the third target makes a sentence wrong.
6. (Replaced by WP1b.) The empty-region rule is switched off for patient C's comparison by patching
   the package namespace inside the example script; a switch in `set_filters()` would be a package
   change after the tag, so it is listed as a future development.
7. Run 8's bounds (70 and 140) stay as they are; no further check (user, 2026-11-15).
8. `progress.md` item 1 (the Paper 1 plan) decides whether run 8 and patients C and D go in the paper;
   here they are in the scripts, Table 1's CSV and the review page.

## Verification

- **Pima**: from `scripts/paper1-pima/`, `Rscript run_all.R`, `Rscript tests.R` (all checks pass,
  including new check 11); `all.equal()` of `numbers.rds$runs$TP`, `$TN`, `$run4_no_hull`,
  `$external` against the scratchpad copy; the FP runs against Appendix B; `Rscript
  05_full-space-gap.R` (ratios identical to the stored `full_space_gap.csv`); `quarto render
  review.qmd`; each new PNG checked by eye.
- **WP1b**: `Rscript run_all.R`, then `Rscript tests.R` (18 of 18: checks 12, 13 and 14 are new),
  then `Rscript 07_rule-effects.R` (6 test rows changed by the ranges, 1 by the empty-region rule);
  `identical()` of runs 1–7 of the three targets, `run4_no_hull` and `external` against a copy of
  `numbers.rds` taken before WP1b; run 8 and patients C and D against Appendix B; `quarto render
  review.qmd` (its coverage claims pass); the four new PNGs checked by eye; `git status` unchanged
  (the folder is gitignored).
- **Package**: `devtools::document()`, `devtools::test()` (578 PASS expected, no change);
  `Rscript vignettes/precompile.R` (the `*-claims` chunks pass); `devtools::check()` 0/0/0.
- **Numbers**: every loan number edited is grepped in the fresh harness log; every Pima number in
  `numbers.rds`, rounded as `review.qmd` does; every conditional number has a scratchpad driver.
- **Text**: grep the repo for "and the target" (hull context), `634.22`, `0.547`, `12,398`,
  `rounding` (Foundation), "eight supported", `known-issues`, "§7" of `local-convex-hull.md`,
  `R/[a-z_]+\.R:[0-9]+`, "eigenvector pair", and "since 2026" / "until 2026" in `documentation/`
  outside the methodology folder: no hit left that is not intended.
- **Future developments**: `future-developments.md` and `methodology-changes/future-improvements-review.md`
  have no hit outside `.claude/reference/`; in CLAUDE.md "future" appears only in the §2 listing, the
  §4 line and §9; the script check of WP4b step 3 passes.
- **Sync**: `diff` of the two `reference_data_prep_functions.md` copies is empty.
- **Site**: `grep -c "CLAUDE\|progress" docs/search.json docs/sitemap.xml` gives 0.
- `git status` against the Step 0 snapshot: only intended files changed; nothing untracked removed.

## Appendix A. Findings by file

Severity: W wrong, S stale, O omission, D duplicate, C cosmetic, H history wording.

**documentation/**

| File | Findings |
|---|---|
| `algorithms.md` | S: C6 "634.22". O: C4 signature and outputs; B1, B3, C7, C8, A4 outputs; C2 switch checks; §2 inventory; S3 argument names. C: C4 note against roxygen; C7 PhD zero rule incomplete |
| `local-convex-hull.md` | W: §3 sparse CF "checked against the constraints". O: §4 run 2. C: "Appendix A" / "A7"; §6 row label |
| `fixed-constraint-restore.md` | S: §5, §6, §9 harness numbers. O: `all_distances_x`, `dist_x`. C: §4 "instead"; "eigenvector pair". H: moderate |
| `sparse-counterfactual.md` | S: §1, §2, §3.3 path and table. C: PhD zero rule; "line 102". H: heavy |
| `nearest-vertex-tie-break.md` | S: §5.1 (6813.80, 14638.39, 0.368), §2. C: §8 dated section numbers; research-code line. H: whole note |
| `boundary-distance-and-robustness.md` | S: §4.2 +15.9%. O: §2 Pima share. C: "eigenvector pair". H: moderate |
| `scaling_axis_relabel_note.md` | O: two print methods. C: "CLAUDE.md Section 9". H: light |
| `future-developments.md` | O: F1, F5, F7. Merged into `future-improvements-review.md` and removed (WP4b) |
| `select_target_class_note.md` | H: heavy (bug record) |
| `biplot-rotation/` note | W: hull "and the target" (intro, §3.3, §3.6), "no enclosure rule". S: §4.4 superseded decision; §3.2, §4.7 (p = 7, 21 pairs); §4.1 (0.547 to 0.509). C: §3.3 order wording, the [0,1] clamp. H: heavy in §3–§4 |
| `biplot-rotation/README.md` | O: `hull = FALSE` pin. C: "eigenvector pair". H: light |
| `future-improvements-review.md` | S: B2 row, Round 4 summary, A5/A6. C: three dangling companion references; 13 line citations. Moved to `documentation/` as the single file (WP4b) |
| `Boundary_Logic_Pima_diabetes_workflow.html` | S: export of a removed vignette |

**Root and package documents**

| File | Findings |
|---|---|
| CLAUDE.md | W: §2 "§7"; §5 `new_title` scope; §6 README loan example, "Filter 4", `hull = FALSE` under PCA. S: §5 check date; §6 634.22. O: §2 (12 archive files, 3 documentation entries, walkthrough table); §3 file map; §5 raw-unit list; §5 pkgdown; §10. D: four rule groups. C: "three places"; §7 "table" |
| `1 Foundation intro documents.txt` | W: §4A signatures and class; §4B eight types; `rounding` throughout; §4D hull; §4H `plot_biplotEZ()`; §5 fields and classes; §10.6. S: §3 step numbers; §11.6; §12. O: §7.6; §9. C: §1.3, §4E, "eigenvector pairs" |
| `2 implementation_summary.txt` | W: §4.8 hull (twice), filter letter. S: "four-stage", 100x100, §3A–§3G references, §6.1, §6.10, §5 `calc_gini()`. O: §5 helpers; §4.5 fields; §4.4.1. C: §4.6.1, §6.11. H: throughout |
| NEWS.md | W: six internal contradictions in 0.2.0. O: three unannounced features; no date |
| README.md | O: `@v0.2.0`, `bl_apply_scaling()`. C: function grouping |
| DESCRIPTION, `.Rbuildignore` | S: date. O: `URL`, `BugReports`. C: a redundant and a missing ignore line (not edited; assumption 4) |
| Vignette `.Rmd.orig` | W: "random forest". C: "Step 5b". O: `"fixed"` range case |
| Walkthrough 9–15 | W: Step 15 text; object flow hull. S: −0.0016. C: column name, `flt` fields, "unconstrained". H: two lines |
| Walkthrough 7–8, 4–6 | H: one line. O: pointer to CLAUDE.md §9 |
| Harness | C: `library(boundary_logic)` comment |
| R/ roxygen and comments | S: `plot.bl_local_result()`, `bl_find_sparse_cf()` title. O: `loan_data`. C: `result.R` steps, console note. H: three comments |
| progress.md | S: items 11 and 13 overstated; a script that does not exist; a scratchpad copy |
| Memory folder | W: `reference_data_prep_functions.md` copies differ. S: 9 fields; `project_boundarylogic.md` header and count |
| `docs/` | O: `search.json` and `sitemap.xml` index CLAUDE.md and progress.md |

**scripts/paper1-pima/ (local)**

| File | Findings |
|---|---|
| `03_local.R` | W: tie-break by rounding; "Six runs"; run 2 comment |
| `review.qmd` | W: hull "and the target"; run 3 in the plausibility paragraph |
| `issues.md` | W: #11 run 3. O: #12 TP run 7; Gate 2 issue 1. S: closed #4, #14, #6 wording |
| `fig-07.R`, `fig-09.R` | S: figure references in comments |
| `05_full-space-gap.R` | C: `in_plane_hull` diagnostic includes the target |
| The two review files | S: superseded (deleted, WP2) |

## Appendix B. Scan results for the three targets

Search distance is √d²_M at the search point; CF is Pregnancies / Glucose / Insulin / Age.

| Target | Run | Result | Pair | Distance | CF | Notes |
|---|---|---|---|---|---|---|
| FP 26 (7/136/135/51) | 1 | found | (2,3) | 0.539 | 7 / 130.95 / 121.66 / 54.97 | 7.93 rounded to 7; 1 "No effect"; sparse 2 |
| | 2 | = run 1 | | | | |
| | 3 | found | (2,3) | 0.591 | 7 / 123.47 / 120.77 / 51 | restore applied (7.22, 50.53); 2 "No effect" |
| | 4 | found | (2,3) | 0.595 | 7 / 123.27 / 119.60 / 51.04 | restore applied (7.33) |
| | 5 | found | (1,2) | 0.639 | 8 / 137.79 / 113.19 / 55.48 | differs from run 1; 2 Contradicts |
| | 6 | = run 1 | | | | |
| | 7 | found | (2,3) | 0.591 | 7 / 123.47 / 120.77 / 50.53 | restore applied (7.22) |
| | 1, 2, 4 without the hull | = with the hull | | | | the hull does not move this target |
| TP 23, TN 25 | all | as `numbers.rds` | | | | unchanged by this plan |

Behaviour coverage with the three targets: hull effect (TP runs 1, 2, 4; patients A, B); binding
direction (TN run 6); fixed success (TN run 7; FP runs 3, 4, 7); fixed failure (TN runs 3, 4; TP run 7);
rounding not to the nearest (TN runs 1, 5; FP run 1); rounding with no valid value (TP run 5); pair gain
(TN, FP run 5); "No effect" (FP run 1); sparse sparser than full (TN 3, FP 2); misclassified target
(FP); outside the training data (patients A, B). The three that were not shown are covered by WP1b:

**Expected results of WP1b** (from the scans of 2026-11-15; search distance / distance to the
counterfactual; CF is Pregnancies / Glucose / Insulin / Age).

| Point | Search | Pair | Distances | CF | Prediction |
|---|---|---|---|---|---|
| TP target, run 8 (Glucose 70 to 140) | default | (2,3) | 2.239 / 2.257 | 5 / 139.61 / 135.31 / 41.12 | 0.475 |
| TP target, run 1 (for comparison) | default | (3,4) | 1.844 / 1.819 | 6 / 148.91 / 106.10 / 41.20 | 0.481 |
| TN and FP targets, run 8 | default | | = run 1 | | |
| Patient C (4/95/400/35, prediction 0.019) | with the rule | (3,4) | 1.442 / 1.605 | 7 / 118.71 / 467.41 / 42.58 | 0.546 |
| | empty-region rule off | (1,3) | 1.387 / 1.399 | 3 / 106.70 / 312.44 / 32.33 | 0.513 |
| Patient D (4/130/30/35, prediction 0.376) | with the ranges | (1,4) | 0.841 / 0.886 | 3 / 131.71 / 111.72 / 34.08 | 0.514 |
| | training ranges off | (3,4) | 0.569 / 0.494 | 4 / 121.36 / −15.41 / 33.15 | 0.502 |

In the TP target's run 8 the search point has Pregnancies 5.91 and the whole value taken is 5.

Among the 79 test rows (run 1's filters): the training ranges change 6 (rows 11, 18, 20, 34, 64 and
68; rows 20 and 64 would have Age 20.92, below the youngest training patient at 21), and the
empty-region rule changes row 16 only. Not chosen: test
rows 16 and 64 (real, but small effects, and row 16 is also outside the training range), patients
16/180/300/50 and 0/170/400/60.
