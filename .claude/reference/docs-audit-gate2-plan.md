# Documentation audit after the Gate 2 fixes

Date created: 2026-11-04
Date last edited: 2026-11-04

> **IMPLEMENTED 2026-11-04** on `paper_nr1_changes` (uncommitted at archive time). The
> implementation log at the end records the results.


Branch: `paper_nr1_changes`. Before any change, this plan is copied to
`.claude/plans/docs-audit-gate2.md`. That copy is the authoritative plan; this file is the
plan-mode working copy (it replaced the Gate 2 fixes plan, which is archived as
`.claude/reference/paper1-gate2-fixes-plan.md`).

## Context

The Gate 2 fixes of 2026-11-04 changed five behaviours:
- `bl_find_boundary()`'s `dist_mahalanobis`;
- the Shapley labels;
- strict `"fixed"`;
- the `"integer"` constraint, which replaces `round_to`;
- the harness `bl_pick_point()` example.

You asked for a check that the notes, the harness and the walkthroughs are still current. Two
read-only audits did that, covering every note in `documentation/`, the implementation summary,
the Foundation document, README, the vignette, CLAUDE.md, the harness and the three walkthroughs.
The walkthroughs were also checked against the last harness log.

**Results.**
- **Every console sample in the walkthroughs matches the harness log.** What is stale is prose,
  some line references, a few unlogged numbers, and older errors in two root documents.
- **Clean:**
  - `fixed-constraint-restore.md`;
  - `sparse-counterfactual.md`;
  - `boundary-distance-and-robustness.md`;
  - `scaling_axis_relabel_note.md`;
  - `nearest-vertex-tie-break.md`;
  - `biplot-rotation/README.md`;
  - `select-target-class/`;
  - `review_section4_to_6.md`.

**Your answers of 2026-11-04:**
- **No further integer examples.** The harness's Step 10 already has
  `credit_score = c("increase", "integer")`, and walkthrough Steps 10–14 show its effect.
- **The print label `Integer rounded : ` stays as it is.**

This plan corrects the rest. No behaviour changes: R code edits are roxygen and comments only.

## Files to read first

- The audit findings, recorded below per file.
- The latest harness log:
  `C:\Users\AZR2110\AppData\Local\Temp\8\claude\c--Users-AZR2110-OneDrive---Liberty-Group-Limited-Documents-9-Werk-Boundary-logic\929a5276-8265-4b5b-a864-24f8f09be099\scratchpad\harness_fix3b.log`
  (the candidate distances 0.5914 and 0.5949 are in `harness_fix3b2.log`, the same run with more
  driver prints).
- `R/local_cf.R` around `bl_find_local_cf()`'s roxygen (`blocking_constraint`) and
  `print.bl_local_result()`'s roxygen; `R/shapley.R` L211.

## 1. R source: roxygen and comments only

- **`R/local_cf.R`:**
  - The `blocking_constraint` field (about L552) also names a failed integer rounding.
  - The `print.bl_local_result()` description (about L1774) mentions the `Integer rounded :`
    line and the integer note.
- **`R/shapley.R` L211:** the comment says class 0 sorts "decreasing", but the code sorts it
  ascending. Correct the comment, not the code.
- Then `devtools::document()` and `devtools::test()`; expect 434 PASS / 0 FAIL / 4 WARN.

## 2. The harness, `development/03_loan_status_Boundary_Logic.R` (comments only)

- **Header (L7, L27–28):**
  - The SVM baseline is commented out.
  - Step 9 is "inspect predictions, choose the target row, view its global counterfactual
    (`bl_pick_point()`)".
  - Step 10 is "select the target and set the constraints".
- **L410–412:** the banner refers to Step 7b's plot arguments, but those calls take none. Point it
  at the later calls that name the four features.
- **L463–465:** add that the chosen `tdp = 13` is a true negative (predicted 0.111), not a false
  negative.
- **L471–473 and L480:** say the points are rows 13 to 10 and are not drawn
  (`plot_points = FALSE`); fix the typo "emply".
- **Step 10 comments:**
  - L514–516, `"fixed"`: a failed restore is one that no longer crosses *or* leaves
    `train_ranges`.
  - L517: `c(min, max)` is in model units (z-scores here), unlike `"integer"`, which is in raw
    units.
- **L593:** `new_title = "xx"` becomes "Sparse counterfactual, test row 13".
- **Steps 7 and 7b (L287–288, L419–420):** the commented `bl_pick_point()` lines add "(Step 9 runs
  it for the target)".
- **Harness driver (scratchpad, not the repo):** prints the numbers the walkthroughs quote that
  the log does not show:
  - the converted Step 16 applicant, and the credit-score training range (raw and z);
  - the unrounded sparse paths of Steps 14 and 15;
  - `bl_sparse_free`'s raw `loan_percent_income`;
  - `bl_shap_free$shapley_df$Contribute`;
  - the 634/635 candidate distances;
  - the `loan_percent_income` vs `loan_amnt / person_income` match count (review_section7_to_8
    L505).

## 3. The walkthroughs

**`development/reference/review_section9_to_15.md`:**
- **Headings and flow:**
  - L363: the heading "The algorithm: 10-pair loop" becomes "the pair loop (up to `max_pairs`
    = 21; 6 pairs here)".
  - L1030: the object-flow heading starts at Step 9.
  - L5: Step 16 repeats Steps 10–14 without `set_filters()`.
- **The object flow:**
  - L1057: `credit_score = c("increase", "integer")`;
  - L1092–1093: `print(bl_local)` also ends with the integer note.
- **Step 11:**
  - The table bullets (L598–603) explain the `Integer rounded :` line and the integer note.
  - Stage H points to `documentation/fixed-constraint-restore.md` §5 for the rounding.
- **L300:** the `"fixed"` constraint-table row says that a reported counterfactual holds the
  feature exactly.
- **Distances:**
  - Where 0.5838 is compared (L152, L1002), say it is the search's distance; the reported
    counterfactual's is 0.5914 (after the restore and the rounding).
  - L1002's "cost of the constraints" adds the whole-number credit score.
- **Step 12:** row 13 is a TN, so its circle is blue; the cross marks the projection of the
  restored and rounded point, which lies off the contour.
- **Missing Parameters tables** (the walkthrough convention):
  - Step 10b `set_filters()` (`tgt` and the three constraints);
  - Step 14 `plot(bl_sparse)` (show the call: `new_title`, `label_offset_var`,
    `label_offset_dist`, `ticks_var`, `ticks_n`);
  - Step 14 `bl_find_sparse_cf()`: its table keeps only `round_to` as deprecated, noting that no
    argument is set.
- **Numbers re-quoted from the driver run:**
  - L219: the converted applicant;
  - L226: the training range;
  - L228–230: the warning, quoted on one line as logged;
  - 0.360 (L954, L1026, L1118);
  - 0.0399 (L1026);
  - the Step 15 "No effect" claim (L1026);
  - 0.5914 and 0.5949 (L499, L1156).

**`development/reference/review_section7_to_8.md`:**
- **Step 7b:** mention its message "97 identified counterfactual(s) are outside the training
  feature ranges, but are retained" (log L241).
- **L390:** "about a sixth longer" is wrong. Make it "the full squared distance is 21.5% larger
  than the plotted one (11754.70 / 9675.25)".
- **L506:** Step 10 constrains `credit_score` to `c("increase", "integer")`.
- **L49–58:** the Stage 1 unpack table gains `metric_inv`, which Stage 9 now uses.
- **L505:** re-quote the 44,992-of-45,000 match from the driver, or remove it if it does not
  reproduce.

## 4. The notes and the root documents

**`documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`:**
- **Line references:** update the four moved references (L442, L452, L471, L596) to the current
  `R/local_cf.R` lines (function start L589, pair loop L649–839, filters L697–758, `dist_mah`
  L780–787, running minimum L801–838) and date them.
- **§3.5 (L624–634):** add the integer rounding after the restore; it replaces `B_x`/`B_z`/`B_pred`
  and can make `solution_found` `FALSE`.
- **Off-plane points:** L543–544 (condition 3), §4.6 (L763–773, now three counterfactuals leave
  the plane) and L1105 add the rounded point.
- **L645:** the §3.6 actionability row adds `"integer"` (applied after the search).
- **L1418:** the source list adds `.round_integers()` and `.project_local()`.

**`2 implementation_summary.txt`:**
- L5: "Last updated" 2026-11-04.
- L52–55 and L66–71: the test tree (13 files) and the `documentation/` tree (the three missing
  notes).
- L414 (`bl_counterfactual`): it is the selected point (restored or rounded); the print order adds
  `Integer rounded :` and the integer note; a failed restore shows on the Reason line.
- L559 (§6.7): the strict failure and the integer step.
- §7 status table (L592–602): the current tests; 13 files; the count taken from the real
  `devtools::test()` run.

**`1 Foundation intro documents.txt`:**
- **Today's changes:**
  - L3: the updated date;
  - L223: the valid types add `"integer"`;
  - L590–591: the Mahalanobis distance to the boundary is implemented, not future work;
  - L592–594: the deferred-tests list matches what is now tested;
  - L445–453: the test tree.
- **Older errors, corrected while here:**
  - L213: `bl_find_boundary()` reports out-of-range counterfactuals and keeps them; it does not
    constrain them;
  - L242 and L437: `bl_robustness()` is in `R/boundary_plot.R`;
  - L45–46 and L455–457: only the loan vignette exists;
  - L230: `bl_find_local_cf()`'s signature has `distance`.

**README and the vignette:**
- README L206, the vignette L454–455 and the summary L367: a restore also fails when the restored
  point leaves `train_ranges`.
- README L198–200: `"integer"` rounds in raw units when a scaling is set, unlike the other
  constraints.
- Vignette L774: Shapley bars can also be "No effect".
- **The vignette's Step 12** fixes `person_age` under strict `"fixed"`. `devtools::check()` knitted
  it today with 0 errors, so its later steps still run. Knit it into the scratchpad and check that
  its prose matches the output.

**`CLAUDE.md`:**
- L22: Phase 2's distances include `dist_mahalanobis`.
- L193, the test-coverage list: now covers `set_filters()` forms, the sparse path values and
  `bl_find_boundary()`'s `dist_mahalanobis`; trim "still untested" to what is.
- L109: the dated `shapley.R:257-258` reference.
- L150: the duplicated "See …" sentence.

## Error cases

| case | action |
|---|---|
| a driver number differs from the walkthrough | re-quote from the run; never edit by hand |
| the 44,992 match does not reproduce | remove the claim, and say so in the report |
| the vignette prose disagrees with its knitted output | correct the prose; leave the code |
| a harness comment edit changes the log beyond echoed code | stop and report (only comments change) |

## Verification

- `devtools::document()` and `devtools::test()`: 434 PASS / 0 FAIL / 4 WARN.
- **Harness run** (scratchpad driver, `pdf()` first): the log differs from `harness_fix3b.log`
  only in echoed comments and timestamps. The walkthrough numbers are re-quoted from it.
- **Knit the vignette** into the scratchpad, not `docs/`.
- `devtools::check(document = FALSE, args = "--no-manual", error_on = "never")`: 0 / 0 / 0.
- **Grep the changed files** for `round_to` as an option, "search's point is kept",
  "observed -> change", "10-pair", and `R/local_cf.R:` line references; none may be stale.
- `git status`: only the listed files change, with no `Rplots.pdf`.

## Closing the session

- Update `CLAUDE.md`: the items above, and §2 with this audit's plan if it is archived.
- Update `progress.md`: what was corrected, the test status and the exact next step (your review
  of `review.html`; Task D).
- No commit or push without your instruction.

## Implementation log (2026-11-04)

All items done as planned; no behaviour change (R edits are roxygen and one comment).

- **Tests:** 434 PASS / 0 FAIL / 4 WARN (13 files). **Harness:** the log is identical to the
  previous run apart from the echoed comments and timestamps.
- **Driver-confirmed numbers:** the converted applicant (0.368 / 0.620 / −1.203 / −5.399), the
  credit-score training range (489–754, z −3.73 to 2.67), the sparse paths (0.360), Step 15's
  0.039858 and its two "No effect" features, √d²_M 0.5838 (search) / 0.5940 (restored) / 0.5914
  (reported), the 634/635 candidates (0.5949 / 0.5914), and the `loan_percent_income` match (44,992
  of 45,000 rows).
- **One correction found:** the walkthrough said `loan_percent_income` differs from
  `round(loan_amnt / person_income, 2)` by at most 0.005; the driver gives 0.01.
- **Vignette:** knitted (into R's `tempdir()`: the scratchpad path is too long for knitr's figure
  files). Its Step 12 `person_age` restore applies (0.556), every later step runs, and its prose
  matches; Step 13 gained one sentence on the restore and the rounding.
- **Also fixed while here:** the NEWS line that still called the `round_to` questions open.
