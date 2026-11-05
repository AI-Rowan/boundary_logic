# Future Improvements Review: Replies on the Open Items

Date created: 2026-11-05
Date last edited: 2026-11-05

**Author of the replies:** Adriaan Rowan
**Prepared by:** Claude, from the package's notes as of commit `12b38e4` (branch `paper_nr1_changes`)
**Companion documents** (paths are relative to the package root):
- `documentation/known-issues.md`: the source of Part A (removed on 2026-11-05 once every item was
  decided; in git history at `c13fe44`);
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`, §4 and §6: the rotation's
  problems and further work;
- `documentation/fixed-constraint-restore.md`, §8: the follow-up searches for `"fixed"` and
  `"integer"`;
- `documentation/methodology-changes/methodology-changes.qmd`: §5 (possible future changes) and
  §6 (interpretation points);
- `CLAUDE.md`, §9: deferred work;
- `.claude/plans/cryptic-snuggling-bee.md`: the plan behind this form.

---

## Purpose

This form collects your decision on every open item, so that none is left as an "issue". Each item
ends in one of three ways:

- **Fix now:** implemented before Paper 1. Each item says whether this changes Paper 1's numbers.
- **Future improvement:** the package is unchanged, and the paper describes the item as an area of
  future improvement.
- **No action:** the item is dropped, and is not in the paper.

On 2026-11-04 you decided that none of the known issues would be implemented for Paper 1. This form
lets you revisit that item by item.

Once every item is resolved and the fixes are implemented, the methodology note
(`methodology-changes.qmd`) becomes the single record of the design decisions from the thesis to the
package. Its section on future work is then rebuilt from the items you mark "Future improvement",
with the paper sentences you agree here.

## How to reply

1. For each item in Part A, tick one decision box (change `[ ]` to `[x]`).
2. Write your recommended solution in the quote block under **Your solution**. A few words are
   enough if you agree with an option as written.
3. Where a draft paper sentence is given, tick "OK as drafted" or write your edit.
4. For Part B, write `Now`, `Future` or `Drop` in the **Your reply** column. Longer comments go in
   the block below the table.
5. Save the file and tell me.

I will then read your replies and add follow-up questions under each item's **Follow-up** heading.
The rounds continue until every item's **Status** reads "Resolved".

**Notes on the numbers.**
- All Pima numbers come from the Paper 1 example: `scripts/paper1-pima/output/table1.csv`, and
  `documentation/methodology-changes/method-change-effects.csv`.
- Phase 3 runs on the GAM refitted on four features (*Pregnancies*, *Glucose*, *Insulin*, *Age*).
- √d²_M is the Mahalanobis distance from the target.
- Line numbers are as of 2026-11-05.

**The runs in Table 1:**

| Run | Constraints |
|---|---|
| 1 | unconstrained (*Pregnancies* `"integer"` only) |
| 3 | *Pregnancies* and *Age* fixed |
| 4 | *Pregnancies* fixed, *Age* cannot decrease |
| 5 | run 1 searched on pair (1,2) only |
| 6 | *Glucose* and *Insulin* cannot increase |

Run 2 is not in Table 1.

---

## Summary

| # | Improvement area | Fixing it moves Paper 1's numbers? | Effort | My recommendation | Your decision |
|---|---|---|---|---|---|
| A1 | Measuring how far the 2D boundary lies from the full boundary | No (it adds a result) | Medium | Fix now (an experiment in the scripts) | Fix now (a) |
| A2 | Keeping every pair's plane and boundary point | No | Low to medium | Future improvement | Future (b) |
| A3 | One colour scheme for every biplot | No | Low | Fix now | Fix now (a) |
| A4 | A configurable tolerance for fixed variables | No, at the default | Low | Fix now | Future (c) |
| A5 | A follow-up search when a fixed variable cannot be restored | Yes: TN runs 3 and 4 | Medium | Future improvement | Future (b) |
| A6 | Searching in the actionable variables only | Yes: every fixed run | High | Future improvement | Future (b) |
| A7 | Letting the free variables compensate after integer rounding | Yes: TN runs 1 and 5 | Medium | Future improvement | Future (b) |
| A8 | A README constraint example that finds a counterfactual | No | Low | Fix now, before the release tag | Fix now (a) |
| A9 | The loan harness's biplot title | No | Low | Fix now | Fix now (a) |
| B1–B24 | The other open items (Part B) | see Part B | | see Part B | see Part B |

---

## Round 4 (2026-11-05): the final decisions

Every item now has a final decision. Round 3's answers were given in the chat (AskUserQuestion)
and are recorded under each item's follow-up. Round 4 changed two things:
- **No tests now.** The package only has to work on your computer, so B21 (tests of computed
  values) and B5 (the CVA identity test) are future work, in the design document only, with B23.
  B2's paper sentence notes the risk (revised in Part B's Round 4 note).
- **B10 is done now.** You decided the methodology note's §3.6 yourself: the local search keeps
  only boundary points inside the convex hull of the training data and the target on each rotated
  plane, as the thesis requires (p.117). This changes Paper 1's Phase 3 numbers. The supervisor has
  not seen the methodology note yet; the version they will receive presents §3.6 as restored.

The implementation follows the plan `.claude/plans/cryptic-snuggling-bee.md` (approved
2026-11-05).

| Step | Items | What it does | Moves Paper 1's numbers? |
|---|---|---|---|
| 1. The hull on each rotated plane | B10 | `hull = TRUE` in `bl_find_local_cf()`; `hull = FALSE` gives the earlier results | Yes: Phase 3 |
| 2. Small package fixes | A3, A8, A9, B3, B9 | One colour map; the README constraint example (run 4, with its output); the harness title; the pair wording; a warning for a target at the centroid | No |
| 3. Re-runs | | The loan harness, the Pima pipeline, the methodology effects (a new M6 for the hull) and the notes | |
| 4. The approximation-gap experiment | A1 | `scripts/paper1-pima/05_full-space-gap.R`, local search only, with the hull checks | No; it adds a result |
| 5. Release housekeeping | B24 | Version 0.2.0, the final thesis page citations, the pkgdown site; merge and tag on your instruction | No |

**Paper only, no package change:** the agreed sentences, grouped as in G1, and B2's note (the
platform statement in the paper, and `sessionInfo()` in the supplement).

**Future improvement, in the paper:** A2, A4, A5, A6, A7, B1, B2, B6, B12, B13, B14, B15, B16, B17
and B18 (B4 goes with B1). A2, B7 and B10 are described in §5 of the paper, as the method stands.

**Future improvement, in the design document only:** B5, B8, B11, B21, B22 and B23.

**No action:** B19 and B20.

### Implemented (2026-11-05)

Steps 1 to 5 are done; nothing is committed. `devtools::test()` gives 465 PASS, 0 FAIL; the Pima
`tests.R` passes 13 of 13; `devtools::check()` gives 0 errors, 0 warnings, 0 notes.

**What moved in Paper 1 (the hull, B10):**
- the TN target: nothing;
- the TP target, run 1: pair (2,3) becomes (3,4); search distance 1.50 to 1.84; counterfactual 1.42 to
  1.82. It now lowers every variable (*Age* 51 to 41.2) instead of raising *Age* to 59.6;
- the TP target, runs 4 and 5: no counterfactual now (figure 7(a) says so in its panel);
- run 7 (new): a counterfactual for the TN target with *Pregnancies* held fixed;
- the loan harness: unchanged.

**Sentences to look at again:**
- **A1:** the measured result and a filled-in sentence are under A1's follow-up.
- **A5:** Table 1's empty runs are now TP runs 3, 4, 5 and 7 and TN runs 3 and 4. The sentence itself
  holds.
- **Pima issue 11 (interpretation):** the "older patient" point no longer holds. Your choice (round
  5, in the chat): "Without the plausibility rule, the nearest counterfactual for the true positive
  raised Age from 51 to 59.6, outside the training data; inside the convex hull every counterfactual
  lowers Age, which is not actionable, so the constrained runs report none."
- **A1:** your choice (round 5) is the upper-bound sentence, recorded under A1's follow-up.
- The A7 numbers (TN, 0.14 to 0.42), B14's 27% and the B2 sentence are unchanged.

**Round 5 (2026-11-05), in the chat:**
- **The methodology note** has no approval boxes. It is organised as "Changes made" (history, core
  methodology, code improvements) and "Future work". The interpretation points sit in the sections
  they concern: issues 5 and 11 in §3.6, and the approximation gap as a third property in §3.12.
- **Outlier removal (the note's §3.11):** the code default of `bl_prepare_data(hull_fraction = )` is
  now 1, no removal. Each application keeps its own setting: 0.9 for the loan harness and vignette,
  1 for Pima. The help page explains how to use it. Tests: 465 PASS, 0 FAIL.

---

## Part A: the known issues

### A1. Measuring how far the 2D boundary lies from the full boundary

**Status:** Resolved: fix now (a), the local search only, with the A1.6 hull checks (plan 1)

**Today.**
- The searches find the decision boundary as a contour of the model's predictions on a grid in a 2D
  plane, not in the full p-dimensional space.
- The validity of a counterfactual is not in doubt: every local counterfactual is re-scored by the
  model and must cross the cutoff.
- Its distance is less certain. A counterfactual is the nearest boundary point *within the planes
  searched*, which may be further away than the nearest point of the full boundary.
- Paper 1's distances are therefore upper bounds of unknown tightness. Nobody has measured the gap.

**Where.** The global contour is at `R/biplot_grid.R:277` (used by `bl_find_boundary()`). The local
contour, built per pair in `bl_find_local_cf()`, is at `R/local_cf.R:736`.

**Effect on Paper 1.** None on the numbers. The paper can only say the distances are upper bounds.
Paper 1 also has no comparison with other methods (paper plan §9), so this is a likely reviewer
question.

**Options.**
- **(a) Measure it now, in the scripts only.** No package change.
  - For the 79 test rows of the four-feature model, find the smallest √d²_M in the full space for
    which the GAM crosses 0.5, within the training ranges. A multi-start constrained optimiser, for
    example `optim()` with a penalty or `nloptr`, would do this.
  - Compare it with the package's search point, taken before the integer rounding so that rounding
    does not confound the ratio.
  - Report the ratio of the package's distance to the full-space distance: its median, its range,
    and how often it is close to 1.
  - Moves Paper 1's numbers: no; it adds one result (a sentence or a small table in §8).
  - Effort: medium, about one session.
- **(b) Future improvement.** The paper states that the distances are upper bounds within the
  searched planes. Moves Paper 1's numbers: no. Effort: none.

**My recommendation: (a).**
- With p = 4 the full-space search is cheap.
- A measured ratio turns a caveat into a result, and goes some way to the missing comparison.
- One caution: the full-space optimum ignores the plausibility rules on the plane (the hull), so it
  can be nearer for reasons that are not a fault of the plane.

**Draft paper sentence.**
- If (b): "The counterfactual is the nearest boundary point within the planes that are searched, so
  its distance is an upper bound on the distance to the nearest point of the decision boundary in
  the full $p$-dimensional space. Measuring this gap, for example against a constrained search in
  the full space, is left for future work."
- If (a), the sentence reports the result: "On the Pima test data, the distance to the biplot
  counterfactual was [median ratio] times the distance found by a constrained search in the full
  four-dimensional space ([range])."

**Your reply**

- [ x] Fix now: option (a)
- [ ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ x] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now, as an experiment in the scripts; the package does not
change. Six design points follow. My default is the first box in each: tick it, or tick another
and add a line.

**A1.1 Which searches are measured?**
- [ ] The local and the global search (my default).
  - The local search: the four-variable GAM, with each of the 79 test rows as the target, under the
    run 1 settings.
  - The global search: Phase 2's seven-variable GAM on the one global plane, for the same 79 rows.
  - The same optimiser serves both. The global gap is likely to be larger, since one plane is
    searched instead of six, and Phase 2's importance values rest on those distances.
  - The paper sentence then gives two ratios.
- [ ] The local search only, which gives the paper's counterfactuals.

**A1.2 What the full-space point must satisfy.**
- It must cross 0.5 from the target's side, and stay within the training ranges, as the package's
  feasibility filter does. No integer rounding, and no hull.
- The package's point is taken before the integer rounding: for the TN target, the search distance
  of 0.14, not the counterfactual's 0.42. The ratio then measures the 2D approximation alone.
- The package's point is also one of the optimiser's starts, so the ratio cannot fall below 1.

- [ ] Agree
- [ ] Change (write below)

**A1.3 The optimiser.** `nloptr`, which is installed; only the script depends on it. It minimises
d²_M subject to f(x) = 0.5, with the training ranges as bounds. Its starts are the package's point
and the boundary crossings found along 50 random directions from the target, with a fixed seed.

- [ ] Agree
- [ ] Change (write below)

**A1.4 Where the result goes.**
- A new script, `scripts/paper1-pima/05_full-space-gap.R`, run after `run_all.R`.
- Its output in `output/full_space_gap.csv`, and the summary in `numbers.rds`.
- A section on `review.qmd`, one sentence or a small table in the paper's §8, and the approximation
  entry of the final design document.

- [ ] Agree
- [ ] Change (write below)

**A1.5 The paper sentence.** The draft you approved names "the biplot counterfactual", but the
comparison uses the boundary point before rounding. Revised:

> On the Pima test data, the distance to the boundary point found on the biplot planes was [median
> ratio] times the distance found by a constrained search in the full four-dimensional space
> ([range]).

The numbers go in after the run. If the result is surprising, for example a median ratio well above
1, I bring it back to you before writing it.

- [ ] OK
- [ ] Edit below

**A1.6 A plausibility check in the same script (optional; linked to B10).**
- Whether each of the package's local boundary points lies inside the convex hull of the training
  data on its own rotated plane. This puts a number on the methodology note's §3.6, which today
  says "not measured", before your supervisor decides on it.
- Whether each full-space point lies inside the convex hull of the training data in the full space
  (with the `geometry` package). The ratio is then also reported for the rows whose full-space point
  is plausible, which answers the caution under my recommendation above.

- [ ] Add both (my default)
- [ ] Add the first only
- [ ] Leave out

Your replies:

>

**Round 3 (2026-11-05, in the chat).** The local search only; otherwise my defaults.
- A1.1: the local search only. The global (Phase 2) search is not measured, so the paper sentence
  gives one ratio.
- A1.2, A1.3 and A1.4: agreed as written.
- A1.5: the revised sentence, with the numbers filled in after the run.
- A1.6: both hull checks.

**Result (2026-11-05, `scripts/paper1-pima/05_full-space-gap.R`, 79 test rows, 50 random directions
per row).** Every row has both a package point and a full-space point, and no ratio falls below 1.

| | Median ratio | Quartiles | Range | At or below 1.05 |
|---|---|---|---|---|
| With the hull (the method) | 1.19 | 1.09 to 1.31 | 1.00 to 2.75 | 15% |
| Without the hull | 1.18 | 1.08 to 1.26 | 1.00 to 2.45 | 16% |

- By outcome group, with the hull: FN 1.13, FP 1.28, TN 1.16, TP 1.24.
- The TP target: 1.84 on the planes against 1.21 in the full space (ratio 1.52; 1.24 without the
  hull). The TN target: 0.140 against 0.134 (ratio 1.05).
- 59 of the 79 full-space points lie inside the convex hull of the training data in four dimensions;
  their median ratio is 1.18, so the gap is not explained by implausible full-space points.

**The paper sentence waits for you.** The median search distance on the biplot planes is about 19%
longer than the nearest boundary point in the full space, and for 15% of the rows it is within 5%.
The revised A1.5 sentence would read: "On the Pima test data, the distance to the boundary point found
on the biplot planes was a median 1.19 times the distance found by a constrained search in the full
four-dimensional space (quartiles 1.09 and 1.31)." Please confirm or edit it.

- [ ] OK, use this sentence
- [ ] Edit below

>

**Round 5 (2026-11-05, in the chat).** You chose the upper-bound wording: "The biplot search finds a
boundary point a median 19% further from the target than the nearest boundary point in the full
four-dimensional space, so the distances it reports are upper bounds (quartiles 9% and 31%
further)." Status: Resolved.

---

### A2. Keeping every pair's plane and boundary point

**Status:** Resolved: future improvement (b), sentence (ii)

**Today.**
- `bl_find_local_cf()` searches every pair of components, but keeps the plane, grid and boundary
  point of the winning pair only. The other pairs keep only their distances.
- So another pair's plane cannot be drawn, or its counterfactual reported, without a second search
  pinned to it. Run 5 (pair (1,2) only) exists for this reason.
- The `"fixed"` restore and the `"integer"` rounding are tried on the winning pair's point only. A
  target can therefore lose a counterfactual that another pair's point would have given (see B12).

**Where.** The per-pair distances are recorded at `R/local_cf.R:684` (`all_distances_mahalanobis`).
`best_result` overwrites the plane and the point at `R/local_cf.R:865`.

**Effect on Paper 1.** None on the numbers. Figure 6 already shows the search as a dot chart of the
per-pair distances.

**Options.**
- **(a) Store each pair's plane and boundary point.**
  - Keep `Vr_rot`, `tVr_rot`, `B_x`, `B_z`, `B_pred` and the distances for each pair, in a list
    returned when `keep_pairs = TRUE` (default `FALSE`).
  - Rebuild a pair's grid on demand for plotting.
  - Moves Paper 1's numbers: no. Effort: low to medium.
  - It enables B12, and a figure of the losing planes.
- **(b) Future improvement.** Moves Paper 1's numbers: no.

**My recommendation: (b).** Nothing in Paper 1 needs the other planes. Choose (a) if you want B12
(retry the restore on the next-best pair) for Paper 1.

**Draft paper sentence.** "Only the plane of the selected pair is retained in order to be visualised if required"

**Your reply**

- [ ] Fix now: option (a)
- [X ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [x ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. Your sentence states what the package does
now. For the paper, I suggest the JCGS spelling ("visualized"). So that it also reads as future
work, I suggest a second sentence that names the improvement:

- [ ] (i) Your sentence, with the spelling only: "Only the plane of the selected pair is retained,
  so that it can be visualized if required."
- [ ] (ii) Your sentence plus the improvement (my recommendation, since it also covers B12): "Only
  the plane of the selected pair is retained, so that it can be visualized if required. Retaining
  the plane and boundary point of every pair would allow the alternative counterfactuals to be
  compared on their own planes, and the restore of a fixed variable to be retried on the next-best
  pair."
- [ ] Other (write below)

>

**Round 3 (2026-11-05, in the chat).** Sentence (ii).

---

### A3. One colour scheme for every biplot

**Status:** Resolved: fix now (a)

**Today.** `plot.bl_projection()` colours the outcome groups TP green, TN steel blue, FP orange and
FN red. Every other biplot uses the thesis scheme: TP red, TN blue, FP purple, FN orange. So red
means FN in one plot and TP in the others.

**Where.** The odd one out is `R/projection.R:320` (`col_map`). The thesis scheme is at
`R/boundary_plot.R:170`, `R/project_points.R:129` and `R/local_cf.R:1682`.

**Effect on Paper 1.** None. Paper 1 does not use `plot.bl_projection()`.

**Options.**
- **(a) One shared colour map** (a private helper) used by every plot method, on the thesis scheme.
  Moves Paper 1's numbers: no. Effort: low.
- **(b) Future improvement.**
- **(c) No action.**

**My recommendation: (a).** It is a small change, and the thesis scheme is the documented
convention.

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [ X] Fix now: option (a)
- [ ] Future improvement: option (b)
- [ ] No action: option (c)
- [ ] Other (describe below)

Your solution:

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now. No question.
- The three places that use the thesis scheme already agree exactly. They use R's named colours:
  TP `"red"`, TN `"blue"`, FP `"purple"` and FN `"orange"`.
- Without true labels, points are coloured `"blue"` or `"red"` by predicted class. A point with no
  label is `"grey50"`.
- One private helper will hold these colours, and `plot.bl_projection()` will use it instead of its
  green, steel blue, orange and red scheme.
- The green target and counterfactual markers do not change.

---

### A4. A configurable tolerance for fixed variables

**Status:** Resolved: future improvement (c), sentence agreed

**Today.**
- The local search accepts a `"fixed"` variable within ±0.5 of its observed value, in model units,
  for every variable. There is no argument to change it.
- Afterwards, the restore sets the variable back to its observed value exactly.
- The meaning of 0.5 depends on each variable's scale. In the raw-unit Pima model it is half a
  pregnancy or half a year. On the standardised loan harness it is half a training SD (±3,001 for
  `loan_amnt`).

**Where.** `R/local_cf.R:94` (`.fixed_tol <- 0.5`).

**Effect on Paper 1.**
- **Too tight leaves no candidates.** TP run 3 (*Pregnancies* and *Age* fixed) has no boundary point
  within the tolerance on any of the 6 pairs.
- **Too loose makes the restore fail more often** (A5).

**Options.**
- **(a) A `fixed_tol` argument to `bl_find_local_cf()`.** A single number, or a named vector per
  variable, in model units, with default 0.5. Moves Paper 1's numbers: no, at the default. Effort:
  low.
- **(b) As (a), with the tolerance in raw units or as a fraction of the within-group SD**
  (`sqrt(diag(W))`). The default must still equal 0.5 model units, or results move. Effort: low to
  medium.
- **(c) Future improvement.**

**My recommendation: (a).** The restore keeps every counterfactual exact whatever the tolerance, so
exposing it only gives the user control over how many candidates the search keeps.

**Draft paper sentence.**
- If (a) or (b): "The search accepts a fixed variable within a tolerance of its observed value, 0.5
  in the units of the model by default and adjustable per variable, and the reported counterfactual
  restores the observed value exactly."
- If (c): "The tolerance for a fixed variable is 0.5 in the units of the model, which is half a
  pregnancy or half a year for the Pima data. A tolerance set per variable, for example as a fraction
  of its within-group standard deviation, would suit variables measured on different scales."

**Your reply**

- [ ] Fix now: option (a)
- [ ] Fix now: option (b)
- [ X] Future improvement: option (c)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ x] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement, with the (c) sentence. No question. The
sentence's first half describes the method as it is, so it goes in the paper's description of the
constraints. Its second half goes in the future-work paragraph on constraints (see G1).

---

### A5. A follow-up search when a fixed variable cannot be restored

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.**
- `"fixed"` is strict (since 2026-11-04). When the search's point, with its fixed variables set
  back to their observed values, no longer crosses the boundary, no counterfactual is reported.
- This happens even when a counterfactual that holds the fixed variables exists.

**Where.** `R/local_cf.R:1067` (`.restore_fixed()`). Also `documentation/fixed-constraint-restore.md`
§4 and §8 (first bullet).

**Effect on Paper 1.**
- The TN target's runs 3 and 4 find boundary points, but the restored points score 0.475 and 0.487,
  so neither run has a counterfactual.
- With TP run 3 (no boundary point, A4), three of the four runs with a fixed variable are empty.
  Only TP run 4 has a counterfactual (Pima issue 12).

**Options.**
- **(a) A follow-up search.**
  - From the restored point, extend the move in the free variables along the search's own direction
    until the model crosses 0.5, within `train_ranges` and the other constraints.
  - It is a 1-D, deterministic search. Its point lies off the plane, as a restored point already
    does.
  - Moves Paper 1's numbers: **yes**. TN runs 3 and 4 may gain counterfactuals, which changes
    Table 1, figure 7, the review page and the methodology note.
  - It is a new method change, so it needs your supervisor's approval.
  - Effort: medium.
- **(b) Future improvement.** The paper reports the empty runs as a result of the strict rule, as
  agreed for issue 12.

**My recommendation: (b), with a trade-off to weigh.**
- (b) keeps Paper 1 on the method your supervisor is reviewing now.
- The cost is that three of the constrained runs in Table 1 stay empty, which weakens the
  demonstration of the constraints.
- If you choose (a), A7 (a) costs little extra, because it uses the same search.

**Draft paper sentence.** "A fixed variable is held strictly: the counterfactual must keep its
observed value, and when the restored point no longer crosses the decision boundary, no
counterfactual is reported. A follow-up search in the free variables, with the fixed variables held
at their observed values, may recover a counterfactual in these cases."

**Your reply**

- [ ] Fix now: option (a)
- [X ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. Two points to confirm:
- The sentence box is not ticked. Is the draft sentence OK?
- Table 1 then reports three runs with no counterfactual, each with its reason, as agreed for Pima
  issue 12: TP run 3 ("no boundary point"), and TN runs 3 and 4 ("the restore fails").

- [ ] Sentence OK as drafted, and Table 1 as described
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted, and Table 1 is as described.

---

### A6. Searching in the actionable variables only

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.** The rotation and the pair search use every variable, fixed ones included. The fixed
variables are held only by the tolerance during the search and by the restore afterwards. With half
the variables fixed, few points on a 2D plane hold them.

**Where.** `.bl_rotate()` is at `R/local_cf.R:45`, and the pair loop at `R/local_cf.R:688`. Also
`documentation/fixed-constraint-restore.md` §8 (second bullet).

**Effect on Paper 1.** The same as A5: with two of the four variables fixed, only TP run 4 has a
counterfactual.

**Options.**
- **(a) Hold the fixed variables before the rotation.**
  - The local biplot and its boundary are then built in the subspace of the actionable variables.
  - Every boundary point holds the fixed variables exactly, so no restore is needed.
  - Moves Paper 1's numbers: **yes**, every fixed run, through a different projection for each set
    of constraints.
  - It is a change of method, and needs your supervisor's approval.
  - Effort: high.
- **(b) Future improvement.**

**My recommendation: (b).** It is the cleanest solution, but too large a change to make during
Paper 1.

**Draft paper sentence.** "The Procrustes rotation and the search over component pairs use all the
variables, including those that are held fixed. Building the local biplot in the subspace of the
actionable variables would ensure that every boundary point keeps the fixed variables at their
observed values, at the cost of a different projection for each set of constraints."

**Your reply**

- [ ] Fix now: option (a)
- [x ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. The sentence box is not ticked.

- [ ] Sentence OK as drafted
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted.

---

### A7. Letting the free variables compensate after integer rounding

**Status:** Resolved: future improvement (b), sentence as drafted

**Today.**
- `.round_integers()` tries only the whole values just below and just above each `"integer"`
  variable, with the other variables left where the search put them.
- When the nearer whole value does not cross the boundary, the further one is taken. This happens
  even when a small change in the free variables would let the nearer one cross.

**Where.** `R/local_cf.R:1161` (`.round_integers()`). Also `documentation/fixed-constraint-restore.md`
§5.

**Effect on Paper 1.**
- The TN target's run 1 search point has *Pregnancies* 10.21 (√d²_M 0.14, prediction 0.500).
- *Pregnancies* 10 does not cross (0.487), so 11 is taken. √d²_M rises to 0.42 and the prediction
  to 0.548.
- Run 5 behaves the same way (√d²_M 0.17 to 0.44).
- The TN counterfactual's distance is therefore mostly the cost of a whole pregnancy (Pima issue
  13).

**Options.**
- **(a) Compensate.** When the nearer whole value does not cross, hold it and run the A5 follow-up
  search along the free variables. Keep the result if it is nearer than the other whole value.
  - Moves Paper 1's numbers: **yes**. The TN run 1 and run 5 distances (0.42 and 0.44) would
    probably fall.
  - It needs A5's search.
  - Effort: medium, or low if A5 (a) is chosen.
- **(b) Future improvement.** The paper explains that a count variable moves in whole steps, as
  agreed for issue 13.

**My recommendation: (b),** unless you choose A5 (a). The two then belong together.

**Draft paper sentence.** "A count variable can only change in whole steps, so the nearest valid
counterfactual may lie well past the decision boundary. For the true negative, rounding *Pregnancies*
from 10.21 to 11 triples the Mahalanobis distance, from 0.14 to 0.42; allowing the continuous
variables to compensate after the rounding may reduce this cost."

**Your reply**

- [ ] Fix now: option (a)
- [x ] Future improvement: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Draft paper sentence: [ ] OK as drafted  [ ] Edit below

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: future improvement. The sentence box is not ticked. I checked
its numbers against `table1.csv`: the distance rises from 0.1405 to 0.4191, a factor of 2.98, so
"triples" holds.

- [ ] Sentence OK as drafted
- [ ] Edit below

>

**Round 3 (2026-11-05, in the chat).** The sentence is OK as drafted.

---

### A8. A README constraint example that finds a counterfactual

**Status:** Resolved: fix now (a), the run 4 constraints with their output (plan 2)

**Today.**
- The README's "Actionability constraints" example takes the first test patient the Pima GAM
  predicts as diabetic.
- It sets *Glucose* to `"decrease"`, and *Age* and *Pregnancies* to `"fixed"`.
- `bl_find_local_cf()` returns no counterfactual: the constraints eliminate every segment on all 15
  pairs (re-run on 2026-11-04).

**Where.** The example is at `README.md:188`, and its set-up (the GAM and `bl_results_pima`) at
`README.md:136`.

**Effect on Paper 1.** None on the numbers. But the paper will point readers and reviewers to the
package, and the README is the first thing they will run.

**Options.**
- **(a) Change the example so that it finds a counterfactual.** For example one fixed variable and
  one direction, chosen by running it, with its output shown. Moves Paper 1's numbers: no. Effort:
  low.
- **(b) Keep the example**, and add a sentence saying that it finds none and why. It then shows that
  constraints can rule out every counterfactual. Effort: low.

**My recommendation: (a),** before the release tag that the paper pins.

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [a ] Fix now: option (a)
- [ ] Fix now: option (b)
- [ ] No action
- [ ] Other (describe below)

Your solution:

>

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now (a). The README's Pima example has its own six-variable
GAM, not Paper 1's model, so its results will differ from Table 1.

**A8.1 Which constraints?**
- [ ] (i) Paper 1's run 4, on the same patient (my recommendation). *Pregnancies* is `"fixed"`,
  *Age* is `"increase"` (it cannot decrease), and *Glucose* is `"decrease"`.
  - It matches the paper's constrained example.
  - If it finds no counterfactual, I drop *Glucose*, then try the next patient the model predicts
    as diabetic. I will tell you which version was used.
- [ ] (ii) The current example with one constraint relaxed. *Glucose* is `"decrease"`, *Age* is
  `"fixed"`, and *Pregnancies* is `"integer"` instead of fixed. It also shows the integer
  constraint.
- [ ] Other (write below)

**A8.2 What the README shows.**
- [ ] The printed result (the Target and Counterfactual table), and one sentence saying that strict
  constraints can leave no counterfactual, with `blocking_constraint` reporting why. This is my
  default; the sentence keeps the point of option (b).
- [ ] The printed result only
- [ ] Neither: the code only, as now

>

**Round 3 (2026-11-05, in the chat).** A8.1 (i): the run 4 constraints on the same patient, with
the fallbacks as written. A8.2: the printed result, plus the sentence on `blocking_constraint`.

---

### A9. The loan harness's biplot title

**Status:** Resolved: fix now (a)

**Today.** The first-pass result in the loan harness is titled "Loan default (prior defaulters) --
XGB, CVA biplot". But Step 3 keeps only the applicants **without** a prior default on file.

**Where.** The title is at `development/03_loan_status_Boundary_Logic.R:251` and the filter at
`:142`. The walkthrough quotes the title at `development/reference/review_section4_to_6.md:358`.

**Effect on Paper 1.** None.

**Options.**
- **(a) Fix the title** to "Loan default (no prior default) -- XGB, CVA biplot", in the script and
  in the walkthrough. Effort: low.
- **(b) No action.**

**My recommendation: (a).**

**Draft paper sentence.** None; this is not a paper item.

**Your reply**

- [x ] Fix now: option (a)
- [ ] No action: option (b)
- [ ] Other (describe below)

Comments:

>

#### Follow-up

**Round 2 (2026-11-05).** Recorded: fix now. The title becomes "Loan default (no prior default) --
XGB, CVA biplot", in the script and in the walkthrough. No question.

---

## Part B: the other open items

These items come from the rotation note (§4, §6), the methodology note (§3, Table 5),
`fixed-constraint-restore.md` (§8), CLAUDE.md (§9) and `progress.md`. I have pre-marked each with
my suggestion.

- **Your reply:** write `Now`, `Future` or `Drop` in the last column, or leave it blank to accept my
  suggestion.
- **Paper sentences:** I will draft them in the follow-up round, for the items that stay as future
  improvements and belong in the paper.

| # | Improvement area | Source | What it would do | Moves Paper 1's numbers? | My suggestion | Your reply |
|---|---|---|---|---|---|---|
| B1 | A minimal Givens rotation with a reproducible search over the planes that contain the target | Rotation note §5, §6.1; methodology note §3.12, Table 5 | Makes the local plane unique and reproducible, and searches beyond fixed pairs. The acceptance tests 4, 5 and 7 (§6.1) are incomplete. | Yes, every Phase 3 number | Future (after Paper 1) |Future (after Paper 1) |
| B2 | A cross-platform reproducibility check | Rotation note §4.1, §6.2 (1) | Re-run the Pima pipeline (`run_all.R`, `tests.R`) on Linux and macOS, for example through GitHub Actions, and compare `numbers.rds`. A check, not a change. This matters because of the plane's under-determination: a 5.0e-17 change in V moves the TP target from pair (2,3) to (1,3), and its distance from 1.42 to 1.70 (K3). | No; it tests whether they hold | Now, before submission |Future (after Paper 1) |
| B3 | Pair terminology | Rotation note §4.2 | Describe the search as planes through the target that are labelled by component pairs, not as eigenplanes. Only 21% of a selected (2,3) plane's second direction lies in the (2,3) eigenplane. | No; wording only | Now (paper wording; draft below) | Now (paper wording; draft below)|
| B4 | A frame convention for the local plot | Rotation note §4.3, §5.2 | Fixes the orientation of the axes within the local plane, which is now LAPACK's choice. | No numbers; the figures' orientation may change | Future (with B1) | Future (with B1) |
| B5 | A unit test that pins the CVA identity | Rotation note §6.2 (4) | A `testthat` test that d²_M = (n − g) · dist_z² for every pair, including two-group CVA. `17_3` §11 asserts it today; `tests/testthat/` does not. | No | Now |  Future (with B1)|
| B6 | Feasibility as a constraint inside the search | Rotation note §4.5, §5.4 | Seeks the nearest feasible point, instead of filtering contour vertices after contouring. On the loan data, `train_ranges` removes 18–56% of the candidates per pair. | Yes | Future | Future|
| B7 | Off-plane counterfactuals | Rotation note §4.6 | Marks the restored, rounded and sparse counterfactuals as off-plane on the plot, or recomputes `dist_*` for the reported point. | Only the printed distances | Future, plus a paper sentence on which distance is quoted (Table 1 reports the search's and the counterfactual's √d²_M) | Future, as recommended|
| B8 | Pair control | Rotation note §4.7 | Adds a `pairs =` argument to force any pair, and ties the `max_pairs` default to choose(p, 2). | No (p ≤ 7 already searches every pair) | Future | Future|
| B9 | A guard for a target at the centroid | Rotation note §4.8 | Checks what happens when the target is at the centre (the Procrustes matrix is zero), then warns. | No | Now | |
| B10 | A hull rebuilt on each rotated plane | Methodology note §3.6, Table 5 | Restores the thesis's plausibility rule (p.117). | Yes: some counterfactuals may move or disappear | Depends on your supervisor's decision on §3.6 | |
| B11 | A per-row feasibility flag in the global search, and the training-range filter as an option | Methodology note §3.3, Table 5 | Flags the global counterfactuals outside the training ranges (10 of 313 training and 5 of 79 test rows on the four-feature model). The option offers the thesis's rule (p.93). | No, at the default | Now |Future (with B1) |
| B12 | Retrying the fixed restore on the next-best pair | `fixed-constraint-restore.md` §8 | Rescues targets whose winning point fails the restore. | Yes, fixed runs | Future (needs A2) | Future|
| B13 | Integer rounding at other step sizes | `fixed-constraint-restore.md` §8 | A step size per variable in `set_filters()`, for example 0.5, or 100 for a loan amount. | No | Future | Future|
| B14 | A full Mahalanobis decomposition in Phase 2 | CLAUDE.md §9; methodology note §3.1 | A Shapley decomposition of the full distance. It removes the gap between the per-variable measure and the full distance: 27.4% on seven variables, 35.2% on four. | Yes: the importance values | Future |Future |
| B15 | Aggregated Shapley values | Thesis §6.2; Table 5 | Global importance built from local attributions. | No (a new output) | Future |Future |
| B16 | A CVA biplot for a singular W (p > n) | Table 5 | The generalised SVD of Ganey and Gardner-Lubbe (2026; package `wideRhino`). | No | Future | Future|
| B17 | LDA for more than two classes | CLAUDE.md §9; Table 5 | Multi-class models, with one boundary per pair of classes. | No | Future | Future, mention extension to Multi-class models being a key area of future research. Currently resolved by only targeting one of the multi classes as 1, the rest as 0, so counterfactual is out of the target class, not towards a specfic class|
| B18 | A different transformation per variable | CLAUDE.md §9; Table 5 | Display only. A non-affine transformation needs spline-calibrated axes. | No | Future | Future areas to incorporate spline axes - see Groenen, P., Le Roux, N.& Gardner-Lubbe, S. (2014). Spline-based nonlinear biplots. Advances in Data Analysis and Classification.https://doi.org/10.1007/s11634-014-0179-1|
| B19 | A non-interactive boundary-arrow layer | CLAUDE.md §9 | Arrows to the global counterfactuals for chosen rows, in one plot. | No | Future | Drop|
| B20 | Label direction per variable | CLAUDE.md §9 | A `label_dir_var` argument in `plot_biplotEZ()`. | No | Future | drop |
| B21 | Tests of computed values | CLAUDE.md §9 | Tests of the values of the local search result, the Shapley values on a real model, `bl_robustness()` and `bl_surrogate()`, and `bl_build_result()`. | No | Now, before the release | Now, before the release  |
| B22 | A vignette on the mathematical method | CLAUDE.md §9 | Derives the projection and the counterfactual geometry from first principles. | No | Future | Future|
| B23 | CRAN preparation | CLAUDE.md §9 | | No | Future |Future |
| B24 | Housekeeping before the release | `progress.md` | Rebuild the pkgdown site. Move the rotation note's citations from the 30 September 2025 draft to the final thesis pages. | No | Now, before the release | Now, before the release |

**Draft paper sentence for B3.** "Each pair of components labels one plane through the target
point; since the Procrustes rotation turns the plane, the selected plane need not coincide with the
plane spanned by the two components that name it."

**Your comments on Part B** (write the item number first, for example "B2: ..."):

>

#### Follow-up (Part B)

**Round 2 (2026-11-05).**

**Recorded without a question:**
- Future improvement: B1, B4, B6, B8, B12, B13, B14, B15, B16, B22 and B23.
- Future improvement, with your comment: B17 and B18 (sentences below).
- Future improvement, with a paper sentence on which distance is quoted: B7 (sentence below).
- No action: B19 and B20. They were removed from CLAUDE.md §9 on 2026-11-05.
- Now, before the release: B21 and B24 (one question each below).
- Now, as paper wording only: B3 (one question below).

**Questions:**

**B2. The cross-platform check (you: Future).** In K3, a change of 5.0e-17 in V moves the TP target
from pair (2,3) to pair (1,3), and its distance from 1.42 to 1.70. So a reader who re-runs the code
on macOS or Linux, or with another linear algebra library, may not reproduce Table 1. JCGS expects
reproducible code. Two practical facts:
- WSL has no Linux distribution installed on this machine, so a local Linux run is not available.
- `scripts/` is in `.gitignore`. GitHub Actions cannot run the Pima scripts until they are in a
  repository, for example the paper repository, as the supplement (task T4).

- [x ] (a) Future, with a mitigation (my recommendation if you keep Future).
  - The paper and the supplement state the platform (R 4.6.0 on Windows, with its linear algebra
    library), and the supplement includes `sessionInfo()`.
  - The paper says that the local plane is under-determined, and that the results reproduce on the
    same platform.
  - Draft sentence: "The Procrustes rotation determines the local plane only up to its second
    direction, which the linear algebra routine settles. The results are therefore reproducible on
    a given platform, but another platform may select a different plane, and hence a different
    counterfactual; the minimal rotation proposed as future work would make the plane unique."
- [ ] (b) Run the check before submission. Once the scripts are in the paper repository, one GitHub
  Actions run on Linux and macOS compares `numbers.rds`. About half a session.
- [ ] (c) Future, with no mitigation.

>

**B3. Pair terminology (Now).**
- [ x] The draft sentence above is OK as drafted
- [ ] Edit below

`?bl_find_local_cf` and the print notes speak of eigenvector pairs. Should the help pages use the
paper's wording? That is a sentence or two, with no code change.
- [ x] Yes, reword the help pages to match (my recommendation)
- [ ] No, the paper only

>

**B5. The CVA identity test (you: Future, with B1).** B5 does not depend on B1.
- The identity d²_M = (n − g) · dist_z² follows from biplotEZ's CVA scaling (V′WV = I/(n − g)).
- The per-pair rotation cancels, so the identity holds for the current rotation and for the
  minimal rotation alike.
- The test is about ten lines, and fits with the value tests of B21, which you marked Now.

- [x ] Fold into B21, now (my recommendation)
- [ ] Keep as Future, with B1

**B9. The centroid guard (blank, so Now by the form's rule).**
- It adds a warning to `.bl_rotate()` for a target at the centroid. There the Procrustes matrix is
  zero, so the plane is arbitrary.
- No result changes.
- You moved B5 and B11 to go with B1, so please confirm where B9 goes:

- [ x] Now (a warning only)
- [ ] Future, with B1

**B10. A hull on each rotated plane (blank: depends on your supervisor's decision on §3.6).**
- If the supervisor approves the departure, B10 becomes future work, with the sentence below.
- If the supervisor rejects it, B10 becomes Now. That moves Phase 3's results, so Table 1,
  figures 5 to 9 and the methodology note would be re-run.
- Until the supervisor decides, which way do you lean? A1.6 would give the supervisor a measured
  count of the package's points that fall outside the hull.

Draft sentence, if B10 stays future work: "The local search does not restrict the counterfactual
to the convex hull of the training data on the rotated plane, so a counterfactual may lie in a
region with little training data. Rebuilding the hull on each rotated plane would restore this
plausibility rule, at the cost of removing or moving some counterfactuals."

- [ x] I lean towards keeping the departure (B10 as future work)
- [ ] I lean towards restoring the hull (B10 now)
- [ ] No view until the supervisor decides
- [ ] Sentence OK as drafted, if B10 stays future work

>

**B11. A feasibility flag in the global search (you: Future, with B1).** B11 concerns Phase 2's
global search, not the rotation, so it does not depend on B1.
- On the seven-variable model that the paper uses for Phase 2, no global counterfactual falls
  outside the training ranges (0 of 313 training and 0 of 79 test rows).
- The 10 and the 5 are on the four-variable model, which the paper uses for Phase 3 only.
- So the paper needs no sentence on B11.

- [x ] Future on its own, not in the paper (my recommendation)
- [ ] Future, with B1
- [ ] Other (write below)

**B21. Tests of computed values (Now, before the release).** The plane is under-determined (B2), so
exact values of the local search may differ on another platform. Tests that pin them may then fail
on CRAN or in continuous integration.
- [x ] (a) My recommendation:
  - property tests everywhere: the counterfactual crosses the cutoff, meets every constraint and
    the training ranges, its distances recompute from its values, and the chosen pair minimises
    the distance;
  - exact values for everything that does not depend on the plane: the projection, the global
    search, the Shapley values for a given counterfactual, `bl_robustness()`, `bl_surrogate()` and
    `bl_build_result()`;
  - exact values of the local search on Windows only (`skip_on_os()`), with a comment pointing to
    B1.
- [ ] (b) Exact values everywhere.

The data: iris as two classes, as the existing tests use, plus one test on the bundled Pima file.

>

**B24. Housekeeping (Now, before the release).** It runs last, after the fixes and the supervisor's
decisions:
1. rebuild the pkgdown site;
2. move the rotation note's citations to the final thesis pages;
3. merge and tag, on your instruction.

Which version should the release that the paper pins carry? `DESCRIPTION` says 0.1.0.
- [x ] 0.2.0 (my recommendation: there are new arguments and behaviour since 0.1.0)
- [ ] 1.0.0
- [ ] Keep 0.1.0

**Draft paper sentences for the items that stay as future work.** They follow the JCGS rules. Tick
or edit each.

- **B1** (the minimal rotation): "The Procrustes rotation fixes only the direction of the target in
  the local plane, and the second direction is arbitrary. A minimal rotation through the target,
  combined with a reproducible search over the planes that contain the target, would make the plane
  unique and would search beyond fixed pairs of components."
  [x ] OK  [ ] Edit below
- **B6** (feasibility inside the search): "The constraints and the training ranges are applied as
  filters on the vertices of the contour, so the nearest feasible point may lie between two
  vertices. Treating feasibility as a constraint inside the search would find the nearest feasible
  counterfactual directly."
  [ x] OK  [ ] Edit below
- **B7** (the distances Table 1 quotes): "Table 1 reports two Mahalanobis distances: that of the
  boundary point found on the plane, and that of the reported counterfactual after the fixed
  variables are restored and the count variable is rounded. The reported counterfactual lies off
  the plane, so the local biplot shows only its projection onto the plane."
  [ ] OK  [ ] Edit below
- **B14** (the full Mahalanobis decomposition): "The per-variable importance standardizes each
  variable by its within-group standard deviation and ignores the correlations between the
  variables, so the squared per-variable measures do not sum to the squared Mahalanobis distance;
  on the seven-variable Pima model they exceed it by 27% in total. A Shapley decomposition of the
  full Mahalanobis distance would remove this gap, at a cost that grows exponentially with the
  number of variables."
  [x ] OK  [ ] Edit below
- **B15** (aggregated Shapley values): "Aggregating the local Shapley values over many observations
  would give a global measure of importance built from the local attributions."
  [x ] OK  [ ] Edit below
- **B16** (a singular W): "The canonical variate analysis biplot requires a nonsingular
  within-group covariance matrix, which fails when there are more variables than observations; the
  generalized singular value decomposition of Ganey and Gardner-Lubbe (2026) would extend the
  method to such data."
  [x ] OK  [ ] Edit below
- **B17** (multi-class models, from your comment): "The method is restricted to binary classifiers.
  A multi-class model can be explained one class at a time, by coding the class of interest as 1
  and the other classes as 0, so that the counterfactual moves the observation out of the class of
  interest rather than towards a specific other class. The canonical variate analysis biplot
  already separates more than two groups (here the four outcome groups), so extending the decision
  boundaries and the counterfactual search to multi-class models, with one boundary for each pair
  of classes, is a key area of future research."
  [ x] OK  [ ] Edit below
- **B18** (nonlinear transformations, from your comment): "The biplot axes are linear, so a variable
  can be displayed in its original units only when its transformation is affine. Spline-based
  nonlinear biplots (Groenen, le Roux, and Gardner-Lubbe 2015) would allow nonlinear
  transformations, such as a logarithm for a skewed variable, to be shown on calibrated nonlinear
  axes."
  [ x] OK  [ ] Edit below

  On the citation: the DOI you gave (10.1007/s11634-014-0179-1) was published online in 2014, but
  the journal issue is 2015. JCGS cites the issue year, so the reference reads: Groenen, P. J. F.,
  le Roux, N. J., and Gardner-Lubbe, S. (2015), "Spline-Based Nonlinear Biplots," *Advances in Data
  Analysis and Classification*, 9, 219–238.

Your edits:

>

**Round 3 (2026-11-05).** Recorded from your ticks, and from the chat for B7 and B10:
- B2: future, with the mitigation in (a); the draft sentence goes in the paper's first future-work
  paragraph (G1).
- B3: the sentence as drafted. The help pages are reworded to match (plan 2).
- B5: folded into B21, now (plan 3).
- B9: now, a warning only (plan 2).
- B10: open until your supervisor decides §3.6. You lean towards keeping the departure (future
  work). The sentence is OK as drafted, if it stays future work.
- B11: future on its own, not in the paper.
- B21: option (a) (plan 3).
- B24: version 0.2.0 (plan 4).
- The sentences for B1, B6, B7, B14, B15, B16, B17 and B18 are OK as drafted.

**Round 4 (2026-11-05).** Your instruction: "I don't want the tests. It must only work on my
computer. There can be a note in the paper that this may be a problem." Then: "add the convex hull
on the rotated plane as something that must be done now".
- **B21:** future, in the design document only. The tests are needed only if the package is used
  beyond your computer, so B21 goes with B23 (CRAN preparation).
- **B5:** future, with B21.
- **B2:** future, with this revised sentence for the paper: "The package was developed and run on
  one platform (R 4.6.0 on Windows). The Procrustes rotation determines the local plane only up to
  its second direction, which the linear algebra routine settles, so the results are reproducible
  on that platform, but another platform may select a different plane, and hence a different
  counterfactual; the minimal rotation proposed as future work would make the plane unique." The
  supplement includes `sessionInfo()`.
- **B10: now.** Your design choices (in the chat):
  - the hull on each rotated plane is built from the training data plus the target, as Appendix A
    of the thesis;
  - it is switchable: `hull = TRUE` is the default, and `hull = FALSE` gives the earlier results;
  - the rotation scripts (`17_1` to `17_3`) pin `hull = FALSE`, so the rotation note's
    measurements stay valid;
  - there is no enclosure rule for closed contours.

  The paper describes the hull in §5 as part of the method, not as future work; the draft future
  sentence above is no longer used. The methodology note's §3.6 records it as restored, decided by
  you on 2026-11-05.

  **The hull's effect, measured before the re-run (2026-11-05):**
  - over the 79 test rows (run 1), every row keeps a counterfactual: 68 are unchanged and 11 move,
    with a median distance 4.5% longer;
  - the TN target is unchanged in every run;
  - the TP target moves from pair (2,3) to pair (3,4), and its search distance from 1.50 to 1.84.
    Inside the hull, every TP counterfactual lowers *Age* from 51 to about 41, so any run that
    holds *Age* (runs 3 and 4) has no counterfactual. Run 5 (pair (1,2) only) has none either;
  - `hull = FALSE` reproduces the earlier results exactly.

  **Your decisions on this (in the chat):**
  - keep the hull on, and keep figure 7(a) on the TP target: the panel shows that no
    counterfactual exists inside the hull;
  - add a **run 7** to Table 1: *Pregnancies* fixed and *Glucose* unable to increase. For the TN
    target, a fixed variable then holds and a counterfactual exists. In the full run (2026-11-05)
    the search distance is 0.30 and the counterfactual's 0.29: *Pregnancies* is restored to 10,
    *Insulin* rises from 122 to 149.25, and *Glucose* (128.97) and *Age* (39.04) barely move. The
    TP target's run 7 has none (the restore fails).

---

## General follow-up

**G1. How the future work is grouped in the paper.** About 15 items stay in the paper as future
work. One sentence each would make a long list, so I suggest four short paragraphs in §9 (the
conclusion and further work), each built from the sentences agreed in this form:

1. **The local plane:** B1 (with B4) and B2.
2. **Constraints and the search:** A4, A5, A6, A7, B6, B12 and B13.
3. **Attribution:** B14 and B15.
4. **Scope:** B16, B17 and B18.

Some items describe the method as it is, so they sit elsewhere:
- A2, B3 and B7 go in §5 (local interpretation);
- B10 goes in §8 (properties and limitations);
- A1 is a result in §8.

- [ x] Four themed paragraphs, as above (my recommendation)
- [ ] One sentence per item
- [ ] Other (write below)

**G2. Items in the design document only, not in the paper.** These are software or documentation
items:
- A3, A8 and A9;
- B4 (folded into B1's sentence), B5, B8, B9 and B11;
- B19 to B24.

- [ x] Agree
- [ ] Move an item into the paper (write below)

>

**Round 4 (2026-11-05).** G1: B10 moves out of §8 (limitations) into §5, as part of the method.
G2: B21 joins the design-document-only items.

---

## Part C: for information, no reply needed

**Already resolved.** These items were fixed before or during the Task D list, and become design
decisions in the final document (the former `documentation/known-issues.md` §2):
- the deterministic nearest boundary vertex, which also removed the dependence on earlier calls;
- the cutoff passed to the grid;
- `rotate_deg` removed from the local plots;
- `xgboost::xgboost()` models accepted;
- empty CVA groups dropped;
- the thesis distance rule by projection (`distance = "auto"`, 2026-11-05);
- `used_in_sparse` for fixed variables;
- the stale-documentation sweep.

**Properties the paper states as limitations of the method, not as future work.** These come from
the rotation note §4.9:
- a 2D biplot is a slice of the decision surface, not a view of it;
- rotating to fit one target worsens the representation of the other points;
- the search finds the nearest boundary within the planes it draws, not in the full space (A1
  measures this);
- the CVA class means leave the plane;
- "optimal" refers only to the smallest distance.

If you think any of these should be future work instead, say so here:

>

**Round 2 (2026-11-05).** No comment, so Part C stands as written.
