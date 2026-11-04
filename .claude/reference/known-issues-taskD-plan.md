> **IMPLEMENTED 2026-11-04.** Archived from `.claude/plans/i-have-reviewed-the-robust-rossum.md`
> (the approved plan, with its implementation log). Part A (the Pima scripts) is local only;
> Part B is `documentation/known-issues.md` and one CLAUDE.md §9 bullet.

# Pima decisions after the pruning, and Task D (the known-issues note)

Date created: 2026-11-04
Date last edited: 2026-11-04

## Context

You answered `scripts/paper1-pima/pruned-model-issues-review.md`:

| # | Issue | Your decision |
|---|---|---|
| 4 | Figure 7 cannot show runs 1, 4 and 6 of one target | (a) three panels from both targets |
| 5 | Run 2 binds for neither target | (b) drop run 2 from Table 1 |
| 9 | The Phase 2 Mahalanobis check is -27.4% | (b) also report the 4-feature model's check |
| 11 | The TP counterfactual raises Age and Pregnancies | (a) paper text |
| 12 | Most "fixed" runs have no counterfactual | (a) paper text; (b) and (c) recorded as future issues |
| 13 | The integer rounding triples the TN distance | (a) paper text; (b) recorded as a future issue |

**None of the package-level options is implemented for the paper.** They are recorded as future
issues in Task D's note, `documentation/known-issues.md`, together with the other Task D issues.

**Your choices for the note (2026-11-04):**
- **Open items in full:** the Task D items still open, the new Pima items, and the two items you
  parked on 29 Oct.
- **Resolved items:** the Task D items fixed since the list was written go in a short "resolved"
  table, with their commits.
- **CLAUDE.md §9:** gets one bullet pointing to the note. It keeps its own future-work items,
  which the note links to rather than repeats.
- **The arrow layer is left out.** CLAUDE.md §9 keeps its existing note.

**What changes:** Part A edits the Pima scripts (gitignored, local). Part B adds one new
documentation file and one CLAUDE.md bullet. There are **no `R/` changes**, so `devtools::test()`
is not re-run (434 PASS stands).

## Files to read first

- `scripts/paper1-pima/`: `fig-07.R`, `03_local.R`, `02b_prune-refit.R`, `review.qmd` and
  `issues.csv` (all read in this session).
- `.claude/plans/boundary-logic-fix-prompt.md`, Task D (the format of the note and its list).
- `documentation/fixed-constraint-restore.md` §8, for the wording of options (c) and (d).
- A recent note's header style, e.g. the top of `documentation/fixed-constraint-restore.md`.
- `progress.md`, for the parked items: the README's Pima example and the harness title.

## Part A: Pima scripts (local)

### Issue 4(a): Figure 7 from both targets
- **`fig-07.R`:**
  - Drop the target-choice rule. The panels are fixed:
    - 7(a): TP run 4 (its run 1 is Figure 5(b));
    - 7(b): TN run 1;
    - 7(c): TN run 6.
  - The titles name the target and the constraint. The ticks and offsets stay per panel:
    - 7(a) is pair (1,2), like Figure 6(b), so it reuses 6(b)'s Pregnancies offset;
    - 7(c) keeps its extra ticks for Age, Glucose and Insulin.
  - A panel whose run has no counterfactual is reported in the log, as now.
- **`review.qmd`:**
  - Remove `t7` and `both_bind()`. Keep `binds()`, which issue 5 uses.
  - Give each panel its own caption, e.g. "Figure 7(a). TP target, Pregnancies fixed and Age
    cannot decrease (run 4); its run 1 is Figure 5(b)."
  - Replace the paragraph about the choice of target with one line: each panel is on its run's own
    winning plane.

### Issue 5(b): run 2 leaves Table 1
- **`03_local.R`:** `table1_cols` becomes `c("run1", "run3", "run4", "run6")`. Run 2 stays in
  `numbers.rds` and in the local-runs table.
- **`review.qmd`:**
  - Drop run 2 from Table 1's `cols`.
  - The text under the Table 1 heading says "runs 1, 3, 4 and 6". A sentence, generated from the
    data, says that run 2 is left out because it equals run 1 for both targets (both unconstrained
    counterfactuals already raise Pregnancies and Age). If run 2 ever binds, the sentence says so
    instead.

### Issue 9(b): the Mahalanobis check for the 4-feature model
- **`02b_prune-refit.R`:**
  - Run `plot(bl_bnd_train, type = "boxplot")` on the 4-feature training-row search, the same
    call as `02_global.R`.
  - Store `mahalanobis_check` and `robustness` in `reduced_info`.
- **`03_local.R`:** `numbers$reduced` carries both.
- **`review.qmd`:** in "Pruning to four features", add a two-row table of the check for the 7- and
  the 4-feature model: plotted, full, and difference (% of full).

### `issues.csv`
- **4:** closed. Fixed (a): Figure 7 shows TP run 4, TN run 1 and TN run 6.
- **5:** open, paper text. (b) is done in the scripts; Task E writes the sentence.
- **9:** open, paper text. (b) is done: the 4-feature check is on the review page, and its number
  is quoted after the run. Task E writes one sentence on both checks.
- **11:** open. Paper text (a), Task E.
- **12:** open. Paper text (a), Task E; options (b) and (c) are recorded in
  `documentation/known-issues.md`.
- **13:** open. Paper text (a), Task E; option (b) is recorded in
  `documentation/known-issues.md`.
- **Review page:** the heading becomes "Open issues (paper text for Task E)".

## Part B: Task D, `documentation/known-issues.md` (new)

**Format of each open entry:** symptom, file:line, impact on results, and suggested fix (Task D's
format). Every file:line is taken from the current code with grep when the note is written.

**Open entries:**

| # | Issue | Source |
|---|---|---|
| 1 | The approximation-gap experiment has not been done: how far the 2D boundary lies from the p-dimensional one | Task D 3 |
| 2 | Only the winning pair's plane is stored (`Vr_rot` and the others); the other pairs keep only `all_distances*`, so Figure 6(b) needed a pinned run | Task D 4 |
| 3 | `plot.bl_projection()` uses a different confusion colour scheme from the other biplots | Task D 5 |
| 4 | The ±0.5 tolerance for "fixed" is hard-coded (`.fixed_tol`, in model units) | Task D 6 |
| 5 | A follow-up search with the fixed features held exactly: Pima TP run 3 and TN runs 3 and 4 have no counterfactual | Pima 12(b); `fixed-constraint-restore.md` §8 (c) |
| 6 | A search in the free features only | Pima 12(c); §8 (d) |
| 7 | After integer rounding, let the other features compensate: the TN distance goes from 0.14 to 0.42 | Pima 13(b) |
| 8 | The README's Pima constraint example finds no counterfactual for its first predicted-diabetic patient | parked 29 Oct |
| 9 | The harness plot title says "Loan default (prior defaulters)" | parked 29 Oct |

**Resolved since the list was written** (a table: item, how it was fixed, commit):
- `best_pair` instability and call-history dependence: the exact nearest vertex (`f0e3dd4`). The
  plane under-determination that remains is in CLAUDE.md §9 and the rotation note §4.1, and the
  note refers to them rather than repeating them.
- The cutoff not passed to `bl_build_grid()`: `f4f3647`.
- `rotate_deg` not applied to the sparse cross: `rotate_deg` was removed from the local plots
  (`f4f3647`).
- Stale documentation: the session-1 sweep (`f2a076b`) and the 2026-11-04 audit. Before it is
  marked resolved, I re-check by grep that `dist_z`, the sparse colours, the Shapley `seed` and the
  "stratified" split are now correct.

**Other future work:** a short section that points to CLAUDE.md §9 and the rotation note §6
(the Givens rotation, the PCA selector default, the full-Mahalanobis Shapley, per-variable scaling,
LDA multi-class). The arrow layer is not mentioned.

**CLAUDE.md §9:** one new bullet: "Known issues: `documentation/known-issues.md` records the open
defects and limitations (Task D, 2026-11-04); read it before starting future work."

## After implementation

- Append a dated entry to `.claude/plans/paper1-pima-example.md` (issue decisions 4, 5 and 9) and
  a Task D line to the implementation log of `.claude/plans/commit-and-push-the-mutable-hearth.md`.
- Archive this plan as `.claude/reference/known-issues-taskD-plan.md` with an IMPLEMENTED banner
  (CLAUDE.md §2: a documentation change in the package).
- At session close, update `progress.md`: Tasks C and D done, Task E next. Your filled-in form is
  not edited.

## Error cases

| case | action |
|---|---|
| a Figure 7 panel's run has no counterfactual after the re-run | not drawn; the caption gives the reason (as now) |
| run 2 binds for a target after the re-run | the Table 1 sentence says so; I report it to you |
| a "resolved" Task D item turns out still open on re-check | record it as open, not resolved |
| a file:line moves before the note is committed | re-grep at the end; line numbers are taken last |

## Verification

1. From `scripts/paper1-pima/`, run `run_all.R`, then `tests.R` (13 checks), then render
   `review.qmd`.
2. Look at fig-07a/b/c: every axis shows at least 2 ticks, and no names overlap.
3. Check the rendered page:
   - Table 1 has no run 2 column, and the run-2 sentence is there;
   - the Pruning section shows both Mahalanobis checks;
   - the open issues are 5, 9, 11, 12 and 13.
4. **`known-issues.md`:**
   - every file:line opens on the code it names (grep);
   - every commit hash exists (`git show --stat <hash>`);
   - every Pima number matches `numbers.rds`.
5. `git status`: the only tracked changes are the new note and the CLAUDE.md bullet. Nothing is
   committed.

## Implementation log (2026-11-04)

- **Part A done as planned.** Figure 7 is TP run 4, TN run 1 and TN run 6; run 2 is out of Table 1
  (the review page says why, from the data); the 4-feature Mahalanobis check is on the review page.
  `tests.R`: 13 of 13 PASS (8a-8c identical). `review.html` re-rendered; the open issues are 5, 9,
  11, 12 and 13, all paper text for Task E.
- **The 4-feature Mahalanobis check is -35.2%** (plotted 335.28, full 247.91), against -27.4% on
  the 7-feature model: the pruning widened the gap. Quoted in issue 9.
- **Figure 7(a) shows the same point as Figure 6(b).** The TP target's run 4 counterfactual is
  identical to run 5's: both are on pair (1,2), and restoring Pregnancies to 7 (run 4) gives the
  same point as rounding it to 7 (run 5). Reported to you; no change made.
- **Part B done as planned.** `documentation/known-issues.md`: 9 open entries, 5 resolved items
  (`f0e3dd4`, `f4f3647`, `f2a076b`), and a pointer to CLAUDE.md §9. Every file:line was printed
  and checked; one (`best_result`, 815 -> 819) was corrected. The README constraint example was
  re-run: still no counterfactual (all 15 pairs eliminated by the constraints).
- **Not in the plan:** `documentation/fixed-constraint-restore.md` §4 quoted the Pima example's
  7-feature run as the current one. Two sentences now say it is the 7-feature Gate 2 run, and that
  on the 4-feature model the TN target's restores fail too (`known-issues.md` §1.5).
- **Added to the resolved table** (from the working plan's Task D section, "the planning
  findings"): `xgboost()` models and empty CVA groups (`f4f3647`), `used_in_sparse` for fixed
  features (`7aa1344`), and the walkthrough's stale `person_age = "fixed"` examples (`7aa1344`).
  The note has 9 open entries and 8 resolved ones.
