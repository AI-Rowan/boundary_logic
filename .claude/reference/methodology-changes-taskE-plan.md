# Task E: the methodology-changes note for the supervisor

Date created: 2026-11-04
Date last edited: 2026-11-04

> **IMPLEMENTED 2026-11-04 (session 8), drafts at Gate 3.** Approved plan, verbatim below the
> implementation notes. Deliverables: `documentation/methodology-changes/methodology-changes.qmd`
> (rendered to `.docx`) and `method-change-effects.csv`, written by
> `scripts/paper1-pima/04_method-change-effects.R` (local) at commit `540f9ad`.
>
> **Implementation notes (changes from the plan):**
> - **Configurations share `bl_results`.** Each configuration builds one fresh `bl_results` and shares
>   it across its runs, as `03_local.R` does. The searches are deterministic, and self-check 1 confirms
>   runs 1-6 reproduce `numbers.rds`.
> - **Self-check 4 is not circular.** It checks that the derived CVA labels, given in the package's
>   level order, reproduce the baseline V exactly; comparing group counts would only compare the labels
>   with themselves.
> - **M3 is narrower than planned.** `bl_find_boundary()` already applies the hull clip and the
>   enclosed-area rule (`R/boundary.R:263-285`). The only difference from the thesis text is that it
>   keeps counterfactuals outside the training ranges.
> - **2.12(a) is corrected.** On the 4-feature model p - K = 1, so component 4 is fixed up to its sign;
>   the completion is arbitrary only when p - K >= 2 (the 7-feature model, the thesis's pair (1,7)).
> - **Interpretation-point numbers are now in the CSV:** the targets' observed values, issue 5's
>   run 2 = run 1, issue 11's GAM edf and the over-50 count. This keeps the note free of hand-typed
>   numbers.
> - **A `claims` chunk** in the `.qmd` asserts every qualitative statement, so the note does not render
>   if one fails.
> - **K3 result:** the thesis level order changes V by 5.0e-17, yet the TP's pair moves (2,3) -> (1,3)
>   and its distance 1.42 -> 1.70. It is cited in §3.12(b) as evidence of the plane's
>   under-determination.
> - **Test status:** `devtools::test()` gives 434 PASS / 0 FAIL / 4 WARN (no `R/` change).
> - **2026-11-05, at the user's request:** §3.3 states the recorded reason for the post-selection
>   range check.
>   - The reason: `ea6a12f`, the PCA rank-2 back-projection, under which a pre-filter picked a farther
>     point.
>   - The commit's "CVA was unaffected" is refuted on Pima.
>   - The section adds the distances-against-feasibility choice, and new CSV rows saying whether the
>     targets' global CFs are in range.
>   - The same reason is now in `review_section7_to_8.md` Stage 8 and the implementation summary
>     §6.11.


## Context

Paper 1 (JCGS) relies on the `boundarylogic` package. Before the paper relies on any change made to the
method since the thesis (Stellenbosch, 27 Feb 2026), the supervisor must approve it. Task E writes a
note the supervisor can comment on in Word. Tasks A–D are done and Gate 2 is closed. Gate 3 follows
this task.

Specification: `.claude/plans/boundary-logic-fix-prompt.md` (Task E, Parts 1–2) and the working plan
`.claude/plans/commit-and-push-the-mutable-hearth.md` (Task E section). **No `R/` code changes.**

**Decisions taken in planning (2026-11-04):**
- **Structure (user, at plan review):**
  - The main focus is a separate section that compares the **core methodology** of the thesis and the
    package. It leaves out the coding improvements used to implement it.
  - The order:
    1. the **history of changes**;
    2. the **core methodology** comparison;
    3. the **code improvements**, summarised: the implementation changes, the coding issues found and
       the bug fixes.

    This replaces the spec's one-line bug-fix appendix.
  - Decision boxes go on each core change. The code-improvements section gets one box for the section.
- **`max_pairs`:** measured on the 4-feature model only. It has 6 pairs, so 10 and 21 both search all of
  them. The pair search's gain is all pairs against pair (1,2) only. The note explains that more
  variables give more pairs than a cap may cover, citing the thesis's 7-feature Pima pair (1,7), which
  is 16th in the package's order.
- **Pair selector:** state the CVA identity and confirm it on Pima. PCA differs; that stays deferred.
- **Paper-text items:**
  - Pima issues 9, 12 and 13 go into the "effect on Pima" of their core changes.
  - Issues 5 and 11 go into a short section of interpretation points, each with a proposed sentence and
    a decision box.

## Planning findings that shape the note

1. **Three baselines.**
   - The thesis **text** (T) defines the thesis version.
   - The thesis's **Appendix A** (A, pp.158–218) prints its code.
   - The branch `origin/original_PhD_code` (B) was committed on **2026-03-09**, after submission
     ("Update starting code for boundary logic", the `*_optimized.R` scripts). B is not always A.

   Report A and B separately where they differ. Two examples:
   - the tie-break: A uses `which.min`, B uses `max.col()` with random ties;
   - the data behind the local hull (verify).
2. **Cite the final thesis.**
   - Printed page = PDF page − 16 in `../Explain XAI code/reference-materials/1 FINAL_ADRIAAN
     ROWAN_27022026.pdf`, with the final section numbers.
   - The `thesis-md` chapters are fine for searching, but their section numbers are one lower from
     Chapter 3 on, and they omit Appendix A.
   - The rotation note and CLAUDE.md §6 cite the **30 Sept 2025 draft**, about 20 pages higher. The note
     says so in one line.
3. **The thesis does not state the (n−g) identity.**
   - Its W is the within-group SSP with no divisor (§3.4 p.48; p.118). So M′WM = I and d²_M = dist_z²,
     with a factor of 1. The factor n−g comes from the package's W/(n−g).
   - The selection rule is on p.86 (§3.9.1); the correspondence is stated qualitatively on p.91.
4. **Ties.** The working plan says "the thesis's searches used the same random tie-break". That is too
   strong:
   - T calls the method deterministic (p.92);
   - A uses `which.min`;
   - only B's post-thesis optimized script uses the random `max.col()`.

## Files to read first

- **Spec:** `.claude/plans/boundary-logic-fix-prompt.md` §Task E; the working plan's Task E section.
- **Thesis:** the final PDF (T and A), and `thesis-md` 07–11 for searching.
- **Branch code:** read only with `git show origin/original_PhD_code:"scripts/<file>"`. The scripts are
  `1 Functions for biplot.R`, `1.1`, `1.2`, `1.3`, `1.4`, `3.1`, `3.2`, `4.1`, `4.3` and `4.4`.
- **Package history:**
  - `git log --format="%h %ad %s" --date=short method_developments paper_nr1_changes`;
  - `NEWS.md`;
  - `2 implementation_summary.txt` §4.2.1, §4.6.1, §4.8, §6;
  - the notes in `documentation/`:
    - `boundary-distance-and-robustness.md`
    - `fixed-constraint-restore.md`
    - `sparse-counterfactual.md`
    - `nearest-vertex-tie-break.md` §5
    - `known-issues.md`
    - `biplot-rotation/biplot_rotation_and_pair_selection.md` §3.4, §4.1, §5, §6
- **Pima scripts:**
  - `scripts/paper1-pima/00_setup.R` (`run_local()`, `extract_numbers()`, `mahal()`,
    `in_train_ranges()`, `global_summary()`, `runs_spec`);
  - `02_global.R`, `02b_prune-refit.R`, `03_local.R` and `review.qmd`, as patterns;
  - `issues.csv`: issues 5, 9, 11, 12 and 13.
- **Toggles:**
  - `plot.bl_boundary(distance =)`: `R/boundary_plot.R:102`
  - `bl_find_local_cf(max_pairs =, distance =)`: `R/local_cf.R:590`
  - `bl_build_result(cva_classes =, outlie =, b_margin =)`: `R/result.R:291`

## Content of the note

### Section 1: History of changes

A dated timeline table: date | source/commit | what changed | category (core / code / fix). Built from
`git log`, `NEWS.md` and the plan archives. The milestones:

| Date | Source | What changed |
|---|---|---|
| 27 Feb 2026 | Thesis (T and A) | The baseline |
| 2026-01-15 | Package repo created | |
| 2026-03-09 | `original_PhD_code` (B) | Post-submission optimized scripts, including the random tie-break |
| 2026-03-10 to 03-24 | Package v1 port | `hull_fraction` 0.9 (`b8d4965`); local search with `max_pairs` 10 and the `.bl_rotate()` porting bug (`4058730`); `outlie` 1 (`f030a64`); `"fixed"` ±0.5 (`bb96efe`); 3-dp predictions (`fe47d75`, `dffbeb8`); the global `train_ranges` check moved after selection (`ea6a12f`) |
| 2026-05-20/23 | Mahalanobis migration (`271b71e`) | W/(n−g); within-class SD in Phase 2; the Mahalanobis selector; `b_margin`; the `.bl_rotate()` fix |
| 2026-05-29 | Accuracy fixes (`1719c73`) | The grid clamp; the Shapley seed |
| 2026-10-29/30 | Paper 1 (`f2a076b`, `7aa1344`, `f0e3dd4`, `f4f3647`) | Tasks A and B; "No effect"; the greedy sparse CF; the `"fixed"` restore; the deterministic tie-break; the Gate 2 code fixes |
| 2026-10-31 | `49c9067` | `max_pairs` 21; the Phase 2 Mahalanobis check |
| 2026-11-04 | `67a1424` | Strict `"fixed"`; `"integer"`; `dist_mahalanobis` |

### Section 2: Core methodology, thesis against package (the main focus)

**2.0 Side-by-side table of the whole method,** step by step: Phase 1 (data, projection, metric,
boundary), Phase 2 (global counterfactual, importance), Phase 3 (rotation, pair search, local
feasibility, constraints, Shapley, sparse CF). Columns: step | thesis (T, page) | package (function) |
changed? It marks the steps that are unchanged as well as those that changed, so the supervisor sees the
whole method once.

**2.1–2.11, one entry per core change.** Each entry has:
- the thesis version (T first, then A/B where they differ);
- the package version (the function);
- the reason;
- the effect on Pima (measured, or "not measured");
- a decision box: "Decision: ☐ Approve ☐ Query ☐ Reject" and "Comment:".

| # | Core change | Thesis → package | Effect on Pima |
|---|---|---|---|
| M1 | Phase 2 variable importance | **T:** T_j = Σ\|d_ij/std_j\| with the overall SD (§4.2.2 p.103; R = ΣT_j, §4.2.3 p.105). Note p.3's different wording. **→** within-class SD √diag(W) (`plot.bl_boundary()`) | **Toggle** `distance = "euclidean"`, 7 and 4 features: totals, ranks, and whether the thesis rule prunes the same three. **Issue 9** (the Mahalanobis check, −27.4% and −35.2%). |
| M2 | The metric W | **T:** within-group SSP, no divisor, four outcome groups (p.48, p.92, p.118). **→** W/(n−g); also used in Phase 2 (M1); under PCA, W from the 0/1 class | Exact: the package's d²_M is (n−g) times the thesis value, and rankings are unchanged. **Stored:** n, g. |
| M3 | Feasibility of the global counterfactual | **T:** boundaries filtered to the training bounds and the hull, empty enclosed areas discarded (§3.10 step 4 p.93; p.125; p.101 looser). **B:** drivers pass `train_ranges = NA` and remove out-of-range test rows first. **→** counterfactuals outside the ranges are kept and reported; test rows are not removed. **Verify** the package's enclosed-area handling. | **Stored:** counterfactuals outside the ranges (train and test; 7 and 4 features), and the test rows that T/B would remove |
| M4 | The pair search space | **T:** all p(p−1)/2 pairs (§4.3.1 p.117; §6.3 p.144). Pima selected (1,7) (p.118, Fig. 4.11); the loans TP example (2,4) (p.137). **→** capped at `max_pairs = 21L` (all pairs for p ≤ 7; 10L until 2026-10-31) | **Toggle** on 4 features: 1, 10 and 21 (10 = 21 = all 6). The text covers larger p. |
| M5 | The pair selector | **T:** Euclidean under PCA, Mahalanobis under CVA, in full dimensions (§3.9.1 p.86; pp.117–118). **→** Mahalanobis for both by default. Under CVA this is identical; under PCA it differs (deferred). | **Toggle** `distance = "euclidean"`, run 1 and unconstrained: the same pair, and the ratio is n−g for every pair. Report whether the integer step (off-plane) also matches. |
| M6 | Feasibility of the local counterfactual | **T:** the counterfactual must lie within the convex hull of the training data (p.117); dropped in the unconstrained mode (p.115). **A/B:** a hull rebuilt in each rotated plane; closed contours must enclose the target. **→** no hull and no enclosure rule; only `train_ranges` and `set_filters` | Not switchable. Not measured. **This departs from T.** The package's reason does not cover a hull rebuilt per plane, so flag it as a likely query. |
| M7 | Immutable features (`"fixed"`) | **T:** two options (§3.9.2 p.88). The examples filter the boundary points and hold Pregnancies through the sparse Shapley step (p.122); p.136's wording conflicts. **→** a ±0.5 search tolerance, then a strict restore | **Stored** (`fixed_check`), runs 3 and 4. **Issue 12.** |
| M8 | Whole-number features (`"integer"`) | **T:** treated as continuous (p.87); whole numbers are not reachable (p.119). **→** rounding in raw units, with a re-check | **Stored** (`integer_check`), runs 1 and 5. **Issue 13.** |
| M9 | Shapley classification | **T:** strict inequalities; a no-impact variable can be ignored (pp.92, 119–120, 123, 139). **A/B:** a zero is "Supports" (class 0) or "Contradicts" (class 1). **→** "No effect" for \|φ\| < 1e-8 (this matches T's intent) | **Observed:** No-effect features per run; would the sparse CF differ? |
| M10 | The sparse counterfactual | **T:** keep the contradicting variables and change the supporting ones (p.120, p.122); no stopping rule. **A/B:** chosen by hand. **→** greedy by \|φ\|, stopping at the first crossing | **Computed:** run 1, greedy against all-"Supports" (`bl_predict()`) |
| M11 | Outlier removal before the projection | **T:** Pima none (p.99); loans a 90% hull in the biplot (p.127). **→** `hull_fraction = 0.9` default, on a 2-PC stand-in, before the model fit | Paper 1 uses `hull_fraction = 1`, as the thesis's Pima example did. Not measured. |

**2.12 Two properties the thesis does not state** (unchanged; one decision box, since the paper relies on
them):
- **(a) Components beyond K.** Four-group CVA has K = 3 non-zero eigenvalues (p.49). Columns K+1..p are
  an arbitrary completion in both. The thesis's Pima pair (1,7) uses one. **Stored:** whether each Pima
  `best_pair` uses a component > 3.
- **(b) The rotation is not unique.** C = Y′Z has rank 1, so the plane's second direction is a numerical
  tie-break. The thesis's §3.5 (pp.52–54) does not discuss uniqueness. This is inherited, and is
  documented in the rotation note §4.1.

### Section 3: Code improvements (summarised)

One table: item | thesis code (A/B) → package | why | effect on Pima. **One decision box for the
section.**

- **K1 Contour band and prediction rounding.** ε 0.01 with 2 dp (A/B) → 0.001 (as in T's Iris example,
  p.78) with 3 dp. **Partial toggle:** `b_margin = 0.01`.
- **K2 The global grid hull.** `outlie` 0.9 (B) → 1. **Toggle.**
- **K3 The CVA group level order.** TN, FN, FP, TP (B) → FN, FP, TN, TP, pinned. **Toggle:** `cva_classes`.
  It can move the plane, through 2.12(b).
- **K4 The inverse of W.** `solve()` → `chol2inv` with a ridge fallback. **Stored:** whether the ridge
  fired.
- **K5 The grid clamp to [0, 1].** **Observed:** whether the warning fires.
- **K6 The nearest-vertex tie-break.** B random → exact nearest; matches T (p.92) and A. Not measured;
  the loan-harness scale is in `nearest-vertex-tie-break.md` §5.
- **K7 The local plane's data and grid extent.** Test data and plot bounds → training data, range ±10%.
  Not measured.
- **K8 The hull for the global test search.** Test-data hull (B) → training hull (as T says). Not
  measured.
- **K9 The constraint syntax and pair order.** Relative forms; pairs ordered by the second index. These
  matter only for ties.
- **K10 New diagnostics.** `dist_mahalanobis` in `bl_find_boundary()`, and the Phase 2 Mahalanobis
  check.
- **K11 Minor items:**
  - grid m;
  - the cutoff source;
  - the empty-boundary fallback;
  - the filter order;
  - the Shapley seed (p > 14);
  - `X_center`/`X_sd` from the filtered data;
  - pairs with no boundary point skipped (verify `1.4:100-113`, `4.3:178`);
  - the PCA local search allowed;
  - empty CVA groups dropped;
  - the Shapley step rounding (2 → 3 dp).
- **Bug fixes, one line each:**
  - the `.bl_rotate()` porting bug (`271b71e`; it restores T/A);
  - **Task A:** the graphics-device restore;
  - **Task B:** `class_obs`;
  - **Gate 2:**
    - `used_in_sparse` for fixed features;
    - `rotate_deg` removed from the local plots;
    - the grid-cutoff check;
    - `xgboost()` models;
    - empty CVA groups named;
    - `bl_apply_scaling()` and the target-range warning;
    - the Shapley labels.
- **Open coding issues:** one line pointing to `known-issues.md` (9 open). The ones with a method
  bearing appear in Section 4.

### Section 4: Possible future changes (spec Part 2)

A summary table (item | why | expected effect | effort | recommendation), then a short paragraph each:
- the minimal Givens rotation with a plane-family search (rotation note §5–§6), which addresses 2.12(b);
- storing every pair's plane;
- the approximation-gap experiment;
- the full-Mahalanobis Shapley in Phase 2;
- aggregated Shapley values;
- a configurable "fixed" tolerance;
- a CVA biplot for a singular W (Ganey & Gardner-Lubbe 2026, *SMA* 35, 825–847, `wideRhino`);
- LDA for more than two classes;
- per-variable standardisation;
- the deferred Task C issues (the follow-up searches for "fixed" and "integer", `known-issues.md`).

**Found in planning:**
- a hull rebuilt in each rotated plane (M6);
- the PCA selector default (M5);
- T's global range and enclosed-area filter (M3).

### Section 5: Interpretation points for the paper

- **Issue 5:** a constraint that does not bind leaves the counterfactual unchanged.
- **Issue 11:** the TP counterfactual ages the patient. It describes the model, not the clinical
  mechanism; the drop in risk after age 55 rests on 19 of 313 training patients over 50.

Each has a proposed sentence and a decision box.

### Front matter, Summary and Appendix

- **Before Section 1:**
  - **Summary:** one paragraph, and the list of core changes that depart from T itself (M3, M6, and M5
    under PCA), not only from its code.
  - **How to read this note:**
    - the three baselines, and the text wins;
    - final-thesis pages;
    - the Pima set-up: 392 rows, 7 features then 4, the GAM, CVA on four outcome groups, seed 121, the
      targets;
    - the commit, from the CSV meta rows;
    - decision boxes.
- **Appendix:** reproduction. The commit, `Rscript run_all.R`, `Rscript 04_method-change-effects.R`,
  and the render command.

**Excluded (spec):** API renames, plotting, the list of models, exact Shapley over 2^p subsets, and
two-group CVA.

## Deliverables

### 1. `scripts/paper1-pima/04_method-change-effects.R` (local, gitignored)

- **Input:** `source("00_setup.R")`, then the outputs of `run_all.R` (`data-model.rds`, `global.rds`,
  `reduced.rds`, `numbers.rds`). Targets: `numbers$targets$TP$row` and `$TN$row`.
- **Output:** `Sys.getenv("METHOD_CHANGES_OUT", "../../documentation/methodology-changes")`.
  `method-change-effects.csv` is in long format:
  - columns: `change` (M1–M11, P2.12, K1–K5), `model` (`7f`/`4f`), `target`, `setting`, `quantity`,
    `value`, `source` (`toggle` / `stored` / `observed` / `computed`), `note`;
  - plus `change = "meta"` rows: version, git SHA, `git_dirty`, R and mgcv versions, time.
- **`measure_config(label, args)`** builds a fresh 7f and 4f `bl_results` with `bl_build_result(bl_dat,
  bl_mod, method = "CVA", <args>)`. It records:
  - the 7f training importance totals and the three weakest;
  - the 4f test global search, and the targets the Task C rule picks;
  - run 1 for the paper's two targets: `best_pair`, √d²_M, `B_x`, `B_pred` and `solution_found`.

  It runs for four configurations:
  - **baseline:** must reproduce `numbers.rds`, which is asserted;
  - `b_margin = 0.01` (K1);
  - `outlie = 0.9` (K2);
  - the thesis group order (K3): `cva_classes = factor(<training confusion from bl_predict()>, levels =
    c("TN","FN","FP","TP"))`.
- **Single toggles on the baseline:**
  - M1: `plot(bl_bnd_train, type = "boxplot", distance = "euclidean")` against the default, 7f and 4f;
  - M4: run 1 with `max_pairs` 1, 10 and 21;
  - M5: `distance = "euclidean"`, run 1 and unconstrained, with the per-pair ratio
    `all_distances_mahalanobis / all_distances^2` (min, max) and n−g.
- **Stored, observed and computed:**
  - M2: n, g, `metric_type`;
  - K4: ridge (`tryCatch(chol(metric))`);
  - K5: warnings around `bl_build_result()`;
  - M3: `global_summary()`, `in_train_ranges(test_data)`;
  - M7: runs 3 and 4, `restore_impact`;
  - M8: runs 1 and 5, `integer_impact`;
  - M9: No-effect counts, runs 1–6;
  - M10: greedy against all-Supports;
  - 2.12(a): `best_pair` > 3, every run.

  Every local run uses `run_local()` on a fresh `bl_results`.
- **Self-checks** (`stop()` on failure):
  1. The baseline run 1 matches `numbers.rds`.
  2. `max_pairs` 10 and 21 give identical results.
  3. The CVA ratio equals n−g within a relative 1e-6 for every finite pair (unconstrained).
  4. The K3 configuration has the same group counts as the baseline.
  5. No CSV row has an empty `value` unless its `note` says why.

### 2. `documentation/methodology-changes/methodology-changes.qmd` → `.docx` (tracked)

- **Format:** `format: docx` (`toc: true`, `number-sections: true`); `execute: echo: false`. It reads only
  `method-change-effects.csv`. Numbers are formatted as in `review.qmd`
  (`formatC(x, format = "f", digits = 2)`, from full precision).
- **Structure:** as in "Content of the note" above. Decision boxes are two-row tables.
- **Prose:** load the `academic-writing-style` skill before drafting. The audience is the supervisor.
- **Render:** from that folder, run
  `"C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe" render methodology-changes.qmd
  --to docx`. Remove any `.quarto/` or `*_files/` it leaves, after checking against `git status`.

### 3. Small corrections found in planning

- **CLAUDE.md §6 (the CVA-identity bullet):**
  - mark "p.111"/"p.106" as pages of the 30 Sept 2025 draft;
  - add the final pages (rule p.86, §3.9.1; correspondence p.91);
  - say that the (n−g) factor comes from W/(n−g).
- **`NEWS.md`, the "No effect" entry:** a zero was "Supports" only for a class-0 target.
- **`.claude/reference/mahalanobis_implementation_plan.md`:** a dated erratum line (B ran `solve(Wmat)`
  with no divisor, on biplotEZ's `bp$Wmat`). The archive text stays verbatim.
- **The working plan's Task E tie-break bullet:** corrected (Finding 4).

### 4. Housekeeping after implementation

- **CLAUDE.md §2:** one line pointing to the note.
- The working plan's implementation log.
- `progress.md`.
- This plan archived as `.claude/reference/methodology-changes-taskE-plan.md`, with an IMPLEMENTED
  banner.
- The memory project file.

## Error cases

| Case | Action |
|---|---|
| The baseline does not reproduce `numbers.rds` | Stop. Re-run `run_all.R` at this commit and investigate before quoting a number. |
| A toggle gives `solution_found = FALSE` | Record it as a result; never re-pick the target. |
| A configuration's rule picks other targets | Record them. Run 1 still uses rows 23 and 25, so the change stays isolated. |
| An empty outcome group in a configuration | `stop()` with a message; mark it "not measurable". |
| T, A and B disagree | Report all three; T defines the thesis version. |
| T is silent | Cite the code and say "not stated in the thesis text". |
| A citation cannot be confirmed in the final PDF | Leave it out or mark it unverified. Never carry over a draft page. |
| An item is hard to classify as core or code | Core if it changes what the method computes (the search space, the feasibility rules, the constraint semantics, the attribution rules); code if it changes how the same thing is computed, a default parameter, or numerical behaviour. Say which in the history table. |
| The Quarto render fails (`rmarkdown` missing) | Report it; install nothing without asking. |

## Verification (test cases)

1. From `scripts/paper1-pima/` (PowerShell, `Rscript.exe`), run `run_all.R`, then
   `04_method-change-effects.R`. All five self-checks pass, and the meta rows name the commit.
2. **Sanity values:**
   - M4: 10 = 21;
   - M5: the same pair, ratio = n−g;
   - M8: TN 0.14 → 0.42;
   - M7: TN restores at 0.475 and 0.487;
   - M1's check: −27.4% and −35.2%.
3. `quarto render` produces `methodology-changes.docx`. Check it (`pandoc -t markdown` on the output)
   for:
   - the section order: history → core → code → future → interpretation;
   - a decision box on every core entry and one on the code section;
   - no hand-typed numbers (every number is an inline `r` call on the CSV);
   - the commit SHA shown;
   - every thesis citation a final-PDF page.
4. `devtools::test()`: 434 PASS / 0 FAIL (no `R/` change; a sanity run).
5. `git status`: only these change:
   - `documentation/methodology-changes/` (qmd, docx, csv);
   - `NEWS.md`;
   - the reference erratum;
   - `CLAUDE.md`;
   - `progress.md`.

**Then Gate 3:** report the `.docx` path and the drafts of D and E, and stop. No commit or push without
your instruction.

## Notes

- `documentation/` is tracked, so once pushed the `.docx` and the CSV are public on GitHub (the spec puts
  them there).
- **Run time:** about 3–5 minutes for the measurement script.
