# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Next session: start here

**Commit the session-14 work first, if the user says so** (item 15), then **update the Paper 1 plan**
(item 1), then write its algorithm boxes (item 2). Start on branch `paper_nr1_changes` (CLAUDE.md
§10). The plan to update is
`../Explain XAI code/reference-materials/Paper_nr1/Exploring_boundaries/.claude/plans/paper-1-plan.md`
(last edited 2026-10-29). Plan the update first; what it must take in is under item 1. The hand-over
for the paper work is the Claude Doc
[boundarylogic v0.2.0: hand-over for Paper 1](https://claude.ai/artifact/Rb8YYhaNf43Xdj2dyZvyHz).
**Open decision for the user** (item 16): the Phase 2 full Shapley decomposition.

---

## Outstanding items (reviewed 2026-11-16)

Items 3–14 are closed (below, one line each). Nothing open is a blocker.

### Paper 1

1. **Update `paper-1-plan.md`.** It must take in:
   - **parameters and wording:** `max_pairs` is 21, not 10 (paper plan lines 11, 38, 76, 242);
     remove "no hull after rotation" and "call-history dependence" (both resolved);
   - **the method:** the convex hull of the training data on each rotated plane, with the empty-region
     rule (B10; `documentation/local-convex-hull.md`); the switches in `set_filters()` (the
     unconstrained counterfactual is one call); the variable importance measures the in-plane gap,
     while the global distances are reported only (`algorithms.md` F3); a target inside every
     training range can lie outside the hull on a plane (`local-convex-hull.md` §3);
   - **the Pima example (changed 2026-11-15):** three targets, TP row 23, TN row 25 and **FP row 26**
     (the misclassified patient; fixed restores that succeed, a "No effect" Shapley value, a sparse
     counterfactual of two features), each its group's median row by a rule that does not rest on
     rounding; **run 8** (*Glucose* in 70 to 140, the range constraint; it binds for the TP target);
     the made-up **patients C** (empty-region rule) **and D** (training ranges), beside A and B. Decide
     which of these go in the paper: Table 1's layout (three targets, runs 1, 3–8), the FP figures
     (`fig-fp-a` to `-d`), the rule-patient figures (`fig-rule-*`), `fig-ext`;
   - **figures and Table 1:** figure 6 is the dot chart only (no 6(b)); figure 7(a) is the TP target's
     run 4 searched without the hull; figures 5(b), 7(a), 7(b), 7(c), 9(a), 9(b) changed;
   - **the TP narrative:** plausibility against actionability (the issue 11 sentence chosen in the
     reply form, now the appendix of `documentation/future-improvements-review.md`): runs 2 and 4,
     not run 3;
   - **§5 (local interpretation):** the hull as part of the method, and the A2, B3 and B7 sentences;
   - **§8 (properties and limitations):** the A1 result (upper-bound sentence, `05_full-space-gap.R`;
     the FP target's ratio is 1.861) and the B2 platform sentence;
   - **§9 (further work):** the four themed paragraphs built from the agreed sentences (G1 in the
     form), drawing on `documentation/future-improvements-review.md` Part 1;
   - **the Pima paper sentences** for issues 5, 9, 11, 12 and 13 (`scripts/paper1-pima/issues.md`;
     issues 15–17 record the new targets and points);
   - **the final thesis's page numbers** (CLAUDE.md §5 anchors);
   - **tasks T3** (paper repo infrastructure, TinyTeX) **and T4** (port `scripts/paper1-pima/` to the
     paper repo, switch `00_setup.R`'s `devtools::load_all("../..")` to `library(boundarylogic)`, check
     `numbers.rds` matches); figures 1–2 (schematics, made in the paper repo).
   Source: the archived Pima and task plans (`.claude/reference/paper1-pima-example-plan.md`,
   `paper1-task-spec.md`), `.claude/reference/pima-targets-docs-consistency-plan.md`, the hand-over
   Claude Doc.
2. **The paper's algorithm boxes** (stage 2 of `documentation/algorithms.md`), straight after item 1:
   four LaTeX boxes for CVA in the paper's notation, with a traceability table. Plan:
   `.claude/plans/before-updating-the-paper-wiggly-wave.md`. When it starts, refresh its stale lines:
   it cites "§7 findings and decisions" (now "Behaviours of the code") and `17d3c18` (now `v0.2.0`).

### This repo

15. **Commit the session-13 rest and session-14 work** when the user says so: on `paper_nr1_changes`
    (the session-13 commit `c866b27` is local, not pushed), push, then fast-forward `main` and
    `method_developments`. Documentation, R/ text (roxygen and comments only) and regenerated help
    pages, so no reinstall and no new tag. Then add the commit SHA to the hand-over Claude Doc
    ("Install and versions") together with the third target and the target rule. The staged index
    already holds three changes (`git mv` of the reply form, `git rm` of `future-developments.md` and
    of the old Pima HTML export): commit them deliberately (CLAUDE.md §5 on `git mv`).
16. **Decide on the Phase 2 full Shapley decomposition** (`future-improvements-review.md` §1.6). With
    d²_M = v'W⁻¹v as the value function the exact Shapley values have a closed form,
    φ_j = v_j (W⁻¹v)_j, which sums to the full distance and needs no model calls; the documents'
    reason for not using it ("O(n·2^p) model scorings") does not hold. It appears in
    `boundary-distance-and-robustness.md` §4.1, `2 implementation_summary.txt` §4.6.1, the
    methodology note §5 and B14's agreed paper sentence. A per-feature φ_j can be negative (a cross
    term), a possible reason not to plot it as a distance. The decision is unchanged until the user
    rules.
17. **Supply reasons for 17 future items marked "Reason: not recorded"** in
    `future-improvements-review.md` (the reply form gave none), if wanted.

### Closed

- 3. Hand-over: the Claude Doc above (2026-11-15).
- 4. Gate 3: the supervisor approved the methodology note, stored as
  `documentation/methodology-changes/methodology-changes.md` (2026-11-15).
- 5–8. Tag `v0.2.0`, `method_developments`, public install, old Claude Doc (2026-11-14).
- 9–10. Rotation-script loading fix committed; Pima review page re-rendered (2026-11-14).
- 11. Root documents corrected against the code (2026-11-16; the 2026-11-15 refresh had covered only
  trees and counts).
- 12. Raw-unit label test for PCA without standardisation added (17 expectations; 2026-11-15).
- 13. `local-convex-hull.md`'s missing §7: CLAUDE.md §2 and the notes now cite §3 and §6 (2026-11-16).
- 14. Every technical note, the implementation summary and the walkthroughs follow the current-code
  rule (2026-11-16).

---

## Standing items / carry-overs

- **Paper 1 work is on branch `paper_nr1_changes`** (branched from `method_developments` at
  `b7aee65` on 2026-10-29, pushed with upstream). Start Paper 1 sessions on that branch
  (CLAUDE.md §10).
- **`paper_nr1_changes` is at `c866b27` locally** (the session-13 close-out, not pushed);
  `origin/paper_nr1_changes`, `main` and `method_developments` are at `e41f22f`. **Tag `v0.2.0` =
  `0689341`.** The session-14 work is uncommitted (item 15). Commit and push only when the user says
  so.
- **The installed `boundarylogic` is 0.2.0** (reinstalled 2026-11-07). Reinstall after each merge to
  `main` that changes package code (CLAUDE.md §5); the session-14 R/ changes are text only.
- **`scripts/` is gitignored** (the user's playground). The Paper 1 Pima example lives in
  `scripts/paper1-pima/` (local only): review page `review.html` (with a coverage table of every
  behaviour of the method and where it is shown), `issues.md` (issues 1–17; the two Gate 2 review
  files were deleted on 2026-11-15), `05_full-space-gap.R`, `06_findings-effects.R`,
  `07_rule-effects.R`, the figure scripts `fig-ext.R`, `fig-fp.R`, `fig-rule.R`. How to run it:
  CLAUDE.md §5. Three targets, eight runs, four made-up patients; `tests.R` has 18 checks.
- **The methodology note** is `documentation/methodology-changes/methodology-changes.md`, plain
  Markdown measured at `e41f22f` (R code = `v0.2.0`), approved by the supervisor; edit it by hand
  (CLAUDE.md §5). `04` and `06` measure runs 1–7 and the TP and TN targets only; re-run on
  2026-11-15, every value identical.
- **Future developments:** `documentation/future-improvements-review.md` is the one place (Part 1,
  35 items with reasons; the 2026-11-05 reply form as its appendix). CLAUDE.md and the notes only
  point to it. None is for Paper 1.
- **The current-code rule** (CLAUDE.md §2): technical notes describe the current code; the history
  lives in the methodology note, `NEWS.md` and the reply-form appendix. Every note now follows it.
- **Test status: 578 PASS / 0 FAIL / 4 WARN / 0 SKIP** (2026-11-16; the 4 WARN are biplotEZ 2-class
  CVA notices). At tag `v0.2.0`: 561 PASS. **`R CMD check`: 0 errors / 0 warnings / 0 notes**
  (2026-11-16). Pima `tests.R`: 18 of 18.
- **Verification:** run R as `.R` files (`Rscript.exe <driver>.R`), with drivers that
  `devtools::load_all()` the repo. The rotation scripts load the working tree themselves.
- **`.claude/plans/` is gitignored** (local only); `.claude/reference/` is tracked.

---

## Session summary (2026-11-15 to 2026-11-16, session 14, latest)

Plan: `.claude/reference/pima-targets-docs-consistency-plan.md` (archived; working copy
`.claude/plans/before-the-next-steps-temporal-storm.md`).

- **Committed** the session-13 work as `c866b27` (not pushed); a `Bash(git commit *)` allow rule was
  added to `.claude/settings.local.json` at the user's request.
- **Pima example** (local): a read-only scan of all 79 test rows under every run showed the TP target
  has a counterfactual in 2 of 6 Table 1 runs. Added the FP target (row 26) and a target rule that
  does not rest on rounding (rank `floor(n/2) + 1`; TP and TN unchanged); run 8 (the range
  constraint); patients C and D (empty-region rule, training ranges; a namespace patch turns the
  empty-region rule off for C's comparison); `07_rule-effects.R` (the ranges change 6 of 79 test
  rows, the empty-region rule 1). Every earlier number is identical; `tests.R` 18 of 18. The two Gate 2
  review files were deleted (their content is in `issues.md`).
- **Documentation audit** (about 100 findings) and corrections, by five then four parallel agents
  plus the coordinator: the Foundation document and the implementation summary against the code;
  NEWS.md 0.2.0 rewritten as one consistent section (dated 2026-11-07) with the missing features;
  README (version, tagged install); the vignette (precompile passed); the walkthroughs; seven notes
  (current harness numbers); `algorithms.md` (output lists, signatures); `local-convex-hull.md`
  (patient C as a property of the method); CLAUDE.md (corrections, missing files and helpers, merged
  rules, §9 only a pointer); R/ roxygen and comments (text only; 30 non-ASCII characters replaced;
  two message wordings), 9 help pages regenerated; the memory copy of
  `reference_data_prep_functions.md` synced; the old Pima HTML export removed.
- **One future-developments file:** `future-improvements-review.md` moved to `documentation/` and
  merged with `future-developments.md` (removed); every note's future section is a pointer.
- **History sweep** of every remaining note and the implementation summary (item 14), with the
  conditional numbers re-measured by one-rule patches (tie-break, empty-region rule, diagonal check,
  CVA level order, PCA selectors); history not recorded elsewhere went to NEWS.md.
- **pkgdown site** rebuilt; CLAUDE.md and progress.md removed from its search index and sitemap.
- **Found for the user:** the closed-form Shapley decomposition (item 16); the importance's CVA
  within-class SD is the package's choice, not the thesis rule (the documents now say so).
- **Tests:** 578 PASS / 0 FAIL / 4 WARN; `R CMD check` 0/0/0. Blockers: none.
- **Next:** commit on the user's word (item 15), then the Paper 1 plan update (item 1).

---

## Session summary (2026-11-14, session 13)

- Release and branch housekeeping (items 5–9): the session-12 changes committed and pushed; tag
  `v0.2.0` on `0689341` pushed; `main` and `method_developments` fast-forwarded and pushed; the
  public `install_github(...@v0.2.0, build_vignettes = TRUE)` checked; the old Pima Claude Doc
  deleted. No package code changed: 561 PASS stands.
- **Finding F4 in the Pima example** (plan `.claude/plans/lazy-growing-token.md`, two revisions):
  - `scripts/paper1-pima/` (local) now searches two made-up external patients, under the four filter
    switches with run 1's filters:
    - **Patient A** (Pregnancies 17): the hull rule. With the hull, pair (2,4), √d²_M 3.20; without
      it, pair (1,4), 2.51, outside the training data. The same before the fix.
    - **Patient B** (Insulin 900): the fix, as the loan applicant shows it. Before the fix, pair
      (2,3), 0.77; now pair (2,4), 2.45, Insulin 641.45 instead of 819.01.
  - Changed files: `00_setup.R` (`external_patients`, `external_settings`, `run_local(train_ranges
    =)`), `03_local.R`, new `fig-ext.R` (in `run_all.R`), `tests.R` (check 10), `06_findings-effects.R`
    (`M13ext` rows in `findings-effects.csv`), `review.qmd` (new section, re-rendered; it now names
    `e41f22f`, which closes item 10).
  - Notes and instructions: `algorithms.md` F4, `local-convex-hull.md` §6, CLAUDE.md §5.
  - Results: Table 1 and every run are unchanged, and `tests.R` passes 14 of 14 (8a–8c included).
  - **Found on the way:**
    - A patient with 16 pregnancies (Glucose 180, Insulin 300, Age 50) changes with the fix only
      through the local **empty-region rule** (F2's local part), not F4. A scan of 92 patients found
      15 clean F4 cases.
    - `local-convex-hull.md` refers to a "§7" (the hull of the training data alone and the empty-region
      rule) that the file does not contain: it ends at §6. Restore or repoint it (not done; new item
      13).
  - **For the Paper 1 plan update:** decide whether `fig-ext` goes in the paper.
- **`algorithms.md` describes the current code only** (user instruction):
  - §7 is now "Behaviours of the code": no decisions, dates, commits or plans. Each F-item states what
    the code does and what would happen without it.
  - The other notes' "since 2026-…" markers are gone. The same rewording was applied to
    `local-convex-hull.md` §6 and to the review page.
  - **Each rule was re-measured on its own** (namespace patches in scratch drivers): adding the
    target to the hull, or switching off the empty-region rule.
    - The loan and Pima numbers in F2 and F4 are unchanged.
    - The old "4 of 79 Pima test rows" split: row 16 depends on the empty-region rule (F2, 1.9%
      further), and rows 46, 57 and 63 on the training-only hull (F4, 7.1% further at the median).
    - No run of the two Paper 1 targets depends on either rule.
- **2026-11-15: the outstanding items closed** (plan `.claude/plans/lazy-growing-token.md`):
  - **Methodology note:** re-measured at `e41f22f` (`04`: six self-checks passed, values unchanged;
    `05`: median ratio 1.187), rendered one last time and stored as `methodology-changes.md`; the
    `.qmd` and `.docx` were removed (`git rm`). CLAUDE.md §2/§5, the notes and the scripts point to
    it.
  - **`documentation/future-developments.md`** (new): every future development (method, code, tests,
    release), merged from the note's Table 5, CLAUDE.md §9, the review form (A2, A4–A7, B1, B2, B4–B8,
    B11–B18, B21–B23) and the notes' future sections. CLAUDE.md §9 slimmed to the rule, the
    `.bl_rotate()` warning and a pointer. Two pieces of guidance from §9's resolved entries moved:
    GAM wrapping (§7) and `bl_pick_point()` with `rotate_deg` (§6).
  - **Current-code rule** in CLAUDE.md §2; `local-convex-hull.md` rewritten to it (and its "§7"
    references repointed).
  - **Root documents:** trees, counts, dates and branches refreshed.
  - **Test:** exact raw-unit labels for a PCA without standardisation (`test-set_scaling.R`, 17
    expectations); 578 PASS / 0 FAIL.
  - **`issues.csv` → `issues.md`** (13 issues, every field checked by script); the review page
    includes it.
  - **Hand-over:** the Claude Doc
    [boundarylogic v0.2.0: hand-over for Paper 1](https://claude.ai/artifact/Rb8YYhaNf43Xdj2dyZvyHz).
  - **Plan archived:** `.claude/reference/outstanding-items-closeout-plan.md` (the session's two
    earlier plans are recorded only in this summary).
- **Session closed 2026-11-15.** Modified (uncommitted): CLAUDE.md, the two root documents,
  `algorithms.md`, `local-convex-hull.md`, `boundary-distance-and-robustness.md`,
  `fixed-constraint-restore.md`, `future-improvements-review.md`, both effect CSVs, `progress.md`,
  `test-set_scaling.R`; new `future-developments.md`, `methodology-changes.md`, the plan archive;
  removed `methodology-changes.qmd`/`.docx`. Local only: `scripts/paper1-pima/` (`00_setup.R`,
  `03_local.R`, `04`–`06`, `fig-ext.R`, `run_all.R`, `tests.R`, `review.qmd`/`.html`, `issues.md`;
  `issues.csv` removed). Tests: 578 PASS / 0 FAIL / 4 WARN; Pima `tests.R` 14 of 14. Blockers: none.
- **Next:** commit on the user's word (item 15), then the Paper 1 plan update (start here, above).

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
