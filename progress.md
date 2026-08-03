# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Standing items / carry-overs

- **`main` was merged and pushed on 2026-08-03** (see session below) — the Mac-tester gate was
  consciously overridden to get the latest code (incl. the loan-default vignette) publicly
  installable. No further merge is pending unless new `method_developments` work accumulates.
- **`scripts/` is now gitignored** (`.gitignore`) — it is the user's personal scratch playground
  and is never published. `17_*`, `18_*`, `LendingClub_PD_Model_R.R`, `loans_full_schema.csv`,
  `01_iris_SVM_PCA_biplot.R`, and `pima diabetes.csv` all live there, untracked, on disk only.
  The former main dev/testing harness (`03_loan_status_Boundary_Logic.R`) and its three
  `review_section*.md` walkthroughs were moved to a new tracked `development/` folder (excluded
  from the built package via `.Rbuildignore`) — see session below.
- **`.bl_rotate()` plane under-determination** — found 2026-07-15, **documented and deliberately
  not actioned**. See `documentation/rotation_plane_underdetermination_note.md`, CLAUDE.md §9.
  Do not "fix" without an explicit instruction.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 136 PASS / 0 FAIL / 4 WARN** (the 4 WARN are pre-existing biplotEZ
  CVA 2-class notices).
- **Verification command:**
  ```r
  "/c/Program Files/R/R-4.6.0/bin/Rscript" -e "devtools::test()"
  ```
- **`.claude/plans/` is gitignored** (`.gitignore:16`) — plan files are local-only and never
  committed. `.claude/reference/` *is* tracked (generic docs + archived plans only, now that the
  script-specific walkthroughs moved to `development/reference/`).

---

## Session summary (2026-08-03, latest) — split scripts/ into playground + tracked development/, merge to main

Goal: publish the latest code to GitHub `main` for easy access (`main` was 19 commits behind
`method_developments`, 0 ahead), while keeping personal scratch scripts out of the published repo.

### What changed

- **`scripts/`** is now fully gitignored (`.gitignore`) — a personal playground, never published.
  It still holds all six files on disk (`01_iris_SVM_PCA_biplot.R`, `17_loan-gam-pca-3d-4var-surface.R`,
  `18_loan-glm-4var-pcabiplot.R`, `LendingClub_PD_Model_R.R`, `loans_full_schema.csv`,
  `pima diabetes.csv`); `git rm --cached` was used (never a raw `rm`), so nothing was deleted from
  disk. `outputs_dev_oos/` (an empty artefact dir from `LendingClub_PD_Model_R.R`) was gitignored
  pre-emptively too.
- **New tracked `development/` folder** holds what `scripts/` is not scratch for: the main
  dev/testing harness `development/03_loan_status_Boundary_Logic.R` (`git mv`'d from `scripts/03_...`,
  history preserved) and its three code walkthroughs, `git mv`'d from `.claude/reference/` to
  `development/reference/` (`review_section4_to_6.md`, `review_section7_to_8.md`,
  `review_section9_to_15.md`). Added `^development$` to `.Rbuildignore` so it never ships in the
  built package (R CMD check NOTE count unaffected).
- **Path references repaired**: `R/datasets.R` roxygen (now points at the vignette instead of the
  no-longer-published script path — `man/loan_data.Rd` regenerated via `devtools::document()`),
  `README.md` repo-structure block, `CLAUDE.md` §2/§5 (script-specific reference table, the
  reference-folder convention table, the `v2 block` pointer), `documentation/mahalanobis_technical_note.md`,
  and `documentation/scaling_axis_relabel_note.md`. The three review docs cross-reference each
  other by bare filename, so those internal links needed no edits. Auto-memory
  `reference_docs.md` synced to match, per the CLAUDE.md memory-sync rule.
- **Committed and merged**: this restructure, plus the previously-pending 2026-07-15 documentation
  work (`CLAUDE.md` §9 rotation-plane entry, `2 implementation_summary.txt` §4.2 caveat, the new
  `documentation/rotation_plane_underdetermination_note.md`), was committed on `method_developments`
  and fast-forward-merged into `main` (0 commits diverged, so no merge commit was needed), then both
  branches pushed to `origin`.

**Test status: 136 PASS / 0 FAIL / 4 WARN — baseline unchanged.** Confirmed after `devtools::document()`
regenerated `man/loan_data.Rd`; no `R/` logic changed, only a roxygen comment.

### Decisions

- **Mac-tester merge gate consciously overridden.** `progress.md` and CLAUDE.md §10 gate the
  `method_developments` -> `main` merge on macOS install confirmation, which has not landed. The
  user weighed "easy access to the latest code" (including the loan-default vignette, which existed
  only on `method_developments`) as the higher priority and explicitly asked to proceed.
- **`03_loan_status_Boundary_Logic.R` is not scratch.** Initially proposed folding all of `scripts/`
  into `.gitignore`; the user corrected this — script 03 is the main tool-logic development/testing
  harness and the subject of required-reading walkthroughs (CLAUDE.md §2), so it and its docs need a
  tracked home, not to be discarded. Resulted in the `development/` split rather than a blanket ignore.
- **Git history left alone.** The six scratch files remain findable in old commits (three were
  previously tracked). No `filter-repo`, no force-push — user's explicit choice.

### Next steps

1. None outstanding from this session — restructure, doc repairs, tests, commit, merge, and push
   are all complete.
2. Carried over: the `.bl_rotate()` minimal-rotation option (CLAUDE.md §9 + technical note §7) —
   still deferred pending an explicit instruction.
3. Optional, not done: rebuild the pkgdown site (`docs/`) — it still shows the old `scripts/` line
   in its rendered `README.md` mirror. Cosmetic; see CLAUDE.md §5 for the rebuild gotchas.

---

## Session summary (2026-07-15) — script 17 target-alignment rotation + rotation-plane finding

Three linked pieces of work on `scripts/17_loan-gam-pca-3d-4var-surface.R` (the user's 3D GAM
decision-surface demo), ending in a significant methodological finding about the package's
local biplot rotation. **No package code was changed** (`git diff -- R/` is empty).

### Modified files

| File | State | What changed |
|---|---|---|
| `scripts/17_loan-gam-pca-3d-4var-surface.R` | **untracked** | rotation added; `proj_pair` knob; retargeted to row 49 / pair (2,3); black curve made exact |
| `documentation/rotation_plane_underdetermination_note.md` | **untracked (new)** | the finding: derivation, measurements, options |
| `CLAUDE.md` | modified | §9 entry for the rotation finding (+ the still-uncommitted Shapley d.p. note carried from the previous session) |
| `2 implementation_summary.txt` | modified | §4.2 caveat under the `.bl_rotate()` description |
| `.claude/plans/rotate_3d_pca_to_target.md`, `rotate_3d_proj_pair_knob.md`, `compare_3d_vs_2d_biplot_target49.md` | gitignored | three approved plans, updated post-implementation |

**Test status: 136 PASS / 0 FAIL / 4 WARN — baseline unchanged, and NOT re-run.** No `R/` file
was touched, so the suite cannot have regressed. **Nothing committed.**

### 1. Target-alignment rotation (plan: `rotate_3d_pca_to_target.md`)

Script 17 now rotates the PCA basis onto the target so it lies exactly on the drawn plane,
simulating what `bl_find_local_cf()` does. A script-local `.bl_rotate_full()` mirrors the
package's private `.bl_rotate()` (`R/local_cf.R:45`) line-for-line but takes the **full 4x4** `V`
and returns the **full** rotated basis — the package's version needs a square `V` (it calls
`solve()`) and returns only 2 columns, while the 3D plot needs 3.

Key property (verified, 13 assertions): the target's rotated coords 3 **and** 4 are both ~1e-17,
so the 4D -> 3D reduction is **lossless for the target** (reconstruction error 8.88e-16) and
lossy for everyone else (0.0456 on an ordinary row). This is why the **full** basis must be
rotated — rotating a truncated 4x3 would merely *project* the target onto the plane.

### 2. `proj_pair` knob (plan: `rotate_3d_proj_pair_knob.md`)

Section 0 config block holds `i` and `proj_pair`. Decisions: **one rotation per run, no pair
search** (unlike `bl_find_local_cf()`); **provenance axis labels** (`c(2,4)` -> "Rotated PC2",
"Rotated PC4", "Rotated residual") since the pair's content always lands in slots 1-2 regardless
of which pair is named; **vertical axis fixed** to rotated column 3; internal `Z$PC1/PC2/PC3`
kept as *slot* names (renaming would ripple through the GAM formula, `acast`, and plotly refs).

### 3. Comparison against script 18 + the finding (plan: `compare_3d_vs_2d_biplot_target49.md`)

Retargeted to `i <- 49`, `proj_pair <- c(2L, 3L)` (the pair `bl_find_local_cf()` selects) and
replaced the black curve with an **exact level set** computed directly on the plane, instead of
contouring a `PC3 ~ s(PC1,PC2,k=60)` smooth fitted through a `|p-0.5|<0.02` cloud. `surf_mat_full`
became dead and was deleted; the green surface keeps its smooth (it is a genuine 3D object).

**The finding:** script 17's plane and the package's are **65.9 degrees apart** (principal angles
`0.00, 65.91`) — they share only the target direction. Root cause: `Y <- rbind(-x, 0, x)` is
**rank 1**, so `M = t(YV) %*% YVr_padded` is rank 1 (singular values `1.63, 6.7e-17, 0, 0`) with a
3-dim exactly degenerate null space. The Procrustes minimiser is pinned only along the target
direction; the plane's second direction is LAPACK's arbitrary tie-break. Bases differing by
**4.6e-15** give planes 66 deg apart. Only **46%** of the package's plane lies in the `span(V2,V3)`
it is named after. **Inherited from the original PhD method** — `origin/original_PhD_code:scripts/1.3
Optimal Rotation.R` does the identical padding; the package is a faithful port, not a porting bug.

Decision (user): **document, change nothing.** Written up in
`documentation/rotation_plane_underdetermination_note.md` with pointers from CLAUDE.md §9 and
`2 implementation_summary.txt` §4.2.

### Dead ends / corrections this session

1. **The planned degenerate `stop()` was justified by a false premise.** The plan specified
   `stop()` at tol `1e-8` because "`w ~ 0` makes the rotation undefined and the assertion fail
   opaquely". Measured: the rotation **never breaks** — probed down to a target with *exactly*
   zero score on both pair members, the basis stayed orthogonal and coords 3/4 stayed at 0. The
   SVD only needs `w`'s *direction*, which float noise supplies. Shipped a **`warning()` at
   `1e-12`** instead: what actually degrades is the target's *angular position within* the plane.
2. **"Rotation angle predicts view quality" — false.** Pair (2,3) rotates 65 deg yet retains
   *more* variance (0.9025) than pair (1,2) at 25 deg (0.8993); pair (2,4) rotates 86 deg and
   still retains 0.8989. Retention tracks how close the rotated 4th direction stays to the
   original PC4. Judge on the SS figure, never the angle. (Also killed the initial worry that
   pair (2,4) would be a bad choice — it is fine.)
3. **The prime suspect for the 17-vs-18 mismatch was wrong.** H2 (the `s(PC1,PC2)` smooth forcing
   single-valuedness) is real — 5.1% of cells hold 2+ separated PC3 roots, median curve error
   0.138 — but is **secondary** to H3, the plane mismatch (median 0.514). The exact-curve fix is
   a genuine improvement but does **not** make the two plots agree.
4. **Suspected a porting bug; there isn't one.** The under-determination is in the original PhD
   method. Checking `origin/original_PhD_code` before recommending a "fix" changed the entire
   nature of the recommendation.
5. **Ruled out early, cheaply:** the GAMs are byte-identical (`max|p17-p18| = 0.000e+00`; only
   `.pred_function()`'s 3-dp floor-round differs, 9.98e-04); biplotEZ `Lmat` **is** `svd(X)$v`
   with no sign normalisation and `e.vects` does **not** permute it; the four local-CF filters
   only shrink the *search set*, never the drawn curve.

### Architecture decisions / conventions

- **Script 17 mirrors the package's rotation; it does not call it.** `.bl_rotate()` cannot be
  reused (square `V` required, only 2 columns returned). The mirror cross-references
  `R/local_cf.R:45` as the authority.
- **Measure before guarding.** Two guards this session were designed from theory and corrected by
  measurement (see dead ends 1 and 2). The surviving guards report *both* diagnostics (angle and
  SS retained) and warn only on the one that is actually predictive.
- **`.bl_rotate()` behaviour is frozen pending an explicit instruction** — changing it moves every
  Phase 3 counterfactual, so it is not covered by the CLAUDE.md §4 "outputs proven equivalent"
  rule.

### Next steps

1. **Decide what to commit — nothing is committed.** Candidates: today's doc changes
   (`CLAUDE.md` §9, `2 implementation_summary.txt` §4.2, the new technical note), the
   *previous* session's still-uncommitted `CLAUDE.md` Shapley d.p. note, and this `progress.md`.
   **`scripts/17_*` and `scripts/18_*` need an explicit `git add`** or they will not be captured.
2. `progress.md` was 2 commits behind on entry (`0f47602`, `355e3e3`); both are now covered under
   History.
3. Carried over: Mac tester confirmation, then merge `method_developments` -> `main`.
4. Deferred by decision: the `.bl_rotate()` minimal-rotation option (CLAUDE.md §9 + technical
   note §7). If ever revisited, prototype **outside** the package and measure reproducibility,
   SS retained, and counterfactual Mahalanobis distances across a range of targets first.

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-06-17** — `355e3e3` round the Shapley plot counterfactual label to 2 d.p. (observed
  stays 3 d.p.); CLAUDE.md §9 note on the resulting d.p. inconsistency left uncommitted.
- **2026-06-16** — `0f47602` fix a stale import and non-ASCII chars in `biplot_grid.R`; doc updates.
- **2026-06-16** — raw-unit scaling support (`00b43a5`): `bl_set_scaling()` +
  `.bl_rescale_biplot_axes()` (raw-unit biplot axes), `.scale_to_raw()` (raw-unit feature values in
  Shapley/sparse/target prints; `kind="level"` vs `"delta"` is the correctness rule), and
  `bl_local$bl_counterfactual` mirroring `bl_target`. All display-only. Script 03 fixed to
  standardise on the training split only (was leaking test data). 77 -> 136 PASS.
- **2026-06-16** — `new_title = NA` extended to every biplot `plot()` method (`dc0c3c4`);
  established as a CLAUDE.md §5 convention. Dead end: the plan's predicted "invalid 'main'"
  error was wrong — `title()` coerces almost anything; only a closure/environment errors.
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
