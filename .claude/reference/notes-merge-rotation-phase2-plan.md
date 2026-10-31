> **IMPLEMENTED 2026-10-31** on branch `paper_nr1_changes`. Archived copy of the approved plan; the notes are `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` and `documentation/boundary-distance-and-robustness.md`.

# Plan: one note for rotation and plane selection; a separate note for Phase 2 distances

> **IMPLEMENTED 2026-10-31.** Uncommitted; see Implementation notes at the end.

Date created: 2026-10-30
Date last edited: 2026-10-31

## Context

The two notes contradict each other on how the best eigenvector pair is chosen.

- **`documentation/mahalanobis_technical_note.md`**
  - **What it argues:**
    - Z-distances are not comparable across pairs, because V is W-orthonormal (§2.2).
    - The pre-2026-05 Euclidean selector was "inconsistent across pairs and with the PhD's
      specification" (§2.3).
    - The two selectors "disagree when …" (§4.1).
  - **What else it holds:** the Phase 2 robustness material.
- **`documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`**
  - §3.3 shows that under CVA the two selectors are **identical**:
    d²_M = (n−g)·dist_z², with a ratio of 15,992 for every pair in the loan harness.
  - §4.1, §5 and §6 cover the rest: the problems, the measurements and the proposal.

What this session established (read-only checks):

1. **The CVA identity holds by construction, not only on the loan data.**
   - In biplotEZ 2.2, `CVA.biplot` sets `W` = the within-class SSP (no divisor) and
     `M = W^{-1/2} V_svd`, so `M'WM = I` exactly.
   - With the package's `W/(n−g)`, this gives `V'WV = I/(n−g)`.
   - It also holds for 2-group CVA. `CVAlowdim(low.dim = "sample.opt")`, the default the package
     uses, builds `cbind(M[,1], M[,-1] %*% fvek.opt)`, where `fvek.opt` is orthogonal.
   - Every rotated basis `V A'` inherits it.
   - The inverse is `V_rho^{-1} = (n−g) V_rho' W`. Hence, for in-plane points,
     `tVr W^{-1} tVr' = (n−g) I`, which gives the ratio.
   - Correction to both notes: biplotEZ 2.2 *does* expose `bp$Wmat`.
2. **The thesis (PDF, 30 Sept 2025) already states this, and fixes the criterion by projection:**
   - **p.106 (printed):** "the nearest counterfactual … corresponds to the shortest change vector
     … using either **Euclidean distance for PCA biplots and Mahalanobis distance for CVA
     biplots**", with `d_i = (z_cf − z_0)V_2'` or `(z_cf − z_0)M_2^-`.
   - **pp.102–103 and p.111:** "In the CVA biplot … the Euclidean distances on the biplot plane
     correspond to Mahalanobis distances in the original feature space."
   - **p.113:** CVA (TP/TN/FP/FN groups) is the projection for local and global interpretation
     "since it ensures that the estimated counterfactual reflects the shortest Mahalanobis
     distance".
   - **pp.139–140 (Pima):**
     - all p!/(2!(p−2)!) component pairs are tried (21 for p = 7), and 11 give constrained CFs;
     - the rotated matrix's first two components are used;
     - the CF must lie on a local boundary *within the convex hull of the training data*;
     - the CF is chosen by `min_i d_i' W^{-1} d_i`, with W the within-groups SSP;
     - the winner is pair (1,7).
   - **p.~148 (conclusions):** a "dual approach": Mahalanobis in feature space plus L2 on the
     plane.
   - **§6.2 (p.172), "Grand Tour instead of component selection":** search a smooth family of
     planes, rotate each to pass through the target, and score by Mahalanobis. This is the
     rotation note's §6.3 proposal.
3. **The research code (`origin/original_PhD_code`, read by an exploration agent) agrees.**
   - The pair search (`4.3 Different Eigenvalues.R`) runs under **CVA only**; under PCA there is
     no `Wmat`, so it errors.
   - It uses biplotEZ's `bp$Wmat`, the SSP with no divisor, so `pc.distz` is exactly ‖Δz‖².
   - The PhD's "Mahalanobis" selection was always the Z-distance selection.
4. **Why a different A per pair does not change this** (the user's question, 2026-10-31).
   - **What differs:** each pair has its own prediction matrix (M Aᵢ')⁻¹ = Aᵢ M⁻¹, so the same Δz
     maps to a different Δx.
   - **What does not:** M⁻¹ = M'W, so Δx = Δz J'AᵢM'W, and
     Δx W⁻¹ Δx' = Δz J'Aᵢ(M'WM)Aᵢ'J Δz' = ‖Δz‖². Aᵢ cancels as AᵢAᵢ' = I.
   - **Checked numerically** (iris, biplotEZ's CVA rebuilt by hand, `Rscript -`, no files
     written; the same Δz = (0.7, −0.4) for all six pairs):
     - the prediction matrices differ by up to 8.96;
     - ‖Δx‖² (Euclidean) ranges from 5.96 to 27.06;
     - Δx W⁻¹ Δx' = 95.55 for every pair, and the ratio to ‖Δz‖² is 147 = n−g exactly.
   - **So:**
     - Δz itself differs between pairs (each plane cuts the surface along a different curve);
     - the Mahalanobis and in-plane criteria rank the pairs identically;
     - the Mahalanobis criterion is correct, and the in-plane distance already measures it.
   - **When it holds:** CVA, the CVA's own W, and an in-plane displacement.
5. **Consequences:**
   - **Under CVA** the package matches the thesis. The Mahalanobis criterion and the in-plane
     Z-distance are the same quantity, up to n−g, so both `distance` options select the same pair.
   - **Under PCA** the thesis criterion is *Euclidean*, but the package default is
     `distance = "mahalanobis"` with `W_binary`, a within-class W by the 0/1 class on unstandardised
     X. So the PCA default departs from the thesis. `distance = "euclidean"` is the thesis's PCA
     rule.
   - The Mahalanobis note's §2.2–2.3 story is wrong for both projections. Z-distances *are*
     comparable across pairs, each in its projection's own metric.
6. **The package's max_pairs default caps the search below the thesis.** `max_pairs = 10L` is
   below choose(p,2) once p ≥ 6; the thesis's Pima example needed 21 pairs.
7. **Phase 2 diagnostic inconsistency (found, not to be fixed here).**
   - `plot.bl_boundary()` plots |v_j|/sqrt(W_jj) (an L1 sum).
   - Its "diagonal vs cross-correlation" diagnostic uses v_j²·(W⁻¹)_jj.
   - The two agree only for diagonal W. Record it as an open issue in the Phase 2 note.

The user asked for:

- **One note on rotation and best-component selection**, laid out as:
  - the original PhD method;
  - the current implementation;
  - current problems and the proposed solution;
  - future developments.
- **A separate note** for the robustness and distance-to-boundary calculation.

## Approach

### A. Combined note: `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` (rewritten in place)

Keep the path, so the folder's scripts and README stay with it. New title: "Local biplot rotation
and selection of the plane". New outline (old section → new):

0. **Header.** Purpose; the files in the folder, with the script table and new "backs" sections;
   a reading guide; a pointer to the Phase 2 note.
1. **Purpose and motivation** (old §1).
2. **The original PhD method.**
   - 2.1 The rotated component biplot: the three-point Procrustes construction, the solution, the
     dilation (defined, not applied), V^new, and the CVA variant (old §2.1–2.4, 2.6). Cite the
     final thesis §3.4 with page numbers. The old note cites "§3.5–3.5.2" of a draft, so
     re-check them against the PDF.
   - 2.2 Worked examples (old §2.5), with the old §5.8 caveat that the toy numbers do not
     reproduce.
   - 2.3 Searching all component pairs (thesis pp.114, 139–140):
     - the J selection matrix and all choose(p,2) pairs;
     - rotate each pair, then use the first two columns;
     - constraints: actionability, feasibility, the convex hull, and boundaries that real data
       support;
     - in-plane nearest point by L2 on the plane;
     - across planes, `min d'W⁻¹d` (W = SSP);
     - the Pima result: 21 pairs, 11 CFs, pair (1,7);
     - if nothing is found, "no constrained CF exists".
   - 2.4 The distance criterion in the thesis:
     - Euclidean for PCA and Mahalanobis for CVA (p.106);
     - under CVA, plane distance is Mahalanobis distance (pp.102–103, 111), so the "dual
       approach" is one metric under CVA;
     - why CVA is the thesis's local projection (p.113).
   - 2.5 The research code (`origin/original_PhD_code`; cite `scripts/<file>:<line>`):
     - **`Biplot_rotation()`** (`1.3 Optimal Rotation.R:32-81`):
       - Y = rbind(−x, 0, x);
       - `svd(t(YV) %*% YVr_padded)`, `A = v u'`;
       - ρ is computed as a diagnostic and never applied;
       - it returns the full `Vrho`/`tVrho`, and callers use columns 1–2.
     - **The pair search** (`4.3 Different Eigenvalues.R`):
       - pairs from `t(combn(p, 2))`, all choose(p,2) with no cap (49–50);
       - search data = the target plus *test* data (27);
       - a new grid and **hull rebuilt in each rotated plane** from the rotated test data, with
         contours clipped to it (`1.1:131-142`, `1.4:121-150`);
       - closed contours kept only if they enclose the target;
       - `set_filters` (range expressions only; no fixed or direction constraints) and
         `train_ranges` applied;
       - model re-score;
       - nearest vertex per contour by 2D Euclidean with `max.col(-d2)`, i.e. random ties
         (`1.4:45`);
       - nearest opposite-class contour.
     - **Selection:**
       - `pc.distz = v Wmat⁻¹ v'`, where v is the Z-difference back-projected by `tVrho[1:2,]`
         (162–166);
       - argmin after dropping zeros (178–186);
       - the winner is re-run.
     - **`Wmat`:** biplotEZ's `bp$Wmat`, the SSP with no divisor, on the TP/TN/FP/FN factor
       (`3.1:183-189`). **CVA only**: under PCA biplotEZ has no `Wmat`, so the search errors.
       Hence `pc.distz` = ‖Δz‖² exactly, the identity of §3.4.
     - **Known defects:**
       - a pair with no surviving contour returns B = 0, so its "distance" is the target's
         distance to the origin and it is not excluded (`1.4:100-113`);
       - the in-loop plot Z comes from uncentred data (`4.3:99`) and the final step from `bp$X`
         (`4.3:225`);
       - the default `3.2` flow does not rotate, and `4.4` Shapley uses the unrotated result.
3. **The current implementation.**
   - 3.1 `.bl_rotate()` code mapping (old §2.7). The CVA remark that "V is W-orthonormal" now
     points to §3.4, not to the Mahalanobis note.
   - 3.2 Pair enumeration and the `max_pairs = 10L` cap (old §3.1).
   - 3.3 The per-pair pipeline, with the deterministic nearest vertex (old §3.2).
   - 3.4 The selector and its metric (old §3.3, merged with the Mahalanobis note's correct parts
     of §2.2, §3.1 and §4.1):
     - `W` in one paragraph, pointing to the Phase 2 note for its definition and Cholesky
       inversion;
     - **the proof of the CVA identity, written to answer "each pair has a different A, so a
       different prediction matrix"** (Context item 4):
       - first state what does differ: the prediction matrix and Δx;
       - then show that Aᵢ cancels: biplotEZ's `M = W^{-1/2}V`, so `M⁻¹ = M'W`, and Δx W⁻¹ Δx'
         = ‖Δz‖²·(n−g);
       - the geometric picture: whitened coordinates, A as a rotation inside them, each plane an
         isometric slice;
       - the 2-group `sample.opt` case;
       - the conditions: CVA, the CVA's W, and an in-plane displacement;
       - the thesis's own statement (p.111);
       - a table of the iris demonstration (the same Δz for all six pairs: different ‖Δx‖²,
         identical Mahalanobis, ratio 147 = n−g) and the loan ratio of 15,992;
     - the PCA case (standardised Euclidean versus `W_binary`), with the measured disagreement
       (51 of 221 rows; the iris fixture);
     - the `distance` argument.
   - 3.5 Winner selection and the "fixed" restore (old §3.4).
   - 3.6 **How the implementation differs from the original** (a table):
     - CVA only in the research code, versus PCA or CVA in the package;
     - search data is test data in the code, versus training data for the package's grid range;
     - the hull is rebuilt in each rotated plane in the code (and required by the thesis), but
       never applied in the package's local search (CLAUDE.md §4; the global polygon is invalid
       after rotation);
     - closed contours must enclose the target in the code; the package has no such rule;
     - the code's B = 0 defect versus the package's `next` on no survivors;
     - `set_filters` and the fixed restore;
     - `train_ranges` is applied in the loop;
     - the exact nearest vertex replaces the random tie-break;
     - the `max_pairs` cap of 10 versus all pairs;
     - the PCA default (Mahalanobis `W_binary`) versus the thesis's Euclidean; the research
       code has no PCA path at all;
     - the W divisor n−g versus SSP (the argmin is unchanged);
     - W comes from biplotEZ-consistent data;
     - columns 1–2 of the rotated matrix in both.
4. **Current problems** (most consequential first). For each problem the main text gives the
   statement, why it happens, one headline number (with its script reference) and a pointer to
   the appendix table that holds the evidence (the user's choice: a condensed main text with an
   appendix):
   - 4.1 The plane is under-determined, and cross-platform reproducibility is untested (old §4.1;
     old §5.3 within-platform determinism; the CVA level-order instance).
   - 4.2 Pair labels are weak (old §5.2; old §5.1's dimensional ladder with the §6.1 correction).
   - 4.3 The local plot's orientation is arbitrary (the in-plane frame; `rotate_deg` is removed).
   - 4.4 The metric:
     - under CVA, the documentation and the `distance` argument suggest two different
       selectors, but they are one: the in-plane distance already measures the Mahalanobis
       distance, so the option cannot change the result (a documentation and API-clarity
       problem, not a methodological error);
     - under PCA, the in-plane and cross-plane metrics differ (old §4.2), and the default
       departs from the thesis (Euclidean).
   - 4.5 Constraints are applied after contouring (the `train_ranges` bias; old §5.7 and the §6.3
     fix 1 evidence).
   - 4.6 Off-plane counterfactuals: the restored "fixed" CF and the sparse CF. Their reported
     distances describe the search point.
   - 4.7 Only pair (1,2) can be forced, and `max_pairs` caps below choose(p,2) for p ≥ 6.
   - 4.8 A target at the centroid (from inspection only, not measured).
   - 4.9 Inherent limitations, documented and not to be "fixed":
     - global fidelity is traded for local accuracy (old §4.3, §5.6);
     - the result is not the nearest boundary in p-space (old §4.4);
     - CVA class means leave the plane (old §4.5);
     - a 2D biplot is a slice, not a view (old §5.4);
     - tilt is not a criterion (old §5.5);
     - adjacent mechanisms not to confuse with this one (old §4.6).
5. **Proposed solution** (old §6.1–6.4; a proposal, not implemented):
   - 5.1 Exactness requires only that the plane contains the target; dim 𝒜 = p−2.
   - 5.2 The minimal Givens rotation, its frame convention, and PCA versus CVA.
   - 5.3 Search the admissible family. State that this realises the thesis's own §6.2 Grand Tour
     proposal.
   - 5.4 One metric end to end, chosen by projection: under CVA it is already one metric; under
     PCA, decide between Euclidean (the thesis) and W_binary. Feasibility becomes a constraint.
   - 5.5 Reporting.
   - 5.6 Adoption path (opt-in `rotation =` / `plane_search =`; defaults unchanged until tested).
6. **Future developments:**
   - the open acceptance tests (old §6.5);
   - cross-platform/BLAS runs;
   - a `pairs =` argument and a `max_pairs` default tied to choose(p,2);
   - a centroid guard;
   - marking off-plane CFs or recomputing their distances;
   - the thesis's other future items that touch the local search;
   - what the proposal does not fix (old §6.6).
7. **Appendix A: measurements.** Every table now in old §4.1, §5.1–§5.6, §5.8 and §6.2–§6.3,
   moved verbatim with its script and section reference. Each table is labelled A.1, A.2, …, and
   the main text cites those labels.
8. **Appendix B: reproducing the package's objects** (old §5.7 practical notes).
9. **References**, with the thesis PDF by filename and page.

Rules for the rewrite:

- **Keep every current measurement and its script reference**, in the main text or in
  Appendix A. Re-word nothing numeric. The proposal's proofs (§5.1–5.2) stay in the main text,
  because they are the argument, not the evidence.
- **Present tense, current state.** History goes in single sentences, as in session 3's
  convention for this folder.
- **Equations:** keep the existing LaTeX style.

**Scope (the user's decision): documentation only.** No package code or tests change. The PCA
default's departure from the thesis is recorded as current problem 4.4, with the fix in §5.4; a
test pinning the CVA identity goes under future developments.

### B. Phase 2 note: rename `documentation/mahalanobis_technical_note.md` → `documentation/boundary-distance-and-robustness.md` (`git mv`, then rewrite; rename confirmed by the user)

Scope: how the package measures distance, and Phase 2's use of it. Outline:

1. **Scope**, and a pointer to the rotation note for Phase 3 pair selection.
2. **The metric W** (old §3.1–3.3, §5.1–5.2):
   - the definition, with divisor n−g;
   - the `metric_type` priority (`W_cva` → `W_binary` → `Sigma`);
   - it is computed on unstandardised X;
   - its relation to biplotEZ's `Wmat` (SSP; exposed in 2.2);
   - the Cholesky inverse and the ridge fallback;
   - V'WV = I/(n−g) under CVA, stated in one line with a link to rotation note §3.4.
3. **Distance to the nearest boundary in `bl_find_boundary()`:**
   - the nearest vertex in the global Z plane (2D Euclidean, exact, deterministic);
   - `dist_z`;
   - `dist_x` (standardised by X_sd);
   - under CVA, (n−g)·dist_z² is the Mahalanobis distance of the in-plane displacement.
4. **Per-variable distance and robustness** (`plot.bl_boundary()`, `bl_robustness()`):
   - why the denominator is within-class (old §2.1);
   - the mechanics and the `distance` table (old §4.2).
5. **Diagonal versus full Mahalanobis** (old §4.3): the Shapley trade-off and the diagnostic,
   plus **open issue 7** above (diag(W⁻¹) versus 1/W_jj; an L1 plot versus a squared
   diagnostic). Not fixed here.
6. **Edge cases** (old §5).
7. **Code references** (old §6).
8. **Future work** (old §7, minus the within-pair item, which moves to the rotation note).

The old cross-pair material (§1 table row, §2.2, §2.3, §4.1) moves to the rotation note,
corrected. The thesis p.174 Platt-scaling idea, which uses this distance, goes under future work.

### C. Ripple updates (mechanical, from an old → new section map)

- **Section references to the rotation note:**
  - the four `17_*.R` scripts' comments;
  - `documentation/biplot-rotation/README.md` (the "backs" column);
  - `documentation/nearest-vertex-tie-break.md`;
  - `documentation/fixed-constraint-restore.md` and `sparse-counterfactual.md`, if they cite it;
  - `development/reference/review_section9_to_15.md:322`.

  Find them with `grep -n "§[0-9]"` and "section [0-9]" in those files. Comments only; no code
  changes in the scripts.
- **Renamed-file references:**
  - `CLAUDE.md`:
    - §2: both entries, and the "Under review" line removed;
    - §5: the Cholesky pointer;
    - §6: the CVA line rewritten as settled, with the PCA-default finding;
    - §9: pointers.
  - `2 implementation_summary.txt`:
    - §2 file tree;
    - §4.2.1: "does not expose Wmat" and "W-orthonormal under exactly this W" corrected;
    - §4.8: the "Not comparable across pairs" row corrected, adding "Under CVA identical; under
      PCA the thesis rule".
  - `NEWS.md:56` (the path only, plus one Documentation bullet);
  - `development/reference/review_section9_to_15.md`;
  - `progress.md`.
- **Leave as they are:** the archived plans in `.claude/reference/` (historical records) and
  `docs/search.json` (generated; rebuilt with pkgdown later).
- **Memory:** update `project_boundarylogic.md` (its "next session" line) and its MEMORY.md
  hook. No memory file mirrors either note, so the sync rule does not apply.

## Verification

0. **Make the iris demonstration reproducible from the bundle.** Add it as a new asserted
   section of `documentation/biplot-rotation/17_3_minimal-rotation-verification.R` (a
   documentation script, not package code):
   - it uses the package's `.bl_rotate()` on an iris CVA `bl_result`;
   - it asserts that the Mahalanobis/‖Δz‖² ratio equals n−g for every pair;
   - it asserts that ‖Δx‖² differs between pairs.

   The note's table quotes this section's output.
1. **Re-measure the CVA claims** in a scratchpad driver after a harness run (CLAUDE.md §5),
   printing:
   - `max|t(V) %*% metric %*% V * (n-g) - I|`;
   - the `all_distances_mahalanobis / all_distances^2` ratio for the three local searches;
   - the same for a **2-group** CVA (`cva_classes = as.factor(train$class)`) on the loan data;
   - the PCA disagreement count (`17_1`'s 51/221 is already measured).

   Quote only the numbers measured there.
2. **Check every section reference resolves:**
   - grep the repo (excluding `.claude/reference/` and `docs/`) for `mahalanobis_technical_note`:
     expect 0 hits;
   - grep for `§` references to the rotation note and check each against the new headings.
3. **Package checks:** no package code changes, so `devtools::test()` should stay at 314 PASS /
   0 FAIL / 4 WARN. Run it once to confirm, and run `17_3` (its assertions are unaffected, but its
   printed section labels change).
   - **The `git mv` needs care:** the old note is tracked, and the rotation folder's changes are
     uncommitted. Use `git mv` only (no `rm`), and check `git status` against the session-start
     snapshot afterwards (CLAUDE.md §5).
4. **Rendering:** both notes render as Markdown, with LaTeX intact and no broken relative links.
   Check the links with a grep of `](` targets against `ls`.
5. Stop after checks; no commit (Parts B/C of session 3 are still uncommitted).

## Implementation notes (2026-10-31)

**Done as planned:** A, B and C, and verification steps 0–4.

**Additions and deviations:**

- **The `bl_find_local_cf()` roxygen section** "Cross-pair selection criterion" (`R/local_cf.R`)
  repeated the "not comparable across pairs" claim.
  - **Changed:** its text was corrected and `man/bl_find_local_cf.Rd` regenerated.
  - **No code changed.** The plan's "no package code" scope is read as no behaviour change.
- **The Phase 2 note also records the thesis's own Phase 2 method** (§4.1.2–4.1.3, pp.123–129):
  - it standardises by the total SD, with an L1 total T_j;
  - so `distance = "euclidean"` reproduces the thesis;
  - the within-class default is a 2026-05 package choice.
- **Measured:**

  | case | basis error | ratio, every pair |
  |---|---|---|
  | loan harness | 1.01e-13 | 15992, in the constrained, free and external searches |
  | loan harness, 2-group | 8.81e-14 | 15994 |

  In iris `17_3` §11, the prediction-matrix differences print as 16.561 / 20.649 / 13.099, at 3
  dp.
- **Found and not fixed:**
  - the harness's Step 16 external applicant is passed in raw units to a z-score model, so
    Z ≈ 115.76;
  - the installed `boundarylogic` is stale (2026-09-10), so drivers must `load_all()`.

  Both are recorded in `progress.md`.
- **The rename is staged** (`git mv`). Nothing is committed.

**Results:**

- **Tests:** 314 PASS / 0 FAIL / 4 WARN / 0 SKIP.
- **`17_3`:** all assertions pass, including §11.
