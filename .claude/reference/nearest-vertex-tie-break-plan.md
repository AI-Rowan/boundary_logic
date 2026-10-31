> **IMPLEMENTED 2026-10-30** on branch `paper_nr1_changes`. Archived copy of the approved plan; the technical note is `documentation/nearest-vertex-tie-break.md`.

# Nearest-vertex tie-break: make the boundary searches deterministic

Date created: 2026-10-30
Date last edited: 2026-10-30

## Implementation log (2026-10-30)

- **Approved by the user**, with my two recommendations:
  - the Gate 2 fixes were committed first, as `f4f3647` (not pushed);
  - in the rotation folder, correct §5.3, §5.9 and §6.2/§6.3 and the script comments, and leave
    the rest of §5 as dated measurements.
- **Done as planned:**
  - the one-line fix and the roxygen;
  - `tests/testthat/test-boundary.R` (6 tests). With the old line restored, 5 of them fail; the
    `which.min` test passes on random data, which has no near-ties;
  - the technical note `documentation/nearest-vertex-tie-break.md`;
  - the rotation note §5 dated line, §5.3, §5.9 (headings marked "explained and fixed"), §6.3,
    the §6.5 pointer, and the intro paragraph;
  - comments and printed text in `17_1`, `17_2` and `17_3`;
  - CLAUDE.md §2 (two entries), §6 (replaced bullet) and §9 (test coverage);
  - `2 implementation_summary.txt` §4.5, NEWS, and the walkthroughs (`review_section7_to_8.md`
    numbers and Stage 5; `review_section9_to_15.md:375`);
  - the Task C plan's Task D/E notes.
- **Not as planned:** the §6.2 claim ("§6.2 alone makes reproducibility a non-issue") sits in
  §6.3's first paragraph; it was corrected there.
- **Verification:**
  - `devtools::test()`: 314 PASS / 0 FAIL / 4 WARN.
  - Harness sourced twice in one session: identical logs. They match the patched-copy run from the
    checks.
  - Every changed walkthrough number comes from that run.
  - `R CMD check`: Status OK, 0 / 0 / 0.
  - `17_3` against the fixed repo: all assertions pass, and §8b gives full 1.147879 and pinned
    1.208670 on every call.

`bl_find_local_cf()`'s reported distance could differ between repeated calls in one session, and
depended on what had been called before it (CLAUDE.md §6; rotation note §5.3 and §5.9). The
mechanism was not established. It is now.

**Cause.** `.nearest_idx_block()` (`R/boundary.R:44`) picks the nearest contour vertex with
`max.col(-d2)`. `max.col()`'s default `ties.method = "random"` treats every entry within a
relative tolerance of 1e-5 *of the row's largest-magnitude entry* as a tie ("the determination of
a tie assumes that the entries are probabilities"; R's `?max.col`). It then draws one tied entry
from R's random-number stream.
- Here the entries are negated squared distances, and the largest magnitude belongs to the
  farthest vertex. So every vertex within 1e-5 × (farthest squared distance) of the nearest is
  "tied".
- The pick depends on the session's random-number state, which every earlier call that uses
  random numbers moves: `bl_find_boundary()`, `bl_surrogate()`, `bl_find_local_cf()`,
  `bl_prepare_data()` (`set.seed`) and `bl_shapley(seed =)`.
- The line is inherited verbatim from the PhD code
  (`origin/original_PhD_code:scripts/1.4 Biplot Boundary Search_optimized.R:45`, commented "index
  of min per row"). The PhD local workflow (`4.3 Different Eigenvalues.R`) uses the same search.

**The decision (user, 2026-10-30): option 1.** Take the exact nearest vertex,
`max.col(-d2, ties.method = "first")`. R documents this as equal to
`unname(apply(m, 1, which.max))`, the first exact maximum with no tolerance. Afterwards, write a
technical note.

## Checks done before deciding (2026-10-30)

All measured in the scratchpad: the patched helper in a copy of the package, and the repo untouched.

1. **Cause isolated (iris, standardised, GLM, PCA).**
   - Inputs, rotation, grid, model scores and contours are bit-identical across repeated calls.
   - `set.seed()` before each call gives 1 distinct result. No seed gives 2. One `runif(1)` changes
     the result; 200 `svd()` calls do not.
   - Deterministic tie-breaking gives 1 distinct result.
   - BLAS is R's 200 KB reference `Rblas.dll`, which is deterministic.
2. **Loan harness, current vs fixed, full console diff.**
   - **Unchanged:** Steps 11–15 for target 13 (constrained 0.547, path 0.360 → 0.547;
     unconstrained 0.5554 / 0.542), and every Phase 1 output.
   - **Phase 2 robustness:**
     - totals 9203.63 → 9203.29 (v1) and 6813.76 → 6813.80 (v2);
     - per-variable totals move in the 2nd decimal;
     - **the pruning ranking is unchanged**, so `feature_cols_v2` and `bl_results_v2` are
       identical.
   - **Surrogate accuracy** vs model: 0.8334 → 0.8337 (v1) and 0.8843 → 0.8845 (v2).
   - **External applicant (Step 16):**
     - the CF prediction moves from 0.488 to 0.368, and the sparse CF from 0.427 to 0.347;
     - √d²_M moves from 14638.44 to 14638.39;
     - the best pair is unchanged.
3. **Seed sweep (old code, 10 seeds, vs the exact nearest).**
   - **Global search:** 51 / 4,429 test rows (1.15%, v1) and 133 / 4,429 (3.0%, v2) get a
     different CF under at least one seed. The 2D distance moves by < 1e-4, but `B_pred` moves by
     up to 0.176 and `B_x` by up to 0.36 training SDs.
   - **Target 13:** identical under every seed.
   - **External applicant:** 4 distinct CFs under 10 seeds.
4. **After the fix:**
   - the global and local results are identical under other seeds and after other calls;
   - exact ties between *distinct* vertices are 0 of 62,024 queries, so `"first"` never falls back
     on vertex order in practice;
   - 0.79% of queries had more than one distinct vertex inside the old tolerance band.
5. **Tests on the patched copy:** 278 PASS / 0 FAIL / 4 WARN, the same as now. No test depends on
   the random pick.
6. **Rotation note scripts.**
   - **`17_3`** passes all its assertions under both versions. Its §8b (the §5.9 table, row 49)
     goes from "1.147879 or 1.153789" for the full search and "1.208670 or 1.202898" pinned, to
     1.147879 and 1.208670 on every call.
   - **§5.3 reproduced** (`17_1` set-up, five calls per row, sourced from the repo without
     changing it). With the current code, row 102 flips between (1,2) at 0.229143 and (1,3) at
     0.228868. With the fix it is (1,2) at 0.229143 every call. Rows 161, 140, 6 and 104 are
     stable under both.
   - **Consequence:** §5.3 attributes the flip to the rotation's plane under-determination (§4.1)
     and says §6.2 "removes the cause". Both are wrong: `.bl_rotate()` is bit-identical across
     calls. §5.3 and §5.9 are both this tie-break. §4.1 itself, sensitivity to input
     perturbation, stands, and so does §6's case for a determinate plane on that ground.

## Files to read first

- `R/boundary.R` (`.nearest_idx_block()`, `bl_find_boundary()`), `R/surrogate.R:160-200`,
  `R/local_cf.R:610-625`.
- `development/reference/review_section7_to_8.md` (Stage 5, the robustness and pruning blocks),
  `review_section9_to_15.md:375`.
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1, §5.3, §5.9.
- CLAUDE.md §6 (the "dist_mahalanobis can differ" bullet) and §9.

## Change

1. **`R/boundary.R:44`:** `min_loc <- max.col(-d2, ties.method = "first")`, with a two-line
   comment. The default `"random"` treats entries within 1e-5 of the row's largest entry as ties
   and draws one at random; `"first"` is the exact minimum.
   - The cross-block update `upd <- d2_min < best_d2` is already strict, so an exact tie across
     B-blocks keeps the earlier block. That matches "first" and needs no change.
2. **`.nearest_idx_block()` roxygen:** `@return` says it returns the first exact nearest row,
   deterministic, using no random numbers.
3. No other code changes. `bl_find_boundary()`, `bl_surrogate()` and `bl_find_local_cf()` all go
   through this helper.

## Error cases

| case | behaviour |
|---|---|
| `d2` row contains `NA` | `max.col()` returns `NA` for that row with either method; the callers' existing `idx[is.na(idx)] <- 1L` / `is.na(idx)` guards are unchanged |
| exact tie between distinct vertices | first in vertex order (deterministic); measured 0 occurrences in 62,024 queries |
| empty `B` | unchanged early return (`rep(NA_integer_, n)`) |

## Tests (new `tests/testthat/test-boundary.R`, plus `test-local_cf.R`)

- `.nearest_idx_block()`:
  - two distinct vertices whose squared distances differ by less than 1e-5 × the farthest (the
    old tie band), plus a far vertex: it returns the true nearest under 20 different seeds;
  - on random `Z`/`B` it equals `apply(d2, 1, which.min)`;
  - an exact tie returns the first index;
  - it leaves `.Random.seed` unchanged.
- `bl_find_local_cf()`: on the standardised-iris fixture that showed the flip (pair (1,3)),
  `all_distances_mahalanobis` and `B_x` are identical under two seeds, and after an intervening
  `bl_find_boundary()` / `bl_surrogate()` call.
- `bl_find_boundary()` and `bl_surrogate()`: `B_z` / `surrogate_pred` identical under two seeds.
  This adds the first tests of Phase 2 search output (CLAUDE.md §9 lists Phase 2 as untested).

## Documents

1. **Technical note (requested):** `documentation/nearest-vertex-tie-break.md`.
   - **The mechanism:** `max.col()` semantics, quoting `?max.col`; why a tolerance relative to the
     farthest vertex is wrong for distances; and that the band widens with distance (the external
     applicant).
   - **The diagnosis:** the bit-identical stages, the seed experiments, and that BLAS was ruled
     out.
   - **Origin:** the PhD code. Any thesis number produced by these searches carries the same
     random pick. Flag this for Task E (methodology changes).
   - **The fix,** and why `"first"` and not a seeded random draw.
   - **Measured impact** on the loan harness and the seed sweep (from the checks above).
   - **What it does not fix:** §4.1's plane under-determination, which is input-perturbation
     sensitivity. It is real, but separate.
   - **Consequences for earlier records:** rotation note §5.3 and §5.9, and CLAUDE.md §6.
2. **Rotation bundle** (`documentation/biplot-rotation/`):
   - **Note §5.3 and §5.9:** an "Explained 2026-10-30" paragraph each, correcting the attribution
     (the random nearest-vertex pick, not the plane) and giving the before/after numbers above.
   - **§6.2:** drop "removes the cause" of §5.3. §6's remaining case is §4.1.
   - **§5's other measurements** stay as measured, with one dated line saying they predate the
     fix.
   - **Script comments that repeat the old attribution:** `17_1:104`, `17_2:160-161`, the printed
     line under `17_3`'s "3) Determinism", and §8b's "practical rule" text. They become present
     tense, with a pointer to the note. No asserted values change.
3. **CLAUDE.md:**
   - §6: replace the "dist_mahalanobis can differ" bullet with one line on the exact,
     deterministic nearest vertex and a pointer to the note;
   - §2: list the new note;
   - §5: drop the "first call on a fresh `bl_result`" advice wherever it appears.
4. **`2 implementation_summary.txt`:** in the `bl_find_boundary()` passage (:289), say that the
   nearest vertex is the exact minimum and deterministic.
5. **Walkthroughs (from a fresh harness run):**
   - `review_section7_to_8.md`: Stage 5 (`max.col(-d2, ties.method = "first")`, with one sentence
     on why); the robustness console block, pruning list, pruning table and object-flow line take
     the new numbers. The ranking is the same, so no prose changes.
   - `review_section9_to_15.md`:375 gets a short note that the nearest is exact.
6. **NEWS:** a bug-fix entry saying results are now deterministic and the measured scale of the
   change.
7. **Plan archive:** `.claude/reference/nearest-vertex-tie-break-plan.md` (IMPLEMENTED banner).
8. **Task C plan:** Task D items 1–2 (the best-pair instability and call-history dependence)
   become "fixed, with a pointer".
9. **`progress.md`.**

## Verification

- `devtools::document()`, then `devtools::test()`: 278 + new PASS, 0 FAIL, the 4 baseline WARNs.
- **Harness:** fresh run, diffed against today's baseline log. The expected differences are
  exactly the lines listed in check 2; anything else is a stop-and-report.
- **Determinism:** the harness run twice in one session (`source()` twice) gives identical logs,
  which was impossible before.
- `R CMD check`: 0 / 0 / 0.
- `17_3` rerun: its assertions pass, and §8b shows one value per call type.
- The row-102 check gives (1,2) on every call.
