# IMPLEMENTED 2026-05-29

> All five fixes implemented and verified (77 PASS, 0 FAIL, 4 WARN).

---

# Plan: Accuracy & Correctness Fixes

## Context

Exploratory audit of the `boundarylogic` package identified five accuracy/correctness issues in the core computation pipeline. None require changes to the mathematical method — all are defensive improvements: guarding against numerical instability, silent wrong-path execution, uncommunicated data exclusions, and non-independent random seeds.

Scope: moderate — fix the identified line(s) and improve directly related logic in the same function.

**Reviewed 2026-05-29 against current source:**
- Fix 1: unchanged; line reference updated (185 → 266).
- Fix 2: unchanged.
- Fix 3: moderate-scope addition (second warning) already implemented; only the first warning (line 263 fallback) remains. Updated plan accordingly.
- Fix 4: unchanged.
- Fix 5: unchanged.

---

## Fix 1 — Condition number check before `solve(V)` (`R/projection.R:266`)

**Problem:** `tV <- solve(V)` inverts the p×p loading matrix with no check for near-singularity. An ill-conditioned V produces a numerically unstable `tV`, which silently corrupts every back-projection (counterfactuals, boundary points, local CFs).

**Change:**
- After `V <- bp$Lmat` (line 265), before `tV <- solve(V)` (line 266), compute the condition number and warn if it is large.
- Added `condition_number = cond_v` to the `structure()` return block and documented in `@return`.

---

## Fix 2 — Clamp grid probabilities to [0, 1] (`R/biplot_grid.R:242`)

**Problem:** `col_vec[floor(grid_prob * 100) + 1L]` assumes `grid_prob in [0, 1]`. Some models can return values slightly outside this range, producing silent `NA` colours in the plot.

**Change:** Added `pmin(pmax(...))` clamp with warning after the chunked prediction loop.

---

## Fix 3 — Warn on boundary pruning fallback (`R/boundary.R:263`)

**Problem:** When no closed contours enclose any training observation, the code silently falls back to using ALL contours. The user sees results that look normal but may be unreliable.

**Change:** Replaced silent fallback with `warning(...)` before `boundary_used <- seq_len(nr_raw)`.

**Note:** The moderate-scope addition (second warning when all boundaries are emptied after consistency pruning) was already implemented at lines 311-322.

---

## Fix 4 — Report surrogate accuracy with hull coverage context (`R/surrogate.R`)

**Problem:** `accuracy_vs_model` silently excludes out-of-hull observations via `na.rm = TRUE`. A surrogate accuracy of 95% when 40% of data is out-of-hull is misleading.

**Changes:**
- Computed `hull_coverage <- length(idx_in) / n` after accuracy calculation.
- Added `n_total`, `n_in_hull`, `hull_coverage` to `structure()` return list.
- Updated `print.bl_surrogate()` to show `(in-hull: X / Y = Z%)` context on the accuracy line.

---

## Fix 5 — Shapley seed: default to `NULL` for independent calls (`R/shapley.R`)

**Problem:** `.shapley_perm_one()` defaulted to `seed = 1L` and called `set.seed()` unconditionally, giving identical permutation sequences across repeated calls.

**Changes:**
- Changed default to `seed = NULL` in both `.shapley_perm_one()` and `bl_shapley()`.
- Added `if (!is.null(seed)) set.seed(seed)` guard.
- Updated roxygen for both functions.
