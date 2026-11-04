# Technical Note: The Sparse Counterfactual — "No effect", Greedy Selection and Whole-Number Features

**Author:** Adriaan Rowan
**Date:** 2026-10-29; §3 revised 2026-11-04
**Package:** `boundarylogic`
**Branch:** `paper_nr1_changes`
**Location:** `documentation/sparse-counterfactual.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.9 -- code mechanics
- `../development/reference/review_section9_to_15.md` Steps 13–14 -- the harness walkthrough
- `fixed-constraint-restore.md` -- the related `"fixed"` restore in `bl_find_local_cf()`

Sections 1 and 2 record changes implemented on 2026-10-29. Section 3 records how `round_to` was
used, what was wrong with it, and the `"integer"` constraint that replaced it on 2026-11-04.

---

## 1. The "No effect" Shapley class (implemented)

**What changed.** `bl_shapley()` used to classify each feature as `"Supports"` or
`"Contradicts"`. For a class-0 target, a Shapley value `>= 0` meant "Supports"; for a class-1
target, `< 0` did. The rule comes from the PhD research script
(`origin/original_PhD_code:scripts/4.4 Local Interpretation Shapley.R`, line 102), and it counts
a value of exactly 0 as "Supports". A third class now comes first:
`abs(shapley_cause) < 1e-8` is `"No effect"`. The other rules are unchanged.

**Why an exact zero occurs.** A feature's Shapley value is a weighted average of its marginal
effects over all coalitions of the other features. It is exactly zero when moving the feature
from its observed to its counterfactual value never changes the prediction. Two package
properties make this common:
- **Tree models are piecewise constant.** A small move that crosses no split leaves the score
  unchanged.
- **Predictions are floor-rounded to 3 d.p.** (`.pred_function()`), which absorbs tiny changes.

The tolerance `1e-8` only absorbs floating-point noise. A non-zero value is built from
differences of at least 0.001, times weights of at least `1/p!`, so it is far larger.

**In the loan harness (target 13):**
- `credit_score` moves from 634 to 635 (an `"integer"` feature, rounded from 634.22 in
  `bl_find_local_cf()`, §3), which crosses no split: Shapley value 0.
- `loan_amnt` does not move at all (it was restored in `bl_find_local_cf()`): Shapley value 0.

Both were "Supports" under the old rule. Both are now "No effect".

**Effects.**
- `bl_find_sparse_cf()` never changes a "No effect" feature, so it stays at its observed value.
- `plot.bl_shapley()` colours it mid grey (`grey45`); "Contradicts" is light grey (`grey70`).

## 2. Greedy selection in `bl_find_sparse_cf()` (implemented)

**What changed.** `bl_find_sparse_cf()` used to change every "Supports" feature (except fixed
ones) to its counterfactual value, then score the result once. It now changes the "Supports"
features one at a time, in order of absolute Shapley value, largest first. It re-scores the model
after each and stops at the first point whose predicted class differs from the target's. The
features not reached stay at their observed values. "Contradicts", "No effect", "Unknown" and
`"fixed"` features never change. If no step crosses, every supporting feature has been changed,
which is the previous behaviour, and `solution_valid` is `FALSE`.

**The path is stored and printed.** `bl_sparse$sparse_path` has one row per step (`step`,
`variable`, `shapley`, `pred`, `crosses`), and `print()` shows it under the comparison table.
`used_in_sparse` is `TRUE` exactly for the features on the path.

**Why re-score instead of summing Shapley values.**
- A Shapley value is a feature's average effect over all coalitions, not its effect when it is
  added first. The running sum `pred_obs + sum(shapley)` therefore only approximates the
  prediction at a partial counterfactual.
- In the harness, `loan_int_rate` alone would reach 0.111 + 0.2915 = 0.403 by the sum. The model
  gives 0.360.
- Re-scoring uses the model's actual prediction. It costs at most one model call per candidate
  (`p` at most), which is negligible next to the local search.

**Not guaranteed minimal.** The order is fixed by the Shapley values. A smaller set of features
outside that order may also cross the boundary. Finding the smallest such set is a combinatorial
search the function does not attempt.

**Harness example (target 13, constraints `loan_amnt = "fixed"`, `loan_int_rate = "increase"`
and `credit_score = c("increase", "integer")`):**

| Step | Feature changed | Shapley | Prediction after | Crosses |
|---|---|---|---|---|
| 1 | `loan_int_rate` (12.69 → 14.146) | 0.2915 | 0.360 | no |
| 2 | `loan_percent_income` (0.05 → 0.0484) | 0.1445 | 0.547 | yes |

The sparse counterfactual changes two features instead of four. `credit_score` ("No effect")
and `loan_amnt` (fixed) stay at their observed values.

## 3. `round_to`, replaced by the `"integer"` constraint (2026-11-04)

### 3.1 How `round_to` was used

`bl_find_sparse_cf(bl_shapley_result, round_to = r)` rounded each feature changed on the greedy
path to the nearest multiple of `r`, as `round(value / r) * r`:

- it acted on the stored values, which are in model units;
- the point was re-scored after each step, so `solution_valid` reflected the rounded point;
- constraints were not re-checked;
- only the sparse counterfactual was rounded, not the counterfactual of `bl_find_local_cf()`.

**Harness example**, as it stood before the change: target 13, with `credit_score = "increase"`.
The model is standardised, so model units are z-scores.
`bl_find_sparse_cf(bl_shapley_values, round_to = 1)` gave:

| Feature | Observed | Counterfactual | Sparse, `round_to = 1` | z-score: counterfactual → rounded |
|---|---|---|---|---|
| `loan_int_rate` | 12.69 | 14.146 | 14.6504 | 0.833 → 1 |
| `loan_percent_income` | 0.05 | 0.0484 | 0.0666 | −1.202 → −1 |
| `credit_score` | 634 | 634.22 | 634 (not changed: "No effect") | — |
| `loan_amnt` | 12,000 | 12,000 | 12,000 (fixed) | — |

The path reached 0.384 after the first step and 0.483 after the second, so it never crossed:
`solution_valid = FALSE` and `pred_sparse = 0.483`. Without rounding the same two features cross
at 0.547. Before the "No effect" class existed, `credit_score` was rounded too: its z-score of
−0.221 rounded to 0, the training mean, and it jumped from 634 to 643.38.

**The problems.**

1. **Units.** For a standardised model the rounding was in z-scores, so `round_to = 1` meant
   "the nearest whole training SD". That is rarely what a user intends.
2. **One value for every feature.** Features have different scales: a whole point of credit
   score, a hundred units of loan amount, a hundredth of a percent of interest rate. A single
   `round_to` cannot suit all of them.
3. **Direction.** Nearest rounding can cancel a small move (back onto the observed value) or
   reverse it. In the example, `loan_percent_income` was meant to fall from 0.05 to 0.0484, but
   it rose to 0.0666.
4. **Constraints.** After rounding, `"increase"` / `"decrease"`, `c(min, max)` and
   `train_ranges` were not re-checked. A rounded value could violate them.
5. **Validity.** Rounding could undo the crossing, as in the example. The greedy path then used
   every candidate and still failed.

### 3.2 What using it does now

`round_to` is kept as an argument only so that old code passing `round_to = NULL` still runs.
Any other value stops with:

> 'round_to' is no longer supported. Mark whole-number features "integer" in set_filters(),
> e.g. set_filters(tgt, credit_score = c("increase", "integer")); bl_find_local_cf() then
> rounds the counterfactual in raw units and re-checks it, and the sparse counterfactual takes
> its values from it.

### 3.3 What replaces it

A whole-number feature is marked `"integer"` in `set_filters()`, alone or with a direction or a
range. The loan harness (Step 10) now uses:

```r
flt <- set_filters(tgt, loan_amnt = "fixed", loan_int_rate = "increase",
                   credit_score = c("increase", "integer"))
```

`bl_find_local_cf()` rounds the counterfactual after the search and the fixed restore
(`fixed-constraint-restore.md` §5). In the harness:

| Feature | Observed | Search's point (after the restore) | Counterfactual | Sparse |
|---|---|---|---|---|
| `loan_amnt` | 12,000 | 12,000 | 12,000 (fixed) | 12,000 |
| `loan_int_rate` | 12.69 | 14.146 | 14.146 | 14.146 |
| `loan_percent_income` | 0.05 | 0.0484 | 0.0484 | 0.0484 |
| `credit_score` | 634 | 634.2206 | 635 | 634 ("No effect") |

Both whole neighbours of 634.2206, 634 and 635, still cross at 0.547 and meet `"increase"`;
635 is nearer the target in the Mahalanobis metric (√d²_M 0.5914 against 0.5949). The
prediction stays 0.547, `credit_score` has Shapley value 0 ("No effect"), and the sparse
counterfactual is unchanged: `loan_int_rate` then `loan_percent_income`, crossing at 0.547.

**How it answers the five problems.**

1. **Units:** it rounds in raw units, converting with `bl_result$scaling` when one is set, so
   `"integer"` means a whole credit point, not a whole SD.
2. **One value for every feature:** it is set per feature, and only the features marked are
   rounded.
3. **Direction:** both the whole value below and the one above are tried, and a candidate that
   breaks a direction constraint is dropped. A move back onto the observed value is kept only if
   it still crosses.
4. **Constraints:** every `set_filters()` constraint and `train_ranges` is re-checked for each
   candidate.
5. **Validity:** each candidate is re-scored and kept only if it crosses; with none valid, no
   counterfactual is reported, with the reason. The Shapley values and the sparse counterfactual
   then start from a valid whole-number point, so the sparse counterfactual can no longer be
   invalidated by rounding.

**What it does not do.** It rounds to whole numbers only. A different step, such as 0.5 or 100,
would need a per-feature step size in `set_filters()` (`fixed-constraint-restore.md` §8).
