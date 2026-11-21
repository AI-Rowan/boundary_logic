# Technical Note: The Sparse Counterfactual — "No effect", Greedy Selection and Whole-Number Features

Date created: 2026-10-29
Date last edited: 2026-11-16

**Author:** Adriaan Rowan
**Package:** `boundarylogic`
**Location:** `documentation/sparse-counterfactual.md`
**Companion documents (paths relative to this file):**
- `../2 implementation_summary.txt` Section 4.9 -- code mechanics
- `../development/reference/review_section9_to_15.md` Steps 13–14 -- the harness walkthrough
- `fixed-constraint-restore.md` -- the related `"fixed"` restore in `bl_find_local_cf()`
- `methodology-changes/methodology-changes.md` §3.8–§3.10 -- the history of these rules

Section 1 covers the "No effect" class of `bl_shapley()`, Section 2 the greedy selection of
`bl_find_sparse_cf()`, and Section 3 how whole-number features are handled.

---

## 1. The "No effect" Shapley class

**The rule.** `bl_shapley()` classifies each feature as `"No effect"`, `"Supports"` or
`"Contradicts"`. A feature with `abs(shapley_cause) < 1e-8` is `"No effect"`. Otherwise, for a
class-0 target a Shapley value `>= 0` means "Supports" and `< 0` "Contradicts"; for a class-1
target, `< 0` means "Supports" and `>= 0` "Contradicts". The thesis's code (Appendix A, p.216,
the `Contribute` labelling) has no "No effect" class: it counts a value of exactly 0 as
"Supports" for a class-0 target and as "Contradicts" for a class-1 target.

**Why an exact zero occurs.** A feature's Shapley value is a weighted average of its marginal
effects over all coalitions of the other features. It is exactly zero when moving the feature
from its observed to its counterfactual value never changes the prediction. Two package
properties make this common:
- **Tree models are piecewise constant.** A small move that crosses no split leaves the score
  unchanged.
- **Predictions are floor-rounded to 3 d.p.** (`.pred_function()`), which absorbs tiny changes.

The tolerance `1e-8` only absorbs floating-point noise. A non-zero value is built from
differences of at least 0.001, times weights of at least `1/p!`, so it is far larger.

**In the loan harness (target 13, Step 13):**
- `credit_score` moves from 634 to 635 (an `"integer"` feature, rounded from 634.0454 in
  `bl_find_local_cf()`, §3), which crosses no split: Shapley value 0.
- `loan_amnt` does not move at all (it is "fixed", and `bl_find_local_cf()` restores it):
  Shapley value 0.

Both are "No effect". Target 13 is a class-0 target, so without the "No effect" class,
Appendix A's rule would label both "Supports".

**Effects.**
- `bl_find_sparse_cf()` never changes a "No effect" feature, so it stays at its observed value.
- `plot.bl_shapley()` colours it mid grey (`grey45`); "Contradicts" is light grey (`grey70`).

## 2. Greedy selection in `bl_find_sparse_cf()`

**The rule.** `bl_find_sparse_cf()` changes the "Supports" features one at a time, in order of
absolute Shapley value, largest first. It re-scores the model after each and stops at the first
point whose predicted class differs from the target's. The features not reached stay at their
observed values. "Contradicts", "No effect", "Unknown" and `"fixed"` features never change. If no
step crosses, every "Supports" feature has been changed, as it would be without the stopping
rule, and `solution_valid` is `FALSE`.

**The path is stored and printed.** `bl_sparse$sparse_path` has one row per step (`step`,
`variable`, `shapley`, `pred`, `crosses`), and `print()` shows it under the comparison table.
`used_in_sparse` is `TRUE` exactly for the features on the path.

**Why re-score instead of summing Shapley values.**
- A Shapley value is a feature's average effect over all coalitions, not its effect when it is
  added first. The running sum `pred_obs + sum(shapley)` therefore only approximates the
  prediction at a partial counterfactual.
- In the harness's search without actionability constraints (Step 15), `loan_int_rate` alone
  would reach 0.111 + 0.2905 = 0.4015 by the sum. The model gives 0.360. With the constraints
  (Step 14) the sum gives 0.111 + 0.4405 = 0.5515, and the model 0.512.
- Re-scoring uses the model's actual prediction. It costs at most one model call per candidate
  (`p` at most), which is negligible next to the local search.

**Not guaranteed minimal.** The order is fixed by the Shapley values. A smaller set of features
outside that order may also cross the boundary. Finding the smallest such set is a combinatorial
search the function does not attempt.

**Harness example (target 13, Step 14; constraints `loan_amnt = "fixed"`,
`loan_int_rate = "increase"` and `credit_score = c("increase", "integer")`):**

| Step | Feature changed | Shapley | Prediction after | Crosses |
|---|---|---|---|---|
| 1 | `loan_int_rate` (12.69 → 14.44) | 0.4405 | 0.512 | yes |

The sparse counterfactual changes one feature, where the full counterfactual (prediction 0.693)
changes three. `loan_percent_income` ("Supports", 0.1415) is not needed; `credit_score`
("No effect") and `loan_amnt` (fixed, "No effect") stay at their observed values.

**Without actionability constraints (Step 15)** the full counterfactual (prediction 0.542) moves
all four features, but `loan_amnt` and `credit_score` have Shapley value 0 ("No effect"), and the
path takes two steps:

| Step | Feature changed | Shapley | Prediction after | Crosses |
|---|---|---|---|---|
| 1 | `loan_int_rate` (12.69 → 14.11) | 0.2905 | 0.360 | no |
| 2 | `loan_percent_income` (0.0500 → 0.0399) | 0.1405 | 0.542 | yes |

## 3. Whole-number features: the `"integer"` constraint

`bl_find_sparse_cf()` does not round: its `round_to` argument is deprecated, and any value other
than `NULL` stops with an error that points to the `"integer"` constraint.

A whole-number feature is marked `"integer"` in `set_filters()`, alone or with a direction or a
range. The loan harness (Step 10) uses:

```r
flt <- set_filters(tgt, loan_amnt = "fixed", loan_int_rate = "increase",
                   credit_score = c("increase", "integer"))
```

`bl_find_local_cf()` rounds the counterfactual after the search and the fixed restore
(`fixed-constraint-restore.md` §5). In the harness (Steps 11 and 14):

| Feature | Observed | Search's point (after the restore) | Counterfactual | Sparse |
|---|---|---|---|---|
| `loan_amnt` | 12,000 | 12,000 | 12,000 (fixed) | 12,000 |
| `loan_int_rate` | 12.69 | 14.44 | 14.44 | 14.44 |
| `loan_percent_income` | 0.0500 | 0.0457 | 0.0457 | 0.0500 (not needed) |
| `credit_score` | 634 | 634.0454 | 635 | 634 ("No effect") |

Both whole neighbours of 634.0454, 634 and 635, still cross at 0.693 and meet `"increase"`;
635 is nearer the target in the Mahalanobis metric (√d²_M 0.6996 against 0.7031). The
prediction stays 0.693, `credit_score` has Shapley value 0 ("No effect"), and the sparse
counterfactual changes `loan_int_rate` alone, crossing at 0.512 (§2).

**Why it rounds this way.**

1. **Units:** it rounds in raw units, converting with `bl_result$scaling` when one is set, so
   `"integer"` means a whole credit point, not a whole SD. Rounding in model units would, for a
   standardised model, round to whole training SDs.
2. **Per feature:** it is set per feature, and only the features marked are rounded. One step for
   every feature could not suit features of different scales, such as a credit score, a loan
   amount and an interest rate.
3. **Direction:** both the whole value below and the one above are tried, and a candidate that
   breaks a direction constraint is dropped. A move back onto the observed value is kept only if
   it still crosses. Rounding to the nearest whole value alone could cancel a small move or
   reverse it.
4. **Constraints:** every `set_filters()` constraint, and `train_ranges` while that switch is on,
   is re-checked for each candidate.
5. **Validity:** each candidate is re-scored and kept only if it crosses; with none valid, no
   counterfactual is reported, with the reason. The Shapley values and the sparse counterfactual
   then start from a valid whole-number point, so rounding cannot invalidate the sparse
   counterfactual.

**What it does not do.** It rounds to whole numbers only. A different step, such as 0.5 or 100,
would need a per-feature step size in `set_filters()` (`future-improvements-review.md` §1.4.2,
"Other step sizes").
