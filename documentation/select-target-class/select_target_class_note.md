# Which observation's class colours the local-biplot target?

Date created: 2026-10-29
Date last edited: 2026-11-16

A note-and-script bundle, following the `documentation/biplot-rotation/` convention.
`select_target_class_example.R` reproduces every number and the figure below. It
writes `select_target_class_results.txt` and `select_target_class_figure.png`.

## The question

`plot.bl_local_result()` colours the target by its confusion group: TP red,
TN blue, FP purple, FN orange. The colour compares the target's **true** class
with its **predicted** class at the cutoff.

`bl_select_target(bl_result, target = i, data = df)` accepts a row index into
any data frame, so the true class must come from the same row as the features:
`df[i, ]`, which need not be a row of `test_data`.

## The rule

`bl_select_target()` finds the true class where it chooses the row, and stores
it in the target as `bl_target$class_obs`.
- **The source row:** the class is read from the one-row source that supplies
  the features `x_obs`:
  - `data[i, ]` for an integer index `i` (`data` defaults to
    `bl_result$test_data`; `row_id` is `i`);
  - `target` itself for a one-row data frame (`row_id` is `NA`).
- **The column:** `class` when the source has one, otherwise the first column
  that is not a feature.
- **The value:** the label is read as an integer, so 0/1, or a factor with
  levels `"0"` and `"1"`, gives 0 or 1. A non-numeric label (such as `"yes"`)
  is unknown, `NA`; so is a source with no column besides the features.

| target | `class_obs` | `row_id` |
|---|---|---|
| `target = i` | class of `test_data[i, ]` | `i` |
| `target = i, data = train_data` | class of `train_data[i, ]` | `i` |
| `target = train_data[i, ]` (with its class column) | class of `train_data[i, ]` | `NA` |
| a one-row frame of the features only, such as an external applicant | `NA` | `NA` |

**What uses it:**
- **`plot.bl_local_result()`** colours the target from `class_obs` against
  `pred_class`:
  - a known class gives the confusion colour;
  - an unknown class (`NA`) gives the colour of the prediction: red for a
    predicted 1 and blue for a predicted 0 (the TP and TN colours);
  - a numeric label other than 0 or 1 matches no group and is drawn grey;
  - a `bl_target` without a `class_obs` field is coloured by prediction.
- **`plot.bl_sparse_result()`** draws the local biplot through
  `plot.bl_local_result()`, so its target has the same colour.
- **The training points** of the local biplot are coloured from their own rows
  of `train_data`.
- **No print method shows `class_obs`:** `print.bl_target()` gives the row, the
  predicted probability and class, and the feature values. The global biplot
  draws a target (`target_point =`) as a yellow highlight, whatever its class.

**Tests:** `tests/testthat/test-local_cf.R` (an index into `train_data`, the
same class from a one-row frame, the `test_data` default, `NA` without a class
column, and the plot of a `bl_target` without the field).

## What a lookup by row number would do, on Pima

Without this rule, a plot that read the class as
`bl_result$test_data$class[row_id]` would take another observation's class
whenever `data` is not `test_data`.
- **Wrong colour:** the target would be drawn in the wrong confusion colour
  whenever that class differs from its own. Nothing would warn, and the label
  would still show `i`, so the plot would look right.
- **Lost colour:** for `i > nrow(test_data)` the lookup would be out of range.
  The target would then lose its confusion colour and be drawn by prediction
  only.

The script measures this on the Paper 1 set-up:
- 392 complete rows with p = 7 (BloodPressure dropped);
- `train_fraction = 0.8`, `seed = 121L`, `hull_fraction = 1`;
- a GAM with `model_type = "GAM"`;
- CVA on the four prediction outcomes.

| quantity | value |
|---|---|
| training rows / test rows | 313 / 79 |
| rows compared (`i <= nrow(test_data)`) | 79 |
| class of `train_data[i]` differs from `test_data[i]` | **36 of 79** |
| training rows with `i > nrow(test_data)` (confusion colour lost) | 234 |

So with the 313 training rows as targets, a lookup by row number would draw:
- 36 targets in the wrong confusion colour;
- 234 without one;
- only 43 correctly.

**Worked case: training row 2,** the first training row with a mismatched class
whose local search finds a counterfactual.
- **Prediction:** probability 0.227, so class 0.
- **True class:** 0, so the target is a **TN** (blue), as `bl_select_target()`
  records it.
- **A lookup by row number:** it would read `test_data[2]`, which is class 1,
  and draw the target as an **FN** (orange).

![Left: the colour from test_data[2] (FN, orange). Right: the colour from the target's own class (TN, blue).](select_target_class_figure.png)

**How the figure is drawn.** Both panels come from one local search, and only
`bl_target$class_obs` differs between them:
- **Left:** `class_obs` is set to `test_data[2]`'s class, the colour a lookup by
  row number would give.
- **Right:** `class_obs` is the target's own class, as `bl_select_target()` sets
  it.

## Reproducing

From this folder:

```
Rscript select_target_class_example.R
```

That command uses the installed build. To measure the working tree, run
`devtools::load_all("<package root>")` first, then
`source("select_target_class_example.R", chdir = TRUE)`. Either way the script
writes its results and figure into this folder.
