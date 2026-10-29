# Which observation's class colours the local-biplot target?

A note-and-script bundle, following the `documentation/biplot-rotation/` convention.
`select_target_class_example.R` reproduces every number and the figure below. It
writes `select_target_class_results.txt` and `select_target_class_figure.png`.

## The question

`plot.bl_local_result()` colours the target by its confusion group: TP red,
TN blue, FP purple, FN orange. The colour compares the target's **true** class
with its **predicted** class at the cutoff. When the true class is unknown, the
target is coloured by its predicted class only.

`bl_select_target(bl_result, target = i, data = df)` accepts a row index into
any data frame, so the true class must come from the same row as the features:
`df[i, ]`.

## The fault (package at `b7aee65`)

`bl_select_target()` stored only `row_id`. `plot.bl_local_result()` then read
the true class as `bl_result$test_data$class[row_id]`. When `data` is not
`test_data`, that is another observation's class.
- **Wrong colour:** the target is drawn in the wrong confusion colour whenever
  that class differs from its own. Nothing warns, and the label still shows `i`,
  so the plot looks right.
- **Lost colour:** for `i > nrow(test_data)` the lookup is out of range. The
  target then loses its confusion colour and is drawn by prediction only.

## Evidence on Pima (Paper 1 set-up)

The set-up:
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

So among the 313 training rows:
- 36 targets were drawn in the wrong confusion colour;
- 234 were drawn without one;
- only 43 were drawn correctly.

The loan data shows the same fault: with 221 training rows (`train_fraction = 0.01`),
the class differed in 120 of 221 rows.

**Worked case: training row 2.**
- **Prediction:** probability 0.227, so class 0.
- **True class:** 0, so the target is a **TN**.
- **The old lookup:** it read `test_data[2]`, which is class 1, and drew the target as
  an **FN** (orange).

![Left: the colour from test_data[2] (FN, orange). Right: the colour from the target's own class (TN, blue).](select_target_class_figure.png)

**How the figure is drawn.** Both panels come from one local search, and only
`bl_target$class_obs` differs between them:
- **Left:** `class_obs` is set to `test_data[2]`'s class, which reproduces the
  colour the old lookup drew.
- **Right:** `class_obs` is the target's own class.

**Checked against the code before the fix.** The same script was run against
`b7aee65`. Both panels were drawn orange, because that code ignored any class
the target carried and always read `test_data[2]`. That confirms the left panel
matches what the old plot drew.

## The fix

`bl_select_target()` resolves the true class where the row is chosen. It
stores the class as `bl_target$class_obs`, an integer, or `NA` when the source
has no class column. The one-row source is:
- `data[i, ]` for an index;
- `target` itself for a one-row data frame.

The column is `class` when present, otherwise the first non-feature column, and
a non-numeric label counts as unknown.

`plot.bl_local_result()` reads `bl_target$class_obs`. A `bl_target` built before
the field existed has none, and is coloured by prediction.

**The results:**
- **An index into `test_data`:** colours as before.
- **An index into any other frame:** uses that row's own class.
- **A one-row data-frame target that carries its class column:** now gets its
  confusion colour.
- **An external target without a class column:** is coloured by prediction.

**Tests:** `tests/testthat/test-local_cf.R`.

## Reproducing

From this folder:

```
Rscript select_target_class_example.R
```

That command uses the installed build. To measure the working tree, run
`devtools::load_all("<package root>")` first, then
`source("select_target_class_example.R", chdir = TRUE)`.
