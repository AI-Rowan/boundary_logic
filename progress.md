# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Standing items / carry-overs

- **Mac tester confirmation pending** — once macOS install is confirmed, merge
  `method_developments` -> `main`.
- **`scripts/03_loan_status_Boundary_Logic.R` working-edit committed** (`f540607`, 2026-06-16):
  `person_age` dropped from the reduced (v2) model, `tdp <- 13`, v2 boundary plot
  `type = "boxplot"`. Committed as-is at user request — not pushed. Known leftovers the user
  chose to keep: Steps 8b/9/15 still pass `person_age` to plot label/tick vectors (harmless
  `variable not found` warnings at plot time), and Step 14 retains a `new_title = "xx"`
  placeholder. Commit used the auto-configured git identity
  (`Rowan <adriaan.rowan@yala.co.za>`); no `user.name`/`user.email` is set.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 136 PASS / 0 FAIL / 4 WARN** (the 4 WARN are pre-existing biplotEZ
  CVA 2-class notices).
- **Verification command:**
  ```r
  "/c/Program Files/R/R-4.6.0/bin/Rscript" -e "devtools::test()"
  ```

---

## Session summary (2026-06-16, latest) — raw-unit scaling support

Three linked deliverables let a user who standardised features before model fitting read the
biplot and counterfactuals in original units. All display-only: stored objects/geometry/model
stay in model (standardised) units. Test baseline 77 -> 136 PASS / 0 FAIL.

1. **`bl_set_scaling()` + raw-unit biplot axes.** New exported `bl_set_scaling(x, center, scale,
   method)` (`R/data_prepare.R`, with private `.align_scaling_vec()`) records the per-feature
   transform `std = (raw - center)/scale` as a `scaling` field on `bl_data`/`bl_filter_result`,
   threaded through `bl_filter_outliers()` -> `bl_assemble()` -> `bl_result$scaling` (and onto
   `bl_projection$scaling`). New private `.bl_rescale_biplot_axes()` (`R/plot_biplot.R`)
   relabels biplotEZ axes into raw units by transforming only `$means`/`$sd` (geometry
   byte-identical); wired into all biplot plot methods (`plot_biplotEZ`, `plot.bl_projection`,
   `plot.bl_local_result`, `plot.bl_surrogate`, `plot.bl_sparse_result`). See
   `documentation/scaling_axis_relabel_note.md` + impl summary §4.4.1.

2. **Raw-unit feature values in Shapley/sparse/target displays.** New private `.scale_to_raw()`
   (`R/data_prepare.R`; `kind="level"` -> `v*scale+center`, `kind="delta"` -> `v*scale`).
   `plot.bl_shapley()` label becomes observed -> counterfactual (raw); `print.bl_shapley()`,
   `print.bl_sparse_result()`, `print.bl_target()` convert (scaling stored on `bl_target`).
   Shapley *contributions* stay in prediction-impact units (not converted). Impl summary §4.9.

3. **`bl_local$bl_counterfactual`.** `bl_find_local_cf()` now returns a `"bl_counterfactual"`
   object (private `.make_bl_counterfactual()`) mirroring `bl_target`: stores model-unit `x_cf`
   (= `B_x`), prints in raw units via `.scale_to_raw()`. `print.bl_local_result()` shows Target
   and Counterfactual blocks together. Impl summary §4.8.

**Script 03**: v2 block now standardises on the **training** rows only and applies those stats
to the test split (was: scaled on full data -> test leakage). Split happens first on raw data.

**Docs**: impl summary §4.4.1/§4.8/§4.9; both review reference docs; technical note;
CLAUDE.md (file map, S5 axis-relabel Always-Do, S9 entries flipped to IMPLEMENTED + new
per-variable-methods future item); memory sync (`reference_data_prep_functions.md`); three
plans archived to `.claude/reference/` with IMPLEMENTED banners.

**Committed + pushed**: `00b43a5` on `method_developments` (24 files, +1515/-28), pushed
`10a81f3..00b43a5` (also carries the earlier `f540607` person_age script edit). Commit used the
auto-configured identity `Rowan <adriaan.rowan@yala.co.za>` (no `user.name`/`user.email` set).

### Dead ends / corrections this session

1. **CVA `scaled=FALSE` worry was a false alarm.** During planning I feared the `$means`/`$sd`
   override wouldn't relabel CVA axes correctly (CVA stores `sd = 1` and skips the sd term in
   biplotEZ's `Xhat` reconstruction). Algebra showed the transform is in fact uniform across
   PCA(T/F) and CVA (`sd_new = scale` supplies the missing `d(raw)/d(std)` factor); confirmed
   empirically (geometry byte-identical, labels on raw scale). No special-casing needed.
2. **Test used a named-vector subset and hit the wrong validation branch.** `d$center[1:2]` is
   still *named*, so `.align_scaling_vec()` took the name-match path and raised "missing entries"
   instead of the expected "length" error. Fix: test with `unname(d$center)[1:2]`.
3. **`.scale_to_raw()` strips names** (it does `unname(values) * ...`). `print.bl_target()` /
   `print.bl_counterfactual()` must re-attach `names()` after calling it, or the printed vector
   loses its variable labels. (The NULL-scaling path returns the value unchanged *with* names, so
   the bug only showed when scaling was set.)
4. **Value-reporting was initially deferred, then pulled back in.** The first round shipped only
   raw *axes* and explicitly deferred raw *values* (Shapley/sparse/target) as CLAUDE.md S9. The
   user then asked for the values too, so it was implemented in a second round — Shapley
   *contributions* stay in prediction-impact units (not feature units, not convertible).
5. **Script 03 standardised on the full dataset (test leakage)** — caught by the user. Fixed to
   fit the scaling on the (filtered) training rows only and apply those stats to the test split;
   the split now happens first on raw data (the hull filter self-standardises, so filtering on
   raw vs standardised is equivalent).

### Architecture decisions / conventions

- **Display-only raw-unit conversion.** Computation, storage, projection geometry, and the model
  all stay in model (standardised) units; only plot/print methods convert, via
  `.bl_rescale_biplot_axes()` (axes) and `.scale_to_raw()` (values). Gated on
  `bl_result$scaling`; no-op for older objects. Added as CLAUDE.md S5 conventions.
- **Level vs delta is the key correctness rule** for `.scale_to_raw()`: a *level* (observed, CF,
  sparse) maps `v*scale + center`; a *delta* (Shapley `data_to_boundary`) maps `v*scale` (no
  centre). Any future feature-value display must pick the right `kind`.
- **`scaling` is a separate field** from `X_center`/`X_sd` (internal PCA standardisation) and the
  `standardise` flag — never conflate them.

### Next steps

1. **Mac-tester confirmation, then merge `method_developments` -> `main`** (standing carry-over).
2. **Optional follow-ups** (only on request): per-variable transform families (non-affine needs
   spline-calibrated axes — CLAUDE.md S9); a convenience accessor if reading `B_x`/CF in raw
   units programmatically becomes common.
3. The untracked `scripts/16_*`/`18_*` are the user's own files — confirm with the user whether
   they should be tracked; not touched by this work.

---

## Session summary (2026-06-16) — `new_title` for all biplot plot functions

### Completed this session

1. **Extended the `new_title` parameter to every biplot `plot()` method.** Previously only
   `plot_biplotEZ()` accepted `new_title` (default `NA`, overriding
   `bl_result$biplot_obj$Title`). A user calling `plot(bl_sparse, new_title = "xx")` had the
   argument silently swallowed by `...`. Now uniform across all four biplot classes:
   - `plot.bl_local_result()` (`R/local_cf.R`) — **code change**: added `new_title = NA`
     param + roxygen; after the hardcoded `sprintf("Local biplot -- target N [pair ...]")`
     title, added `if (!is.na(new_title)) biplot_plot$Title <- new_title` (override the
     auto-generated default only when supplied).
   - `plot.bl_surrogate()` (`R/surrogate.R`) — **code change**: added `new_title = NA` param
     + roxygen; applied `if (!is.na(new_title)) biplot_obj$Title <- new_title` right after
     `biplot_obj <- bl_result$biplot_obj`, before the rotation call.
   - `plot.bl_sparse_result()` (`R/shapley.R`) — **no code change**: forwards `...` to
     `plot.bl_local_result()`, so it works automatically once the above landed. Only the
     roxygen `@param ...` example list was updated to mention `new_title`.
   - `plot.bl_result()` (`R/plot_biplot.R`) — **no change**: already
     `function(x, ...) plot_biplotEZ(x, ...)`; `new_title` flows straight through.

2. **Reference doc** `.claude/reference/review_section9_to_15.md` updated (per the CLAUDE.md
   Section 2 update-trigger rule — it covers `local_cf.R` + `shapley.R`): Step 12 documents
   `new_title` for the local biplot, Step 14 documents the sparse pass-through, plus a
   **base-graphics title-error note** (see dead-end #1). Not mirrored in auto-memory, so no
   memory sync needed.

3. **Verification (all green):**
   - `devtools::document()` — clean; regenerated `man/plot.bl_local_result.Rd`,
     `man/plot.bl_sparse_result.Rd`, `man/plot.bl_surrogate.Rd`.
   - `devtools::test()` — **77 PASS, 0 FAIL, 4 WARN** (baseline unchanged; additive change).
   - Visual check via throwaway `Rscript verify_title.R` (deleted after use, per the
     biplotEZ Windows-segfault rule): all four `plot(..., new_title = ...)` calls render the
     custom title; omitting it reproduces the existing defaults.

4. **Committed and pushed** — commit `dc0c3c4` on `method_developments`, pushed to
   `origin` (`d19a331..dc0c3c4`). 7 files (3 `R/`, 3 `man/`, 1 reference doc).

### Dead ends / corrections this session

1. **The plan's predicted title-error was wrong; corrected after empirical check.** The plan
   asserted that passing a list/data.frame as `new_title` would raise `"invalid 'main'
   argument"`. Direct testing of `graphics::title(main = ...)` showed it is far more tolerant:
   numeric, multi-element character vectors, lists, and data.frames are all coerced via
   `as.character()` and do **not** error. The only failure mode is a value base R cannot
   coerce to a character vector — a **function/closure** or **environment** — which raises
   `"cannot coerce type 'closure' to vector of type 'character'"` from inside the biplotEZ
   `plot()` flush, not from the plot function's own arg handling (no validation is performed,
   matching `plot_biplotEZ()`). Both the reference doc and the archived plan were corrected to
   reflect the verified behaviour.

2. **`Rscript -e` quoting failed in Bash on Windows** when testing the `title()` coercion
   one-liner (the inner double-quotes broke the cmd-level parse: `'num:" , chk(3.14)' is not
   recognized`). Same class of issue noted in earlier sessions. Fix: wrote the probe to
   `chk_title.R` and ran `Rscript chk_title.R`. Reinforces the existing "write a `.R` file"
   rule — it applies to any non-trivial `-e` payload on Windows, not just biplotEZ calls.

### Architecture decisions / new conventions

- **`new_title = NA` is now the uniform convention across all biplot plot methods.** No new
  pattern was invented — the existing `plot_biplotEZ()` contract (`new_title = NA` default;
  `if (!is.na(new_title)) <biplot_obj>$Title <- new_title`; biplotEZ renders `$Title`) was
  simply propagated. Added as a one-line "Always Do" bullet in CLAUDE.md Section 5 so any
  future biplot plot function includes it.

### Next steps

1. **User to commit `scripts/03_loan_status_Boundary_Logic.R`** when satisfied with the
   experimental edits (or revert the `new_title = "xx"` placeholder to a real title first).
2. Carried over: Mac tester confirmation, then merge `method_developments` -> `main`.
3. Outstanding deferred work unchanged — see CLAUDE.md Section 9 (Phase 3 unit tests, full
   Mahalanobis Shapley in Phase 2, CRAN prep, etc.).

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

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
