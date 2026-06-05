# IMPLEMENTED 2026-05-29

> Fixes 1–5 and Fix 7 implemented and verified (77 PASS, 0 FAIL, 4 WARN).
> Fix 6 (non-0.5 cutoff warning) was removed from scope before implementation at user request.

---

# Plan: Usability & Bug Fixes

## Context

Audit identified several crashes and user-facing friction points in the package. Three are outright crashes (silent `NULL` returns or missing `NULL` checks that throw cryptic errors). The rest are interface inconsistencies.

---

## Fix 1 — `bl_assemble()` crashes when `bl_model = NULL` (`R/result.R`)

**Problem:** Lines unconditionally dereferenced `bl_model$model`, `bl_model$model_type`, etc. when `bl_model = NULL`.

**Change:** Guarded all model field extractions with `if (!is.null(bl_model)) ... else NULL`. Updated `@param bl_model` roxygen to state NULL is allowed.

---

## Fix 2 — `bl_build_result()` returns `invisible(NULL)` instead of stopping (`R/result.R`)

**Problem:** On validation failure, the function called `invisible(NULL)` with a message. Callers that stored the result got a silent `NULL`.

**Change:** Replaced two `invisible(NULL)` error-path returns with `stop(..., call. = FALSE)`. The `bl_model = NULL` notification path was deliberately left as `message()` (not an error).

---

## Fix 3 — `print.bl_points()` crashes when no model present (`R/project_points.R`)

**Problem:** Unconditionally accessed `pred_prob`, which is `NULL` when no model was supplied.

**Change:** Wrapped the pred_prob print lines in `if (!is.null(x$pred_prob))` with a `"(no model)"` fallback.

---

## Fix 4 — Typo: `"RForrest"` → `"RForest"` (multiple files)

**Problem:** Random Forest model type was spelled `"RForrest"` (double-r) throughout.

**Changes:**
- Bulk renamed in `R/model_fit.R`, `R/model_utils.R`, `R/predict_utils.R` using `replace_all`.
- Added backward-compat shim in `bl_fit_model()`, `bl_wrap_model()`, `.fit_model()`, `.pred_function()` that warns and redirects the old spelling.
- Updated `parsnip_types` in `print.bl_model()`.

---

## Fix 5 — Silent ignored hyperparameters in `bl_fit_model()` (`R/model_utils.R`)

**Problem:** Any key in `model_params` not recognised by the fitting branch was silently dropped.

**Change:** After `utils::modifyList()`, compare `names(model_params)` against `names(defaults)` and warn on unknowns with `sprintf()`.

---

## Fix 6 — REMOVED

Non-0.5 cutoff warning was removed from scope at user request before implementation.

---

## Fix 7 — Rounding discrepancy in `bl_predict()` docs (`R/project_points.R`)

**Problem:** `@return` said "rounded to 4 d.p." but code does 3 d.p.

**Change:** Updated `@return` description to say "3 decimal places".
