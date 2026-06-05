# IMPLEMENTED 2026-05-29

> All remaining gaps resolved. devtools::document() regenerated clean .Rd files.
> 77 PASS, 0 FAIL, 4 WARN.

---

# Plan: Documentation Gaps

## Context

Audit found ~11 exported functions had roxygen comments but no generated `.Rd` files. Several existing `.Rd` files had minor inaccuracies. Three orphaned `rlang` imports were suspected.

---

## Step 1 — Fix known roxygen inaccuracies

### 1a — `bl_fit_model()`: stale `rounding` param
Already resolved before this plan was implemented.

### 1b — `bl_predict()`: `@return` 4 d.p. → 3 d.p.
Fixed as Fix 7 of usability-bug-fixes plan. Already done.

### 1c — `bl_surrogate()`: missing hull coverage fields
Added `n_total`, `n_in_hull`, `hull_coverage` `\item{}` entries to `@return` `\describe{}` block.
Backtick-wrapped `[0, 1]` to prevent roxygen2 treating it as a cross-reference link.

### 1d — `bl_build_projection()`: missing `condition_number` field
Added `condition_number` `\item{}` to `@return` `\describe{}` block with 1e10 threshold description.

---

## Step 2 — Regenerate all `.Rd` files

`devtools::document()` run. All 21 exports confirmed to have `.Rd` files.

---

## Step 3 — Orphaned `rlang` imports

Confirmed: NAMESPACE and DESCRIPTION contained no `rlang` references. No action needed.

---

## Affected `.Rd` files

- `man/bl_surrogate.Rd` — regenerated with hull coverage fields
- `man/bl_build_projection.Rd` — regenerated with `condition_number` field
