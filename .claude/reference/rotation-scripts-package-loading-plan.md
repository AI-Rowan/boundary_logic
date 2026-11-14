> **IMPLEMENTED 2026-11-07.** The rotation scripts 17_1-17_3 load the repo's working tree inside the
> repo and check the installed version elsewhere; the package was reinstalled (0.2.0).

# Plan: make the rotation scripts load the right package, and reinstall it

Date created: 2026-11-07
Date last edited: 2026-11-07

> **Housekeeping (done 2026-11-07):** written in plan mode into
> `before-updating-the-paper-wiggly-wave.md` (the algorithms-document plan); step 0 moved it here and
> restored that file.

## Implementation notes (2026-11-07)

All steps are done, and nothing is committed.
- **The block** is identical in 17_1, 17_2 and 17_3. It was written once and inserted by script.
- **Version check:** tested against the installed 0.1.0 before the reinstall; it stopped with the
  message.
- **Reinstall:** this devtools rejects `upgrade = "never"`, so it ran as
  `devtools::install(upgrade = FALSE, build_vignettes = TRUE)`. The installed package is now 0.2.0,
  with `set_filters(hull =, train_ranges =)` and the vignette.
- **Runs:** each script `source()`d from `tempdir()`, and 17_3 also from the repo root. Every run
  printed "boundarylogic 0.2.0 (working tree)" and ran to the end.
- **Logs:** they equal the `load_all()` diagnostic runs, apart from 17_3's timing lines.

## Context

`documentation/biplot-rotation/17_1_loan-gam-pca-3d-3var-rotation.R` and `17_2_loan-gam-pca-4var-slices.R`
fail when run as they stand, for example in RStudio or with `Rscript`.

**The error.** Both start with `library(boundarylogic)`, which loads the *installed* package. That is
version 0.1.0, built on 2026-09-10 (`%LOCALAPPDATA%/R/win-library/4.6/boundarylogic`). Its
`set_filters(bl_target, ...)` has no `hull` switch, so the scripts' `set_filters(tgt, hull = FALSE)`
(17_1 line 524, 17_2 line 266) treats `hull` as a feature and stops:

> Unknown variable name(s) in constraints: hull. Available: …

The installed 0.1.0 also differs in ways that would change the numbers if it got further: `max_pairs = 10L`;
`distance` defaults to `"mahalanobis"`, with no `"auto"`; and its nearest-vertex tie-break is random (fixed
on 2026-10-30).

**Checked (2026-11-07):** against the repo's code (`load_all()`), all three scripts that use the package
(17_1, 17_2 and 17_3, which loads it the same way) run to the end, with every assertion passing. So the
scripts are correct; only their package loading is fragile. `17_1a` does not use the package.

**Wider effect.** The Paper 1 repo's scripts also call `library(boundarylogic)` (its `00_setup.R`), so
they get the same stale 0.1.0.

**User decision (2026-11-07):** fix the three scripts' loading, and reinstall the package.

## Design

**One loading block, the same in 17_1, 17_2 and 17_3**, replacing `library(boundarylogic)` (in 17_3,
the line inside `suppressPackageStartupMessages()`):

1. If the development copy is already loaded (`pkgload::is_dev_package("boundarylogic")`, as the
   CLAUDE.md driver does), keep it.
2. Otherwise, look for the package repo:
   - start from the script's own folder when it is known (`source()`'s `ofile`, `Rscript`'s `--file=`)
     and from `getwd()`;
   - walk up to a `DESCRIPTION` whose `Package` is `boundarylogic`;
   - if one is found, `pkgload::load_all()` it (quietly), so the scripts always describe the code
     beside them.
3. Otherwise `library(boundarylogic)`, the installed package.
4. **Version check:** stop unless `packageVersion("boundarylogic") >= "0.2.0"`, with a message that
   names the version found and the fix (`devtools::install()` from the repo, or run inside the repo).
5. Print one line: the version, and whether it is the working tree or installed.

The helpers stay inside `local()`, so `rm(list = ls())` and the scripts' own objects are unaffected. No
other line of the scripts changes.

**Reinstall:** `devtools::install(upgrade = "never")` from the repo's working tree. Its `R/` equals
`main` (`0689341`); `documentation/` is build-ignored.

## Steps

0. **Housekeeping:** restore the algorithms plan to this file, word for word, and save this plan as
   the plan-mode working copy (removed 2026-11-07).
1. **Scripts:** add the loading block to 17_1, 17_2 and 17_3; update the header line of 17_3 that names
   its requirements if needed.
2. **README** (`documentation/biplot-rotation/README.md`, "Running them"): the scripts load the
   repo's working tree themselves inside the repo; elsewhere they need the installed package at 0.2.0 or
   later, and stop with a message otherwise. Drop the now-redundant advice to `load_all()` first.
3. **Test the version check first**, while 0.1.0 is still installed: run the loading block on its own in
   a fresh R, with the working directory outside the repo and no script path, so it falls back to
   `library()`. It must stop with the message. Then **reinstall** the package
   (`devtools::install(upgrade = "never")`).
4. **CLAUDE.md §5:**
   - the "Running the `biplot-rotation/` scripts" bullet: they load the working tree themselves; the
     driver still works;
   - one new line: reinstall the package (`devtools::install()`) after each merge to `main`, because
     the Paper 1 repo and any script run outside the repo use `library(boundarylogic)`. It was 0.1.0
     until 2026-11-07.
5. **Verify** (below).
6. **Stop and report.** Do not commit.

## Error cases and risks

- **The script's own path is unknown** (lines run one at a time in RStudio): the `getwd()` search still
  finds the repo when the project (`boundary_logic.Rproj`) is open, since its working directory is the
  repo root. If neither is found, the installed package is used, and the version check stops with a
  clear message if it is too old.
- **`pkgload` missing:** it comes with `devtools`, which these workflows already use; without it the
  block falls back to `library()` and the version check.
- **Installing fails, for example because another R session has the package loaded and locks the DLL:**
  report it and ask the user to close the other session.
- **Numbers:** the loading change alters nothing in the computation (the same code is loaded), so the
  note's numbers stand; step 5 checks it.

## Verification

1. **Run from inside the repo:** each script run through `source()` from the repo root, with no
   `load_all()` beforehand, prints "boundarylogic 0.2.0 (working tree)" and runs to the end. Its console
   output must equal today's diagnostic logs (`load_all()` + `source()`), apart from timings.
2. **Run from elsewhere:** the same scripts `source()`d with the working directory set to `tempdir()`
   find the repo from their own path.
3. **The installed copy:** after the reinstall, `library(boundarylogic)` in a fresh R gives 0.2.0, and
   `args(set_filters)` shows `hull` and `train_ranges`.
4. **The version check:** tested in step 3, against the real 0.1.0, before the reinstall.
5. **Nothing else changes:** `git status` shows only the three scripts, the README and CLAUDE.md.
