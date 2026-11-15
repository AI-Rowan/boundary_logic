> **IMPLEMENTED 2026-11-15.** Archived from the plan-mode working copy. Session 13 closed the
> outstanding items apart from the Paper 1 plan update: the methodology note stored as Markdown
> (`methodology-changes.md`; the `.qmd`/`.docx` removed), `documentation/future-developments.md`
> (CLAUDE.md §9 slimmed to a pointer), the current-code rule (CLAUDE.md §2) applied to
> `local-convex-hull.md`, the root documents refreshed, the raw-unit label test (578 PASS),
> `issues.csv` moved to `issues.md`, and the hand-over Claude Doc. The same session's two earlier
> plans (the Pima external patients for Finding F4; `algorithms.md` describing the current code only)
> were overwritten in the working copy; `progress.md`'s session-13 summary records them.

# Close the outstanding items (apart from the Paper 1 plan update)

Date created: 2026-11-15
Date last edited: 2026-11-15

## Context

The user asked for the outstanding items apart from the Paper 1 plan update (`progress.md` item 1),
with ways to close them, and has decided each one (2026-11-15):
- **Methodology note:** re-render it. The supervisor is happy with the changes, and the note is to
  be stored as the documentation of the changes (Gate 3 passed). It is to be **stored as a plain
  `.md`** that needs neither the `.qmd` nor the `.docx`, so both can be removed.
- **`scripts/paper1-pima/issues.csv`:** move it to an `.md` file.
- **Hand-over:** do it now, as a **Claude Doc**.
- **Algorithm boxes:** straight after the Paper 1 plan update, with the plan's stale lines refreshed
  then (not now).
- **The current-code-only rule:** make it a standing rule for the technical notes and implement it:
  record the rule now, and sweep one note per session, starting with `local-convex-hull.md`.
- **Future work:** none of the CLAUDE.md §9 value tests or other deferred work is for Paper 1.
  Instead, create **`documentation/future-developments.md`**, a summary of every future development
  (code and method) for later projects. **CLAUDE.md §9 is slimmed to a pointer.**

The Paper 1 plan update itself (with T3/T4 and the `fig-ext` decision) is excluded.

## Steps

### 0. Start the long re-render runs in the background

They take about 25 minutes in total, so they run while the documentation steps proceed:
- From `scripts/paper1-pima/`: `Rscript 04_method-change-effects.R` (about 10 min), then
  `Rscript 05_full-space-gap.R` (about 15 min).
- `numbers.rds` and `findings-effects.csv` were already written today at `HEAD` `e41f22f`, so
  `run_all.R` and `06` are not re-run.
- R/ is unchanged, so the note will name `e41f22f`, the measured code, which equals tag `v0.2.0`'s R
  code.
- No R script is edited while they run (CLAUDE.md §5).

### 1. The rule (C1) in CLAUDE.md §2, "Documentation folder convention"

One line:
- Technical notes in `documentation/` describe the current code and method in the present tense.
- The package's own history (earlier versions, commits, dates, decisions, plans) goes only in the
  methodology note and `NEWS.md`.
- An alternative is written in the conditional ("without this rule, … would …"), with numbers
  measured by removing that rule alone.
- Thesis citations stay. Research-code citations stay only where the note's subject is the original
  method (rotation note §2).
- Walkthroughs and `algorithms.md` already follow it.

### 2. Sweep `documentation/local-convex-hull.md` (C1, first note; also closes item 13)

- Rewrite it to the rule: §2, "The package before 2026-11-05", every "since 2026-…", the "Design
  choices (… 2026-11-05; … revised on 2026-11-06, §7)" line, and §6's "Checks (2026-11-06)".
- Repoint each **"§7"** reference (introduction, §3, §6) to what already describes the rules: §3's
  table rows, and `algorithms.md` F2 and F4. No commit ever had a §7.
- §2's content becomes "what a search without the hull does" (`set_filters(hull = FALSE)`), not
  "what the package did".
- Every number already stated as an `hull = FALSE` result stays. A number that compared two
  versions differing in more than one rule is re-measured with that one rule removed (the
  `patches.R` drivers in the scratchpad), or left out and reported to the user.
- The `progress.md` list of the remaining notes, one per session: rotation note, the tie-break
  note, `sparse-counterfactual.md`, `boundary-distance-and-robustness.md`,
  `fixed-constraint-restore.md`, `select-target-class` note, `biplot-rotation/README.md`,
  `scaling_axis_relabel_note.md`.

### 3. Root documents (item 11), in the present tense

- **`2 implementation_summary.txt`:** the "Last updated" date; the test tree (15 `test-*.R` and
  `helper.R`); the documentation tree (add `local-convex-hull.md`, `algorithms.md`,
  `methodology-changes/`, `future-developments.md`; drop the deleted presentation PDF); the
  expectation count from step 4's real run; and any "findings and decisions" wording now called
  "behaviours of the code".
- **`1 Foundation intro documents.txt`:** the "last revised" date; §8's file tree, rebuilt from `git
  ls-files`; §9's branches (add `paper_nr1_changes` and tag `v0.2.0`).

### 4. Test gap (item 12)

- One `test_that()` in `tests/testthat/test-set_scaling.R`: an iris PCA with `standardise = FALSE`
  and a `bl_set_scaling()`. Call `.bl_rescale_biplot_axes()` on the biplotEZ object and assert the
  exact labels (`std * scale + center`) at matching ticks.
- Then run `devtools::test()` (expected: 561 + n PASS, 0 FAIL).
- First read `R/plot_biplot.R`'s helper and the existing relabel test, to reuse their set-up.

### 5. `documentation/future-developments.md` (new) and CLAUDE.md §9

Present tense, grouped as **Method**, **Code and API**, **Tests**, **Documentation and release**.
Each item gives:
- what it is and why;
- the functions affected;
- its expected effect on results;
- the effort;
- the source.

Sources to merge, without duplicates:
- the methodology note's "Future work" (Table 5, 12 items) and its text;
- CLAUDE.md §9's open entries: the mathematical-method vignette; full-Shapley Phase 2; per-variable
  transformations; the rotation concerns and the §5 deterministic rotation (with the `.bl_rotate()`
  under-determination facts); the Phase 3 and `bl_build_result()` value tests; LDA multi-class; CRAN;
- `future-improvements-review.md`'s design-document-only items: B5 (the CVA identity test), B8
  (`pairs =`), B11, B21 (value tests), B22 and B23, plus its other "future improvement" items;
- the notes' own future sections: rotation note §5, §6.2 and §6.3; `fixed-constraint-restore.md` §8;
  `boundary-distance-and-robustness.md` §7;
- `sparse-counterfactual.md` (rounding resolutions other than 1);
- the implementation summary's LDA design note.

**CLAUDE.md §9** becomes:
- the rule "do not implement these without a new instruction";
- the `.bl_rotate()` warning (do not "fix" the under-determined plane or implement §5 without
  instruction; every Phase 3 number moves);
- a pointer to `documentation/future-developments.md`.

The struck-through resolved entries and the item list move out. Add `future-developments.md` to
CLAUDE.md §2's reference list.

### 6. Freeze the methodology note as Markdown (B1), after step 0

1. **Render once more.** From `documentation/methodology-changes/`, run `quarto.exe render
   methodology-changes.qmd --to gfm -o methodology-changes.md`, with `QUARTO_R` set. This writes
   all 126 inline values as text. The claims chunk must pass on this last render; if a claim fails,
   report it rather than change the check. The note has no figures.
2. **Clean up the `.md`:**
   - turn the YAML title, subtitle, author and date into a title block;
   - keep the table of contents;
   - keep the section numbers in the headings, because the text cites "Section 3.12" and so on
     (add them by hand if gfm drops them);
   - check that the tables, their captions ("Table 5. …") and the `$…$` math read correctly;
   - remove any leftover Quarto `<div>` wrappers.
3. **Add a header line** stating that this is the record of the changes up to `boundarylogic` 0.2.0
   (tag `v0.2.0`), measured at commit `e41f22f`. The values come from `method-change-effects.csv`
   and `findings-effects.csv` in this folder and from the Pima example's `full_space_gap.rds`.
4. **Reproduction section:** replace step 5 (`quarto render …`) with "the values in this note are
   those CSVs' at the commit above". The note is now edited by hand.
5. **Remove the old files:** `git rm` `methodology-changes.qmd` and `methodology-changes.docx`
   (both are tracked, so git keeps them), and any `_files`, `.gitignore` or `mediabag` leftovers.
6. **Keep the CSVs.** They are the raw measurements, and `review.qmd` reads `findings-effects.csv`.
7. **Repoint every reference** to the new `methodology-changes.md`:
   - CLAUDE.md §2 (".qmd rendered to .docx" becomes ".md") and the §5 bullet "Updating the
     methodology-changes note" (no render step: run `04`/`05`/`06` for the CSVs, then edit the
     `.md` by hand, quoting the CSV values);
   - `boundary-distance-and-robustness.md`, `fixed-constraint-restore.md`,
     `future-improvements-review.md`;
   - the comments in `04`, `05` and `06`;
   - `progress.md`;
   - the documentation tree in the implementation summary (step 3).
   - Leave `.claude/reference/` archives and the generated `docs/` site as they are.

### 6b. `issues.csv` to `scripts/paper1-pima/issues.md`

1. **Write `issues.md`** from the 13 rows: a heading per issue (`id`, short title), under "Open"
   (5, 9, 11, 12, 13) and "Closed" (2, 3, 4, 6, 7, 8, 10, 14). Under each heading, the fields
   symptom, where, effect on the paper, proposed fix and decision, verbatim.
2. **Check it** with a script that compares every CSV field with the `.md` text.
3. **Update `review.qmd`:** drop the `read.csv("issues.csv")`, replace the open-issues table with
   `{{< include issues.md >}}` under a section titled "Issues (open and closed)", and reword its
   prose (line 564).
4. Mention `issues.md` instead in `pruned-model-issues-review.md`, CLAUDE.md and `progress.md`.
5. **Remove `issues.csv`.** The folder is gitignored, so there is no git copy: keep one in the
   scratchpad first.
6. Re-render `review.qmd` and check the issues section.

### 7. The hand-over Claude Doc (B2)

Written for the Paper 1 writing work (the author, and anyone they share it with). First load the
`anthropic-skills:docs` skill, then create the doc. Its contents:
- the install line `devtools::install_github("AI-Rowan/boundary_logic@v0.2.0")`;
- v0.2.0: 561 PASS and `R CMD check` 0/0/0;
- where the outputs are: `scripts/paper1-pima/review.html`, `output/numbers.rds`, `table1.csv`,
  the figures;
- Table 1's runs (1, 3, 4, 5, 6, 7);
- the review-only figure `fig-ext` (whether it goes in the paper is still open);
- `documentation/methodology-changes/methodology-changes.md`;
- the open issues, in `scripts/paper1-pima/issues.md`;
- `future-developments.md`;
- how to reproduce;
- the documentation commit SHA, added once the user commits.

Its link goes in `progress.md`.

### 8. Close out

- **`progress.md`:**
  - items 3, 4, 11, 12 and 13 closed;
  - the note-sweep list (step 2);
  - item 2 kept (straight after the Paper 1 plan update, with its plan's stale lines refreshed
    then);
  - the session summary;
  - the test count.
- Update the project memory's status line.
- **Stop and report.** Commit only when the user says "commit". The commit would then hold every
  step's files: `paper_nr1_changes` pushed, `main` and `method_developments` fast-forwarded. The
  changes are documentation and one test, so no reinstall and no new tag.

## Error cases

- **04 or 05 stops on a self-check, or the note's claims chunk fails on the last render:** report it,
  and do not edit a check. Do not remove the `.qmd` until the `.md` is complete.
- **gfm output loses the section numbers or the table captions:** add them by hand in the `.md`,
  checking each against the `.qmd`.
- **The `.docx` is open in Word,** so the `git rm` fails: ask the user to close it.
- **`issues.md` differs from the CSV in any field:** fix it before removing the CSV.
- **A sweep number cannot be re-measured as a single-rule effect:** leave it out of the rewritten
  note and list it in the report.
- **The new test fails:** that means a real relabel bug. Report it; do not change package code
  without approval.
- **CRLF:** use the Edit tool.

## Verification

- `devtools::test()`: 0 FAIL, with the PASS count recorded.
- A grep of `local-convex-hull.md` for `§7|since 20|before 20|ef014ca` finds nothing.
- Every path in the root-document trees exists.
- `methodology-changes.md` names `e41f22f`, keeps its numbered sections and its Tables (Table 5 and
  so on), and contains no inline-R left over (`` `r `` absent). A grep of the repo (excluding
  `.claude/reference/` and `docs/`) for `methodology-changes.qmd` or `.docx` finds nothing.
- `issues.md` holds all 13 issues with every field (checked by script), and the review page's
  issues section renders from it.
- The Claude Doc opens, and its link is in `progress.md`.
- CLAUDE.md §9 is reduced to the rule, the warning and the pointer, and every item it dropped
  appears in `future-developments.md` (checked item by item).
