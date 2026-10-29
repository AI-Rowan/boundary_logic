# Progress

This is the session-handoff doc: current state, decisions, test status, blockers, next step.
Full per-session detail lives in git history (`git log -p progress.md`); durable conventions
live in `CLAUDE.md`. Older sessions are condensed to one line each under **History** below.

---

## Standing items / carry-overs

- **`main` was merged and pushed on 2026-08-03** (see session below) — the Mac-tester gate was
  consciously overridden to get the latest code (incl. the loan-default vignette) publicly
  installable. No further merge is pending unless new `method_developments` work accumulates.
- **`scripts/` is now gitignored** (`.gitignore`) — it is the user's personal scratch playground
  and is never published. `17_*`, `18_*`, `LendingClub_PD_Model_R.R`, `loans_full_schema.csv`,
  `01_iris_SVM_PCA_biplot.R`, and `pima diabetes.csv` all live there, untracked, on disk only.
  The former main dev/testing harness (`03_loan_status_Boundary_Logic.R`) and its three
  `review_section*.md` walkthroughs were moved to a new tracked `development/` folder (excluded
  from the built package via `.Rbuildignore`) — see session below.
- **`.bl_rotate()` plane under-determination** — found 2026-07-15, **documented and deliberately
  not actioned**. See `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md` §4.1,
  CLAUDE.md §9. A candidate deterministic replacement is written up as a proposal in the same
  note's §6 (partial verification harness alongside it) — **not implemented**. Do not "fix" or
  implement §6 without an explicit instruction.
- **Deferred / future work** — see `CLAUDE.md` Section 9 (not duplicated here).
- **Baseline test status: 136 PASS / 0 FAIL / 4 WARN** (the 4 WARN are pre-existing biplotEZ
  CVA 2-class notices).
- **Verification command:**
  ```r
  "/c/Program Files/R/R-4.6.0/bin/Rscript" -e "devtools::test()"
  ```
- **`.claude/plans/` is gitignored** (`.gitignore:16`) — plan files are local-only and never
  committed. `.claude/reference/` *is* tracked (generic docs + archived plans only, now that the
  script-specific walkthroughs moved to `development/reference/`).

---

## Session summary (2026-09-10, latest) — integrated the `documentation/biplot-rotation/` bundle; fixed what the deletions broke

Goal: the user added `documentation/biplot-rotation/` — a much larger, consolidated technical note
plus the four R scripts that reproduce every measurement in it — snapshotted from a separate
writing project (`Explain XAI code`), and deleted three files this repo had built up on the same
topic: the PhD thesis `.docx`, an untracked draft `biplot_rotation_and_pair_selection.md` (never
committed, so no git loss), and the tracked `rotation_plane_underdetermination_note.md` (staged as
a `git mv` into `development/` by the 2026-09-01 session but never committed, so also gone from the
working tree). This is **documentation housekeeping only** — no `R/` change, nothing from the new
note's §6 proposal implemented.

### What changed

- **Ownership.** All four scripts (`17_1a_toy-2d-to-1d-rotation.R`,
  `17_1_loan-gam-pca-3d-3var-rotation.R`, `17_2_loan-gam-pca-4var-slices.R`,
  `17_3_minimal-rotation-verification.R`) and `README.md` carried a "REFERENCE COPY... the sandbox
  copy is the source of truth" banner; replaced with a short origin note (`Explain XAI code`
  sandbox, snapshot 2026-09-10) recording that this repo is now their home.
- **Nothing unique lost to the deletions.** Compared the deleted `rotation_plane_underdetermination_note.md`
  against the new note's §4.1 (using this session's own prior knowledge of the deleted file's
  content — no `git show`/recovery command was used, per explicit instruction). Everything
  substantive was already reproduced in §4.1, except two scope-limiting bullets ("Phases 1 and 2
  never call `.bl_rotate()`"; "`bl_set_scaling()` is display-only") — folded those into §4.1's "What
  is *not* affected" paragraph, then dropped the note's now-redundant citations of the deleted file
  (its source-material paragraph and §7 References).
- **Line-citation check.** Verified every `R/local_cf.R` line range the note cites against the file
  as it exists now: 5 of 6 exact; one (`432-621`, describing the per-pair loop inside
  `bl_find_local_cf()`) was off by one line — corrected to `432-622` (the loop's true bounds).
- **Repointed every stale reference** the deletions left in `CLAUDE.md` (§2 reading list, §6 new
  bullet below, §9 deferred-work entry, and reverted a documentation-folder-convention table edit
  from 2026-09-01 that no longer applies now that neither file lives in `development/`),
  `2 implementation_summary.txt` §4.2's caveat, and this file's Standing items — all now point at
  `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md`.
- **Surfaced two new findings in `CLAUDE.md` §6 Key Architectural Facts and §9 Deferred Work:**
  the note's §5.9 (`bl_find_local_cf()`'s reported `dist_mahalanobis` depends on prior calls in the
  same session — quote a counterfactual from a first call on a fresh `bl_result`), and its §6 (a
  written-up, **not implemented**, closed-form deterministic-rotation + reproducible-plane-search
  proposal, with a partial acceptance-test harness in `17_3_minimal-rotation-verification.R`).
- **Made the scripts runnable here**: added a short "Running these in the `boundary_logic` package
  repo" note to `README.md` (`devtools::load_all()` or install the branch under test first).

**Test status: unchanged, 136 PASS / 0 FAIL / 4 WARN baseline — not re-run this session (no `R/`
file touched; confirmed `git diff --stat -- R/` is empty).**

### Decisions

- **The 2026-09-01 session's "Follow-up" (relocating into `development/`) is left as accurate
  history, not rewritten** — it correctly describes what was true then. Today's changes supersede
  it via a new entry rather than editing the old one, per this file's own convention.
- **Historical session-summary text is not touched even where it names now-dead paths** (the
  2026-09-01 entry, and the History one-liners below) — only the live "Standing items" pointer was
  updated. This follows the project's existing convention that session summaries are point-in-time
  narrative.

### Next steps

1. None outstanding on this specific task.
2. **Nothing is committed.** If the user wants this committed (including the two deletions'
   staged/unstaged git state), that needs an explicit commit instruction.
3. Recommend running `documentation/biplot-rotation/17_3_minimal-rotation-verification.R` against
   this repo's `method_developments` build (`devtools::load_all()` or `devtools::install()` first)
   to confirm it still passes here — not yet done as part of this housekeeping pass.
4. Carried over, untouched: the `.bl_rotate()` deterministic-rotation proposal (CLAUDE.md §9, now
   pointing at the note's §6) and the optional pkgdown site rebuild (Standing items).

---

## Session summary (2026-09-01) — biplot rotation + optimal pair selection: detailed technical write-up (external deliverable, no code changes)

Goal: the user attached the PhD thesis excerpt `documentation/Optimal rotation of the biplot Phd notes.docx` (§3.5–3.5.2, "Optimal rotation of the biplot") and asked for a detailed technical explanation of (1) the Procrustes-based rotation `.bl_rotate()` performs and (2) how `bl_find_local_cf()` selects which eigenvector pair to rotate to, suitable for reuse in a *different* project. Per the user's explicit request, content was planned before writing.

**No package code changed.** `git status` was clean at session start and remained so throughout except for the one untracked file the user added (the `.docx`); `git diff -- R/` is empty. This was a pure documentation-synthesis task.

### What changed

- **New deliverable, now tracked at `development/biplot_rotation_and_pair_selection.md`** (moved there from the scratchpad — see Follow-up below): full derivation of the Procrustes rotation (the `Y`/`Z` construction, SVD solution `A = V*U'`, `V^new = VA'`, the thesis's worked examples), the code mapped line-by-line to `.bl_rotate()`/`bl_find_local_cf()`, and a "known issues" section (rotation-plane under-determination, within-pair/cross-pair metric mismatch, the explicit global-fidelity trade-off, the CVA class-mean cost). Verified against every source file before delivery (formulas, code snippets, and all numeric claims).
- **`.claude/plans/rippling-herding-salamander.md`** (gitignored) — the approved content plan for the write-up, written and approved via the plan-mode workflow before the document itself was written.
- **`CLAUDE.md`** updated: §2 now cross-references the write-up as required reading before touching `.bl_rotate()` or the pair search; §5 gained a new gotcha entry for the docx/doc-format dead end below.
- `documentation/Optimal rotation of the biplot Phd notes.docx` — read in full this session, then **removed from the repo by the user** once its content was consolidated into the write-up above (see Follow-up).

### Follow-up (same session) — relocated into `development/`, source docx removed

The user removed `documentation/Optimal rotation of the biplot Phd notes.docx` (its content is now
fully consolidated in the write-up) and asked for `documentation/rotation_plane_underdetermination_note.md`
and the scratchpad write-up to both move into `development/`, alongside the main dev/testing
harness. This **supersedes** the original "scratchpad only" decision below.

- `git mv documentation/rotation_plane_underdetermination_note.md development/rotation_plane_underdetermination_note.md` — history preserved, staged (not committed).
- `development/biplot_rotation_and_pair_selection.md` written as a new tracked file (was previously scratchpad-only; not yet `git add`ed).
- Path references repaired everywhere they appeared: `CLAUDE.md` §2 (bullet now points at `development/biplot_rotation_and_pair_selection.md`, docx reference removed), `CLAUDE.md` §2 documentation-folder-convention table (the `development/` row now explicitly allows generic technical notes, not just the harness + its walkthroughs), `CLAUDE.md` §9 (rotation-plane deferred-work entry), and `2 implementation_summary.txt` §4.2 caveat — all three `documentation/rotation_plane_underdetermination_note.md` mentions updated to `development/...`. Checked the auto-memory `reference_docs.md` for the memory-sync rule — it doesn't mirror either file, so no memory update needed.
- **Nothing staged this follow-up is committed** — the `git mv` stages itself (git's normal behaviour), but no `git add`/`git commit` was run.

### Dead ends / issues hit

1. **File lock.** The first read attempt on the `.docx` failed with a `PermissionError` (`withBinaryFile: permission denied`) — the file was open in Word (or mid-OneDrive-sync) at that moment. Resolved once the user closed it.
2. **Wrong file format under a `.docx` extension.** `file <path>` reported `Composite Document File V2 Document` (legacy binary `.doc`), not `Microsoft Word 2007+`. `pandoc -t markdown` failed with `couldn't unpack docx container: Did not find end of central directory signature` — pandoc cannot read legacy binary `.doc` regardless of the file's extension.
3. **Word COM automation as a scripted workaround — tried, failed, abandoned.** `New-Object -ComObject Word.Application` + `$doc.Documents.Open()` succeeded (72 paragraphs read back), but `$doc.SaveAs2(path, wdFormatText)` failed with an opaque `Command failed` COM exception. Suspected but unconfirmed cause: the destination path under the session scratchpad's long GUID-based directory (~220 characters) tripping a legacy path-length limit inside Word's COM layer — not confirmed, because the user resolved the actual root cause directly (re-saved the file as real `.docx` in Word) before the hypothesis needed testing. Once re-saved, plain `pandoc -t markdown --extract-media=<dir>` converted it immediately. **Recorded as a new gotcha in `CLAUDE.md` §5** so a future session doesn't re-attempt Word COM as a repair path.
4. Not a dead end but a process note: the harness switched into Plan Mode mid-task for what was fundamentally a writing/synthesis task rather than a code change. Handled by treating "the document's structure and sources" as the thing being planned, then writing the deliverable after explicit approval.

### Decisions

- ~~**The deliverable stays in the session scratchpad, not committed to `documentation/`.**~~ — **superseded by the Follow-up above.** The user explicitly asked to move it (and the rotation-plane note) into `development/` instead, once the source docx was consolidated and removed.
- **Format: single Markdown file, full derivation + code + caveats depth** — both confirmed with the user via `AskUserQuestion` against a Word-doc-output alternative and a lighter conceptual-overview alternative.
- **CLAUDE.md gained several additions** (see above and Follow-up) rather than a new standalone note, since none of them is itself a methodological change to the package — reading-list cross-references, a documentation-folder-convention update, and a tooling gotcha.
- **`development/` now also holds generic (non-script-specific) technical notes**, not just the harness script and its per-script walkthroughs — a deliberate widening of that folder's role, done at the user's explicit request. `CLAUDE.md` §2's convention table was updated to describe this rather than leaving it undocumented.

**Test status: unchanged, 136 PASS / 0 FAIL / 4 WARN baseline — not re-run this session (no `R/` file touched).**

### Next steps

1. None outstanding on this specific task — the write-up was written, verified against sources, relocated to `development/` alongside the rotation-plane note, and all path references repaired.
2. **Nothing is committed.** `git mv` of the rotation-plane note is staged; the new `development/biplot_rotation_and_pair_selection.md` and the `CLAUDE.md`/`2 implementation_summary.txt` edits are unstaged. If the user wants this committed, that needs an explicit commit instruction.
3. Carried over, untouched this session: the `.bl_rotate()` minimal-rotation option (CLAUDE.md §9) and the optional pkgdown site rebuild (progress.md standing items).

---

## History (one line per session, newest first)

Full detail for any entry: `git log -p progress.md` (or `git show <hash>`).

- **2026-08-03** — split `scripts/` into a gitignored personal playground and a new tracked
  `development/` folder (main dev/testing harness + its three code walkthroughs); repaired path
  references; committed and fast-forward-merged `method_developments` into `main` (Mac-tester gate
  consciously overridden), then pushed both branches. 136 PASS / 0 FAIL / 4 WARN unchanged.

- **2026-07-15** — script 17 (`scripts/17_loan-gam-pca-3d-4var-surface.R`) target-alignment
  rotation experiments surfaced the `.bl_rotate()` plane under-determination finding: bases
  differing by 4.6e-15 give rotation planes 65.9 deg apart, and only 46% of the plane lies in the
  `proj_pair` eigenplane it is named after. Documented in
  `documentation/rotation_plane_underdetermination_note.md` and CLAUDE.md §9; no package code
  changed; folded into commit `c6ca185` alongside the 2026-08-03 restructure below. Full
  derivation, the ruled-out false leads, and the minimal-rotation alternative considered: see
  `git log -p progress.md` for this entry's original full write-up.
- **2026-06-17** — `355e3e3` round the Shapley plot counterfactual label to 2 d.p. (observed
  stays 3 d.p.); CLAUDE.md §9 note on the resulting d.p. inconsistency left uncommitted.
- **2026-06-16** — `0f47602` fix a stale import and non-ASCII chars in `biplot_grid.R`; doc updates.
- **2026-06-16** — raw-unit scaling support (`00b43a5`): `bl_set_scaling()` +
  `.bl_rescale_biplot_axes()` (raw-unit biplot axes), `.scale_to_raw()` (raw-unit feature values in
  Shapley/sparse/target prints; `kind="level"` vs `"delta"` is the correctness rule), and
  `bl_local$bl_counterfactual` mirroring `bl_target`. All display-only. Script 03 fixed to
  standardise on the training split only (was leaking test data). 77 -> 136 PASS.
- **2026-06-16** — `new_title = NA` extended to every biplot `plot()` method (`dc0c3c4`);
  established as a CLAUDE.md §5 convention. Dead end: the plan's predicted "invalid 'main'"
  error was wrong — `title()` coerces almost anything; only a closure/environment errors.
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
