# The final thesis in the repository; citations in line; session close

Date created: 2026-11-06
Date last edited: 2026-11-06

> **IMPLEMENTED 2026-11-06.** See "Outcome" at the end.

(The previous plans in this file are implemented and archived in `.claude/reference/`:
`poly-clip-single-definition-plan.md`, `filter-switches-in-set-filters-plan.md`,
`local-plot-hull-trim-plan.md`. On approval, this plan is archived as
`.claude/reference/final-thesis-citations-plan.md`.)

## Context

The user replaced the thesis PDF in `documentation/`: the 30 September 2025 draft
(`Adriaan Rowan - 13551361 - Thesis - 30 Sept 2025.pdf`, tracked) is deleted, and the final version
`documentation/FINAL_ADRIAAN ROWAN_27022026.pdf` (untracked) is added. They also deleted
`documentation/XAI PhD presentation - March 2026.pdf` (tracked; nothing cites it). They ask that every
reference to the September 2025 version be removed, that code and documentation agree with the final
thesis, and that CLAUDE.md and progress.md close the session.

**The new file is the one already cited.** Its MD5 (`14e1353989d9ea4fc3096387cb49c6bd`) equals that of
`../Explain XAI code/reference-materials/1 FINAL_ADRIAAN ROWAN_27022026.pdf`, which the notes have cited
since 2026-11-05 (234 PDF pages; printed page = PDF page − 16; Appendix A from p.158, A7 p.207).

**Decisions (user, 2026-11-06):** history records keep their mentions of the draft (the archived plans
in `.claude/reference/` and the decision form `future-improvements-review.md`); the documents name the
thesis "the final version of 27 February 2026".

**The audit** (two read-only passes, every thesis citation in `R/`, NEWS, README, vignette, tests,
the rotation scripts, CLAUDE.md, the root documents, the walkthroughs and the notes; findings checked
against the PDF): 27 citations in code and about 60 in the documents match. The mismatches are below,
all in documents and one script label; none in `R/`.

## Changes

### 1. The draft and the thesis's location (living documents)

- `CLAUDE.md:177`: delete "; p.106 in the 30 Sept 2025 draft".
- `documentation/biplot-rotation/biplot_rotation_and_pair_selection.md:32-35`: the thesis is the final
  version of 27 February 2026, `../FINAL_ADRIAAN ROWAN_27022026.pdf` (in `documentation/`); delete the
  sentence on the draft. `:1471`: the same file name and location.
- `documentation/boundary-distance-and-robustness.md:375-376`: the same.

### 2. Wrong section numbers (the draft's §3.9.3 is the final's §3.10, "Using biplots to explain machine learning models", from p.92)

- `CLAUDE.md:23`: the rotation note's §2 covers "thesis §3.5–3.5.2, §3.8.3–3.10, §4.3.1" (was "§3.4,
  §3.7.3–3.9, §4.2.1", which are the CVA biplot, grid search and the global biplot in the final).
- Rotation note: `:105` "§3.8.3–3.9.3 (pp.82–93)" → "§3.8.3–3.10"; `:204` "(§3.9.3, p.92)" → "(§3.10,
  p.92)"; `:234` heading "thesis §3.9.3 p.93" → "thesis §3.10 p.93"; `:284` "p.92 (§3.9.3)" → "p.92
  (§3.10)"; `:1473` "§3.8.3, §3.9.1 and §3.9.3 (pp.82–93)" → "§3.8.3, §3.9.1, §3.9.3 and §3.10". A §3.9.3
  citing p.91 content stays. Then grep the repository for any other "§3.9.3" citing pp.92–93.
- `documentation/biplot-rotation/17_3_minimal-rotation-verification.R:525`: the printed label "thesis
  (§2.3)" → "thesis (§3.5, p.55)" (the 28° worked example; §2.3 is the note's section).

### 3. Quotes made exact

- Rotation note `:203`: "either the CVA and PCA biplots" → "either the CVA or the PCA biplots"; `:259`:
  "the total within-group SSP" (p.118); `:277`: drop "either" from the p.86 quote; `:284`: "since it
  ensures" → "as it ensures".
- `boundary-distance-and-robustness.md:185`: divided by "the overall standard deviation of that
  variable" (p.103).
- `local-convex-hull.md:121`: the thesis writes "unconstraint"; give it as a paraphrase without
  quotation marks ("drops the hull only when an unconstrained solution is acceptable (p.115)").

### 4. Where the thesis defines the three rules

The thesis defines actionable, feasible and plausible in §3.9.2, "Constraints for a biplot-generated
counterfactual value" (pp.88–90), and applies them on p.116. Add "§3.9.2, pp.88–90" next to the p.116
citation in `set_filters()`'s roxygen section (`R/local_cf.R`), `local-convex-hull.md` §6 and the
walkthrough's "The three rules" paragraph (`review_section9_to_15.md`).

### 5. Stale links to the rotation note (not thesis citations, but out of line)

The rotation note was renumbered on 2026-10-31; its old §5.3/§5.9 records of the tie-break drift are now
its Appendix A.2.
- `CLAUDE.md:178` ("recorded in the rotation note §5.3 and §5.9"): point to Appendix A.2.
- `documentation/nearest-vertex-tie-break.md:72` ("the §5.9 sequence"): "the call-order sequence the
  rotation note then recorded"; §8 "Records corrected" keeps its list, with one line that its section
  numbers are those of 2026-10-30 and the material is now Appendix A.2.

## Session close

- **CLAUDE.md** (through `/revise-claude-md`, one line each, narrowest place):
  - §5 "Checking a thesis citation": the PDF is `documentation/FINAL_ADRIAAN ROWAN_27022026.pdf`
    (234 pages; Appendix A from p.158; §3.10 from p.92);
  - §5 "Running R": `Rscript.exe -e` with `|` or `||` in the expression fails from Git Bash; use a `.R`
    file;
  - §5 Pima: Figure 7(a) is drawn from run 4 searched without the hull (`03_local.R`, `run4_no_hull`);
  - §6: the local biplot colours its grid only inside the plane's hull (`trim_to_hull`,
    `.grid_in_hull()`; `hull_polygon` stored also with the hull off).
- **progress.md**: session 11 finished: the thesis swap and audit (what was corrected), test and check
  status, files changed; the "Not mine" note replaced (the user made both deletions); the commit must
  `git add` the new PDF and stage both deletions; next step unchanged (the Paper 1 plan).
- **Memory** `project_boundarylogic.md` and its `MEMORY.md` line: status as of 2026-11-06 (sessions 10–11
  uncommitted; next step).

## Not changed

The archived plans and `future-improvements-review.md` (history, by decision); the methodology note
(its spot-check found no mismatch, so no re-render); `docs/` (pkgdown, rebuilt before the release).

## Verification

1. grep the living documents for `30 Sept`, `September 2025`, `13551361`, `reference-materials`,
   `outside the repository` (thesis sense): none left; grep `§3.9.3`: only p.91 citations.
2. `devtools::document()`; `devtools::test()`: 527 PASS, 0 FAIL, 4 WARN; `R CMD check`: 0 / 0 / 0.
3. `17_3` run against the last log: only the label line and timings differ.
4. Archive this plan; stop and report. No commit or push without an instruction.

## Outcome (2026-11-06)

Implemented as planned.
- grep of the living documents: no draft or old-location mention left; the two remaining §3.9.3
  citations point to p.91 (correct).
- `devtools::test()`: 527 PASS / 0 FAIL / 4 WARN; `R CMD check`: 0 / 0 / 0 (`set_filters.Rd`
  regenerated).
- `17_3`: only the label line ("thesis (§3.5, p.55)") and timings differ.
- CLAUDE.md (four learnings via `/revise-claude-md`), progress.md and the memory updated.
- Not committed.
