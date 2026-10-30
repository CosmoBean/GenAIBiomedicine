# Submission Polish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the homework submission easier to read and run without altering analysis results.

**Architecture:** Keep the existing ordered analysis scripts as the reproducible implementation. Update only report callout metadata and the README, then verify their consistency with a portable shell test.

**Tech Stack:** LaTeX, Markdown, POSIX shell.

**Spec:** `docs/superpowers/specs/2026-10-29-submission-polish-design.md`

## Global Constraints

- Do not change Python analysis code, parameters, input data, or generated scientific outputs.
- Preserve `bash submission/run_all.sh` as the single reproducibility command.
- Use the title `Independent student check` for all completed verification callouts.

## Review Focus

- Completed-check language must not say that a student still needs to complete work.
- The README command must be runnable from the repository root.
- The README must identify the ordered scripts and output locations.

### Task 1: Validate and polish submission metadata

**Files:**
- Create: `submission/tests/test_submission_metadata.sh`
- Modify: `submission/report/report.tex:23,31`
- Modify: `submission/README.md`

**Interfaces:**
- Consumes: report callout macro and README entry-point instructions.
- Produces: consistent, grader-facing submission metadata.

- [ ] **Step 1: Write the failing metadata test**

Assert that the report uses `Independent student check`, contains no stale
unfinished-check wording, and that the README includes `bash submission/run_all.sh`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `bash submission/tests/test_submission_metadata.sh`
Expected: FAIL because the report still says `Student verification` and `needs student verification`.

- [ ] **Step 3: Update the report labels and simplify the README**

Rename the callout title and evidence index. Replace the README with a concise
run guide that retains reproducibility details required by the assignment.

- [ ] **Step 4: Run the metadata test and build the report**

Run: `bash submission/tests/test_submission_metadata.sh && bash submission/report/build_report.sh`
Expected: PASS and a successfully compiled PDF.

- [ ] **Step 5: Commit**

Commit the documentation, metadata test, and regenerated report.
