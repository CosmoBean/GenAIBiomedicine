# Submission polish design

## Intent

Make the completed homework easy for a grader to navigate without changing any
scientific analysis, reported result, or generated output.

## Goals

- Rename completed verification callouts to **Independent student check**.
- Remove instructions implying that those checks are still unfinished.
- Make the README a short, direct submission guide: one command, ordered
  analyses, sources, and folder map.

## Non-goals

- Do not change the Python analysis logic, parameters, data, or outputs.
- Do not alter the report's scientific claims.

## Approach

Retain the existing checkpoint-specific scripts because their names and order
provide the clearest reproducibility trail. Simplify only the entry-point
documentation and stale completion wording.

## Validation

A lightweight shell test will assert the new callout title, absence of stale
unfinished-check language, and the README's reproducible run command.
