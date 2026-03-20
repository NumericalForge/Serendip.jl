# Styling

This file records lightweight style preferences for Serendip.

## Naming

- Prefer clear, descriptive names over short or cryptic names.
- Prefer `_` as the separator between words in function and variable names instead of joining words together.
- Prefer avoiding a trailing `!` in mutating function names. The fact that a function mutates its arguments should be made clear in its documentation.
- Prefer keyword names that describe behavior directly, especially for public APIs.

## Public APIs

- Keep public function signatures compact and focused on the options that materially affect behavior.
- Prefer removing stale or low-value options instead of carrying compatibility-only keywords indefinitely.
- Keep docstrings aligned with the actual behavior and defaults of the implementation.
- Whether a function mutates state should be stated clearly in its documentation and inferred from that documentation.

## Changes

- When simplifying or refactoring an API, update exports, tests, and documentation together.
- Prefer small helper functions when they clarify behavior, but avoid keeping dead internal machinery after a public API is reduced.
