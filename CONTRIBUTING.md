# Contributing to Skills for Molecular Simulation

Thanks for your interest in improving this project. This guide covers the
most common contribution types and the quality bar PRs must meet.

## Project context

This repo publishes **structured derivatives** (cheatsheet / patterns /
glossary / formulas / derivations) from 6 books on molecular simulation,
CCUS, and adsorption, plus the Chinese translation of Frenkel & Smit's
*Understanding Molecular Simulation* (3rd ed.). We deliberately do not
host original book text or figures (copyright).

Before contributing, please read:

- [`README.md`](README.md) — project background and structure
- [`docs/skills/index.md`](docs/skills/index.md) — current skill coverage
- [`docs/skills/molsim-cn/derivations.md`](docs/skills/molsim-cn/derivations.md) —
  the quality benchmark for new content

## How to contribute

### 1. Open an issue first

For any non-trivial change (new book, new derivation, terminology
unification across chapters), open an issue to align on scope before
opening a PR. Trivial fixes (typos, single-formula corrections) can go
straight to PR.

Use one of the issue templates:
- **Bug report** — site build broken, broken link, math rendering issue
- **Content fix** — translation error, wrong formula, missing boundary
- **New book request** — propose adding a new source book

### 2. Local development

```bash
git clone https://github.com/KASSARhahaha/skills-4-Molecular-Simulation.git
cd skills-4-Molecular-Simulation
pip install mkdocs-material
mkdocs serve   # http://localhost:8000
```

Required for content edits:
- `mkdocs build --strict` must pass with **zero warnings** before commit.
  This catches broken links, orphaned files, malformed frontmatter.

### 3. Content quality bar

For translation fixes:
- Match the existing terminology (see `docs/skills/molsim-cn/glossary.md`).
  If you propose a terminology change, it must be applied consistently
  across all affected chapters.
- Preserve LaTeX math: inline `$...$`, display `$$...$$`. ASCII
  pseudo-formulas like `acc = min(1, exp(-bDU))` are not acceptable.

For formula / derivation additions:
- Copy equations verbatim from the source book when possible. If
  re-deriving, show the algebraic steps.
- Every formula must have: meaning, boundary conditions, units
  (when applicable). See `formulas.md` for the format.

For new books:
- Must be distributable (your own work, public domain, or with permission).
- Use [book-to-skill](https://github.com/virgiliojr94/book-to-skill) to
  generate the initial scaffold, then refine by hand.
- Provide at minimum: `SKILL.md`, `cheatsheet.md`, `patterns.md`,
  `glossary.md`. `formulas.md` / `derivations.md` optional but encouraged
  for technical books.

### 4. Commit conventions

Use [Conventional Commits](https://www.conventionalcommits.org/):

```
feat(skills/<slug>): <what you added>
fix(chapters/<slug>): <what you fixed>
docs: <documentation change>
chore: <build, ci, tooling>
```

Keep commits small and focused. One logical change per commit.

### 5. Pull request

- Reference the issue (e.g. `Closes #42`).
- Include before/after snippets for content changes.
- Confirm `mkdocs build --strict` passes locally.
- For formula-heavy changes, paste the rendered formula (MathJax output)
  in the PR description so reviewers can verify without pulling.

## Code of conduct

Be respectful. Disagreements about terminology or pedagogy are normal —
argue from evidence (citations, source passages) rather than authority.

## License

By contributing, you agree your contributions are licensed under the
project's [MIT license](LICENSE).

## Maintainer

**Kai Yang** (@KASSARhahaha) — sole maintainer. Review turnaround is
typically 1–2 weeks. If a PR sits longer, ping via issue.
