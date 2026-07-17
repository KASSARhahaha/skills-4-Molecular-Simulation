# Changelog

All notable changes to this project are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Extend deep upgrade (formulas + derivations) to remaining 5 books
- Merged full-book `.docx` output (cross-chapter `\ref` resolution)
- Bilingual search (zh↔en term alignment, cross-book citations)
- EPFL/UC Berkeley Day 0–9 exercise integration

## [0.3.0] - 2026-07-17

### Added
- `LICENSE` file (MIT) with content-licensing note
- `CONTRIBUTING.md` with PR / issue workflow
- `CHANGELOG.md` (this file)
- Issue templates (bug report, content fix, new book request)
- CodeQL security scanning workflow
- 11 repository topics for discoverability
- README badges (license / deploy / site / lang / mkdocs)
- README roadmap, contribution guidelines, acknowledgments

### Changed
- README strengthened from 18 to 161 lines (6-book table, three-tier
  knowledge structure, project tree)

## [0.2.0] - 2026-07-06

### Added
- `skills/molsim-cn/derivations.md` — first-principles derivations of 10
  core algorithms (Metropolis, Verlet, Nosé-Hoover, LJ tail correction,
  Ewald, FEP, TI, BAR, Crooks/Jarzynski, Green-Kubo); 612 lines,
  94 display formulas
- `skills/molsim-cn/formulas.md` — hand-curated 25+ core equations with
  physical meaning, boundary conditions, and units
- "公式" and "推导" columns in `docs/skills/index.md`

### Changed
- `strip_latex.py` rewritten to preserve math verbatim (sentinel-based
  protection); recovered 5525 inline + 1361 display formulas (was 0)
- All 26 molsim-cn chapter summaries regenerated with proper LaTeX math
  (ASCII pseudo-formulas → real `$...$` / `$$...$$`)
- `skills/molsim-cn/cheatsheet.md` and `patterns.md` upgraded to LaTeX

## [0.1.0] - 2026-07-05

### Added
- Initial public release
- 6-book skill derivatives published: `molsim-cn`, `molsim-tutorial-exercises`,
  `ccs-intro`, `adsorption-sim`, `adsorption-sim-notes`, `aspen-adsorption`
- 26 Chinese chapter translations (`docs/chapters/`) from Frenkel & Smit
  3rd edition (482-page XeLaTeX build)
- `docs/skills/index.md` navigation page
- `.github/workflows/ci.yml` — GitHub Actions auto-deploy to Pages
- `mkdocs.yml` with MathJax / arithmatex configuration
- 53 broken `../../images/` → `../images/` image path fixes in chapters

### Notes
- Repo created via mirror-push from `molsim-online` (archived).
- Public URL: https://kassarhahaha.github.io/skills-4-Molecular-Simulation/
