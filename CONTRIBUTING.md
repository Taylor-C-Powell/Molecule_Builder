# Contributing to MolBuilder

Thanks for your interest in contributing to MolBuilder! This document covers the development setup and guidelines for submitting changes.

## Development Setup

```bash
git clone https://github.com/Taylor-C-Powell/Molecule_Builder.git
cd Molecule_Builder
pip install -e ".[dev]"
```

This installs MolBuilder in editable mode along with development dependencies (pytest, pytest-cov, ruff).

## Running Tests

```bash
# Full test suite
python -m pytest tests/ -v

# Quick run
python -m pytest tests/ -q

# With coverage
python -m pytest tests/ --cov=molbuilder --cov-report=term-missing
```

All tests must pass before submitting a pull request.

## Code Style

This project uses [Ruff](https://docs.astral.sh/ruff/) for linting:

```bash
ruff check molbuilder/
ruff format molbuilder/
```

Key conventions:
- Line length: 100 characters
- Target: Python 3.11+
- Use type hints for public API functions
- Docstrings for all public classes and functions

## Project Structure

```
molbuilder/
  core/         Constants, elements, geometry, bond data
  atomic/       Bohr model, quantum numbers, wavefunctions
  bonding/      Lewis structures, VSEPR, covalent bond analysis
  molecule/     Molecular graph, conformations, stereochemistry, builders
  smiles/       SMILES tokenizer, parser, writer
  io/           File I/O (XYZ, MOL/SDF, PDB, JSON, SMILES)
  reactions/    Reaction templates, retrosynthesis engine
  process/      Reactor selection, costing, safety, scale-up
  dynamics/     Molecular dynamics, mechanism visualization
  reports/      ASCII report generators
  visualization/ 3D rendering, quantum orbital plots
  gui/          Tkinter-based molecular editor
  cli/          Interactive menus and demos
```

## Adding Reaction Templates

Reaction templates live in `molbuilder/reactions/templates.py`. Each template needs:
- Functional group pattern (reactant and product)
- Reagents, solvents, catalysts
- Temperature and yield ranges
- Retrosynthetic transform logic
- Safety notes

See existing templates for the expected format.

## Scientific Standards

MolBuilder values scientific accuracy:
- Cite sources for physical constants and empirical data
- Validate computed values against known experimental results
- Include tests that verify scientific correctness, not just code behavior
- Document known limitations (e.g., Bohr model accuracy for multi-electron systems)

## Submitting Changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make your changes with tests
4. Run the full test suite: `python -m pytest tests/ -v`
5. Run the linter: `ruff check molbuilder/`
6. Submit a pull request with a clear description of the change

## Reporting Issues

Use GitHub Issues for bug reports and feature requests. For bugs, include:
- Python version and OS
- Minimal code to reproduce the issue
- Expected vs. actual behavior

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
