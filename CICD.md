# CI/CD Pipeline for MONSDA

This document describes the Continuous Integration/Continuous Deployment (CI/CD) pipeline for MONSDA.

## Overview

The CI/CD pipeline is implemented using GitHub Actions and includes the following stages:

1. **Linting** - Code quality checks with flake8 and isort
2. **Testing** - Unit tests with pytest and coverage reporting
3. **Building** - Package building and artifact verification
4. **Integration Testing** - Full pipeline integration tests (optional)
5. **Publishing** - Automatic deployment to PyPI on releases

## Pipeline Stages

### 1. Linting (Runs on all PRs and pushes)

**Job:** `lint`

Automatically checks code quality using:
- **isort** - Validates import statement ordering
- **flake8** - Checks PEP 8 compliance and code style

**Configuration:**
- Line length: 110 characters
- Ignores: E203, E266, E501, W503
- Continues on error (warnings don't block)

Run locally:
```bash
make lint-check      # Check without fixing
make lint            # Check and fix issues
make format          # Format imports only
```

### 2. Unit Testing (Runs after linting)

**Job:** `tests`

Runs pytest test suite with coverage reporting:
- Tests located in `tests/test_Utils.py`
- Generates coverage reports (XML and HTML)
- Uploads coverage to Codecov

**Dependencies:**
- Conda environment from `environment.yml`
- Python 3.12.2
- pytest and pytest-cov

Run locally:
```bash
make test            # Run unit tests
make test-cov        # Run with coverage report
```

### 3. Building (Runs after tests)

**Job:** `build`

Builds distribution packages:
- Creates wheel and source distributions
- Validates build artifacts with twine
- Stores artifacts for 7 days

**Tools:**
- setuptools
- wheel
- versioneer (for version management)
- twine (for package validation)

Run locally:
```bash
make build           # Build distribution packages
make clean-build     # Remove build artifacts
```

### 4. Full Pipeline Integration Test (Optional)

**Job:** `full-pipeline`

Runs comprehensive integration tests of the Snakemake/Nextflow pipeline functionality.

**Trigger conditions:**
- Manual workflow dispatch with `run_full_pipeline=true`
- PR with `full-pipeline` label

**Process:**
1. Sets up Conda environment from `environment.yml`
2. Installs MONSDA package
3. Updates test configuration with current version
4. Runs test workflow with sample data
5. Uploads debug artifacts on failure

Run locally:
```bash
RUN_INTEGRATION_TESTS=1 bash tests/cicd_test.sh
```

### 5. Publishing to PyPI (Runs on version tags)

**Job:** `publish`

Automatically publishes releases to PyPI when tags are pushed.

**Trigger:** Tags matching pattern `v*` (e.g., `v1.0.0`)

**Requirements:**
- Build job must succeed
- PyPI environment configured with trusted publishing
- Token in GitHub secrets (automated via Trusted Publisher)

**Manual deployment:**
```bash
make deploy          # Build and upload to PyPI
```

## Workflow Triggers

The pipeline runs automatically on:

- **Pull requests** - All PR types (opened, synchronize, reopened, labeled, unlabeled)
- **Push to main/develop** - On code pushes to main and develop branches
- **Version tags** - On tags matching `v*` (triggers PyPI publishing)
- **Manual dispatch** - Triggered via GitHub Actions UI with options

## Local Development Setup

### Prerequisites

- Python 3.12+
- Conda/Mamba (for integration tests)
- Git

### Quick Setup

```bash
# Install development dependencies
make install-dev

# Set up pre-commit hooks
make pre-commit
```

### Available Commands

```bash
# Code quality
make lint            # Run isort and flake8 (auto-fix)
make lint-check      # Check linting (don't fix)
make format          # Format imports with isort

# Testing
make test            # Run unit tests
make test-cov        # Run tests with coverage
make test-integration # Run full pipeline tests

# Building
make build           # Build packages
make clean           # Clean all artifacts
make clean-pyc       # Remove Python cache files

# Versioning
make sync-config-versions   # Stamp config VERSION fields with the git-tag version
make check-config-versions  # Verify config VERSION fields are in sync

# Useful
make help            # Show all commands
```

## Pre-commit Hooks

Pre-commit hooks are configured to automatically check code before commits.

### Setup

```bash
make pre-commit
```

### Hooks Included

- Trailing whitespace removal
- End-of-file fixer
- YAML validation
- Large file detection (>1MB)
- Merge conflict detection
- JSON validation
- Line ending normalization
- Import sorting (isort)
- Code linting (flake8)
- Type checking (mypy, optional)

Run manually:
```bash
pre-commit run --all-files
```

## Configuration Files

### GitHub Actions Workflow
- **File:** `.github/workflows/python-app.yaml`
- Defines all CI/CD pipeline stages
- Configures job dependencies and triggers

### Pre-commit Config
- **File:** `.pre-commit-config.yaml`
- Defines hooks run before commits
- Includes linting and formatting tools

### Makefile
- **File:** `Makefile`
- Provides convenient commands for local development
- Mirrors CI/CD pipeline steps

### Package Configuration
- **File:** `pyproject.toml` - Modern Python project metadata
- **File:** `setup.py` - Package installation script
- **File:** `setup.cfg` - Tool configuration
  - isort settings
  - flake8 settings
  - pytest configuration
  - coverage settings

### Environment
- **File:** `environment.yml` - Conda environment for testing
- Includes all runtime and development dependencies

## Troubleshooting

### Linting Failures

If linting fails locally or in CI:

```bash
# Fix imports
make format

# Fix other issues manually (requires code changes)
make lint-check  # See what needs fixing
```

### Test Failures

```bash
# Run tests with verbose output
make test-cov

# Check coverage gaps
# Open htmlcov/index.html in browser
```

### Build Issues

```bash
# Clean previous builds
make clean

# Rebuild
make build

# Check for issues
twine check dist/*
```

### Integration Test Failures

```bash
# Debug locally
cd tests
RUN_INTEGRATION_TESTS=1 bash cicd_test.sh

# Check logs in tests/LOGS directory
```

## Version Management

MONSDA uses **versioneer** for automatic version management:

1. Version is defined in git tags
2. Versioneer automatically derives version from git state
3. Package version matches git tag (e.g., tag `v1.0.0` → package version `1.0.0`)

**Release process:**
1. Create release tag: `git tag v1.0.0`
2. Push tag: `git push --tags`
3. GitHub Actions automatically builds and publishes to PyPI

### Config Version Synchronization

The `VERSION` field in the config files under `configs/` and `tests/data/`
is kept in sync with the versioneer-derived version automatically, so that
`monsda`'s runtime version check (`monsda_check_version`) always passes.

This is handled by `scripts/update_config_versions.py`, which mirrors the same
version source used at runtime (`MONSDA.__version__`, falling back to
`python setup.py --version`). The script is comment-safe (config files contain
`#` comments and are not strict JSON) and only rewrites empty, `FIXME`, or
version-like `VERSION` values, leaving documentation strings untouched.

It runs automatically in CI:
- In the **build** job, before `python -m build`, so the published package ships
  templates matching the release version.
- In the **full-pipeline** job, before the integration test, replacing the old
  `sed`-based patching.

Both jobs use `fetch-depth: 0` on checkout so versioneer can resolve git tags.

Run locally:
```bash
make sync-config-versions    # Stamp config VERSION fields with the git-tag version
make check-config-versions   # Verify config VERSION fields are in sync (CI mode)

# Or invoke the script directly with an explicit version
python scripts/update_config_versions.py --version 1.5.0 configs/*.json
```

## Best Practices

1. **Always run linting before committing:**
   ```bash
   make pre-commit  # or rely on pre-commit hooks
   ```

2. **Run tests locally before pushing:**
   ```bash
   make test-cov
   ```

3. **Use meaningful commit messages** that reference issue numbers

4. **Add tests for new features** to maintain coverage

5. **Keep dependencies updated** by regularly running `pip install --upgrade`

6. **Review CI logs** if builds fail to understand what needs fixing

## Contact & Support

For issues or questions about the CI/CD pipeline:
- Check GitHub Actions logs for specific failures
- Review this documentation
- Open an issue on GitHub
- Contact the maintainers

## References

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [pytest Documentation](https://docs.pytest.org/)
- [flake8 Documentation](https://flake8.pycqa.org/)
- [isort Documentation](https://pycqa.github.io/isort/)
- [Pre-commit Framework](https://pre-commit.com/)
- [Setuptools Documentation](https://setuptools.pypa.io/)
- [Versioneer](https://github.com/python-versioneer/python-versioneer)
