.PHONY: help install install-dev lint lint-check format test test-cov build clean clean-build clean-pyc pre-commit deploy

help:
	@echo "MONSDA Development Commands:"
	@echo ""
	@echo "Setup:"
	@echo "  make install          Install package in development mode"
	@echo "  make install-dev      Install with development dependencies"
	@echo ""
	@echo "Code Quality:"
	@echo "  make lint             Run isort and flake8 linters"
	@echo "  make lint-check       Check linting without fixing (CI mode)"
	@echo "  make format           Format code with isort"
	@echo "  make pre-commit       Run pre-commit hooks"
	@echo ""
	@echo "Testing:"
	@echo "  make test             Run unit tests"
	@echo "  make test-cov         Run tests with coverage report"
	@echo "  make test-integration Run integration tests (requires conda)"
	@echo ""
	@echo "Building:"
	@echo "  make build            Build distribution packages"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean            Remove all build/test artifacts"
	@echo "  make clean-build      Remove build artifacts"
	@echo "  make clean-pyc        Remove Python cache files"
	@echo ""

install:
	python -m pip install --upgrade pip
	python -m pip install -e .

install-dev: install
	python -m pip install pytest pytest-cov flake8 isort pre-commit build wheel twine

lint:
	isort MONSDA/ tests/
	flake8 MONSDA/ tests/

lint-check:
	isort --check-only --diff MONSDA/ tests/
	flake8 MONSDA/ tests/ --count --statistics --show-source

format:
	isort MONSDA/ tests/

test:
	pytest tests/test_Utils.py -v

test-cov:
	pytest tests/test_Utils.py -v --cov=MONSDA --cov-report=term-missing --cov-report=html

test-integration:
	cd tests && RUN_INTEGRATION_TESTS=1 bash cicd_test.sh

build: clean
	python -m pip install --upgrade pip
	python -m pip install build wheel setuptools versioneer
	python -m build
	python -m pip install twine
	twine check dist/*

clean: clean-build clean-pyc

clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .eggs/
	rm -rf *.egg
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/

clean-pyc:
	find . -type f -name '*.py[cod]' -delete
	find . -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name '*.egg-info' -exec rm -rf {} + 2>/dev/null || true

pre-commit:
	pre-commit install
	pre-commit run --all-files

deploy: build
	python -m twine upload dist/*
