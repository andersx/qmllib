.PHONY: install install-dev test check format typing clean help

install:
	pip install -e .[test] --verbose

install-dev:
	pip install -e .[test,dev] --verbose
	pre-commit install

test:
	pytest

check: format typing

format:
	ruff format python/ tests/
	ruff check --fix python/ tests/

typing:
	mypy python/ tests/

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true

environment:
	conda env create -f environments/environment-dev.yaml
