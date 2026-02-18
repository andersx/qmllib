.PHONY: install install-dev test check format typing clean help

install:
	pip install -e .[test] --verbose

install-dev:
	pip install -e .[test,dev] --verbose
	pre-commit install

test:
	pytest

check: format typing
	@echo "✅ All code quality checks passed!"

format:
	@echo "Running ruff format..."
	ruff format python/ tests/
	@echo "Running ruff lint with auto-fix..."
	ruff check --fix python/ tests/
	@echo "✅ Code formatting complete!"

typing:
	@echo "Running mypy type checking..."
	mypy python/ tests/
	@echo "✅ Type checking complete!"

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

help:
	@echo "Available targets:"
	@echo "  install      - Install package in editable mode with test dependencies"
	@echo "  install-dev  - Install with dev dependencies and setup pre-commit hooks"
	@echo "  test         - Run pytest test suite"
	@echo "  check        - Run all code quality checks (format + typing)"
	@echo "  format       - Format code with ruff"
	@echo "  typing       - Run mypy type checking"
	@echo "  clean        - Remove Python cache files and build artifacts"
	@echo "  environment  - Create conda environment"
	@echo "  help         - Show this help message"
