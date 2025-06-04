# =============================================================================
# Makefile for multicolorfits documentation testing
# =============================================================================

.PHONY: docs docs-clean docs-test docs-github-test docs-serve docs-linkcheck

# Variables
PYTHON := python3.11
DOCS_DIR := ./
BUILD_DIR := $(DOCS_DIR)/_build
SOURCE_DIR := $(DOCS_DIR)/source
VENV_DIR := .docs_test_env

# =============================================================================
# STANDARD DOCUMENTATION COMMANDS
# =============================================================================

docs: ## Build documentation normally
	@echo "📚 Building documentation..."
	cd $(DOCS_DIR) && sphinx-build -b html $(SOURCE_DIR) $(BUILD_DIR)/html

docs-clean: ## Clean documentation build
	@echo "🧹 Cleaning documentation build..."
	rm -rf $(BUILD_DIR)

docs-serve: docs ## Build and serve documentation locally
	@echo "🌐 Serving documentation at http://localhost:8000"
	cd $(BUILD_DIR)/html && $(PYTHON) -m http.server 8000

# =============================================================================
# GITHUB-LIKE TESTING
# =============================================================================

docs-github-test: ## Test documentation build exactly like GitHub Actions
docs-linkcheck: ## Check for broken links
	@echo "🔗 Checking for broken links..."
	cd $(DOCS_DIR) && sphinx-build -b linkcheck $(SOURCE_DIR) $(BUILD_DIR)/linkcheck

docs-test-all: docs-clean docs-github-test docs-linkcheck ## Run all documentation tests
	@echo "🎉 All documentation tests passed!"

# =============================================================================
# DEBUGGING COMMANDS
# =============================================================================

docs-debug: ## Build documentation with maximum verbosity for debugging
	@echo "🐛 Building documentation with debug output..."
	cd $(DOCS_DIR) && sphinx-build -b html $(SOURCE_DIR) $(BUILD_DIR)/html -v -T -E -a

docs-env-info: ## Show environment information for debugging
	@echo "🔍 Environment Information:"
	@echo "Python version: $$($(PYTHON) --version)"
	@echo "Sphinx version: $$($(PYTHON) -c 'import sphinx; print(sphinx.__version__)' 2>/dev/null || echo 'Not installed')"
	@echo "Working directory: $$(pwd)"
	@echo "Python path: $$($(PYTHON) -c 'import sys; print(sys.path)')"
	@if [ -f "pyproject.toml" ]; then echo "Using pyproject.toml"; fi
	@if [ -f "$(SOURCE_DIR)/requirements.txt" ]; then echo "Using requirements.txt"; fi

# =============================================================================
# GITHUB ACTIONS SIMULATION
# =============================================================================

docs-act-test: ## Test using act (GitHub Actions locally) - requires act to be installed
	@echo "🎭 Testing with act (GitHub Actions runner locally)..."
	@if command -v act > /dev/null; then \
		act -j build-docs; \
	else \
		echo "❌ act not installed. Install with:"; \
		echo "   macOS: brew install act"; \
		echo "   Linux: curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash"; \
	fi

# =============================================================================
# DOCKER TESTING
# =============================================================================

docs-docker-test: ## Test documentation build in Docker (exact GitHub environment)
	@echo "🐳 Testing documentation with Docker..."
	@docker run --rm -v "$$(pwd):/workspace" -w /workspace ubuntu:22.04 bash -c '\
		apt-get update -qq && apt-get install -y -qq $(PYTHON) $(PYTHON)-venv python3-pip git && \
		$(PYTHON) -m venv /tmp/docs_env && \
		source /tmp/docs_env/bin/activate && \
		python -m pip install --upgrade pip -q && \
		if [ -f "pyproject.toml" ]; then \
			pip install -e .[docs] -q; \
		else \
			pip install -r $(SOURCE_DIR)/requirements.txt -q; \
		fi && \
		cd $(DOCS_DIR) && \
		sphinx-build -b html $(SOURCE_DIR) $(BUILD_DIR)/html -W --keep-going'
	@echo "✅ Docker test complete!"

# =============================================================================
# CONTINUOUS INTEGRATION HELPERS
# =============================================================================

docs-ci-install: ## Install documentation dependencies (for CI)
	@$(PYTHON) -m pip install --upgrade pip
	@if [ -f "pyproject.toml" ]; then \
		pip install -e .[docs]; \
	else \
		pip install -r $(SOURCE_DIR)/requirements.txt; \
	fi

docs-ci-build: ## Build documentation (for CI)
	@cd $(DOCS_DIR) && sphinx-build -b html $(SOURCE_DIR) $(BUILD_DIR)/html -W --keep-going

# =============================================================================
# HELP
# =============================================================================

help: ## Show this help message
	@echo "📖 multicolorfits Documentation Commands"
	@echo ""
	@echo "🚀 Quick Start:"
	@echo "  make docs-github-test  # Test like GitHub would"
	@echo "  make docs-serve        # Build and serve locally"
	@echo ""
	@echo "📋 Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-20s %s\n", $$1, $$2}'
	@echo ""
	@echo "💡 Pro tip: Run 'make docs-github-test' before every commit!"

# Default target
.DEFAULT_GOAL := help

docs-github-test:
	@echo "🧪 Testing documentation build like GitHub CI..."
	@echo "📦 Creating fresh virtual environment..."
	@rm -rf $(DOCS_DIR)/.docs_test_env
	@python3 -m venv $(DOCS_DIR)/.docs_test_env
	@echo "⬆️  Upgrading pip..."
	@$(DOCS_DIR)/.docs_test_env/bin/pip install --upgrade pip > /dev/null 2>&1
	@echo "📚 Installing documentation dependencies..."
	@$(DOCS_DIR)/.docs_test_env/bin/pip install -r $(SOURCE_DIR)/requirements.txt > /dev/null 2>&1
	@echo "🏗️  Building documentation with strict error checking..."
	@cd $(DOCS_DIR) && ./.docs_test_env/bin/sphinx-build -W -b html source _build/html; \
	BUILD_EXIT_CODE=$$?; \
	if [ "$$BUILD_EXIT_CODE" -eq 0 ]; then \
		echo "✅ Documentation build succeeded (GitHub CI will pass)"; \
		rm -rf .docs_test_env; \
		exit 0; \
	else \
		echo "❌ Documentation build failed (same as it would on GitHub)"; \
		rm -rf .docs_test_env; \
		exit 1; \
	fi

docs-github-test-simple:
	@echo "🧪 Testing documentation build..."
	@rm -rf $(DOCS_DIR)/.docs_test_env
	@python3 -m venv $(DOCS_DIR)/.docs_test_env
	@$(DOCS_DIR)/.docs_test_env/bin/pip install --upgrade pip > /dev/null
	@$(DOCS_DIR)/.docs_test_env/bin/pip install -r $(SOURCE_DIR)/requirements.txt > /dev/null
	@cd $(DOCS_DIR) && ./.docs_test_env/bin/sphinx-build -W -b html source _build/html
	@echo "✅ Build completed - check output above for any warnings/errors"
	@rm -rf $(DOCS_DIR)/.docs_test_env
