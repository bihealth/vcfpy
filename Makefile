.PHONY: default
default: help

.PHONY: help
help:
	@echo "Available targets:"
	@echo "  help                   Show this help message"
	@echo "  fix                    Format source code"
	@echo "  check                  Run checks"
	@echo "  test                   Format source code"
	@echo "  docs                   Build documentation"
	@echo "  docs-clean             Clean documentation build"
	@echo "  sdist                  Build source distribution"

.PHONY: fix
fix:
	uv run hatch run quality:format

.PHONY: check
check:
	uv run hatch run quality:check
	uv run hatch run quality:typecheck

.PHONY: test
test:
	uv run hatch run tests:run

.PHONY: docs
docs:
	uv run hatch run docs:build

.PHONY: docs-clean
docs-clean:
	uv run hatch run docs:clean

.PHONY: sdist
sdist:
	uv run hatch run build:run
