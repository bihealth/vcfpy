.PHONY: default
default: help

.PHONY: help
help:
	@echo "Available targets:"
	@echo "  help                   Show this help message"
	@echo "  fix                    Format source code"
	@echo "  check                  Run checks"
	@echo "  test                   Format source code"

.PHONY: fix
fix:
	uv run hatch run quality:format

.PHONY: check
check:
	uv run hatch run quality:check
# 	uv run hatch run quality:typecheck

.PHONY: test
test:
	uv run hatch run tests:run -- --record-mode=none
