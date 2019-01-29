.PHONY: default black flake8 test test-v test-vv

default: black flake8

black:
	black -l 100 .

black-check:
	black -l 100 --check .

flake8:
	flake8 .

test:
	pytest

test-v:
	pytest -v

test-vv:
	pytest -vv
