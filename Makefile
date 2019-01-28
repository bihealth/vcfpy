.PHONY: default black flake8

default: black flake8

black:
	black -l 100 .

flake8:
	flake8 .

