.PHONY: all
all: black

.PHONY: black
black:
	black -l 100 .

.PHONY: test
test:
	black -l 100 --check .
