.PHONY: help test test-r test-python fixture image run clean

SERVICE_NAME ?= processor-diffcyt
IMAGE ?= processor-diffcyt:dev
PYTHON ?= python3
# Rscript with diffcyt installed runs the R tests directly; otherwise they run in the image.
HAVE_R := $(shell Rscript -e 'quit(status = !requireNamespace("diffcyt", quietly = TRUE))' >/dev/null 2>&1 && echo yes)

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - R end-to-end tests on the synthetic fixture (local Rscript, else in the image) + Lambda handler tests"
	@echo "make fixture  - write the synthetic fixture (planted DA/DS effects) to data/fixture/{linked,flat}/"
	@echo "make image    - build the image (linux/amd64)"
	@echo "make run      - build the image and run it with docker compose over data/fixture/linked/ into data/output/"
	@echo "                (INPUT=/path/to/input to run on other data, e.g. flowsom + prep outputs side by side)"
	@echo "make clean    - remove data/output/"

test: test-python test-r

test-python:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v

test-r:
ifeq ($(HAVE_R),yes)
	Rscript tests/test_main.R
else
	$(MAKE) image
	docker run --rm --platform linux/amd64 --memory 3g -v "$(CURDIR):/src:ro" -w /src --entrypoint Rscript $(IMAGE) tests/test_main.R
endif

fixture:
ifeq ($(HAVE_R),yes)
	Rscript tests/make_fixture.R data/fixture both
else
	$(MAKE) image
	mkdir -p data
	docker run --rm --platform linux/amd64 --user "$$(id -u):$$(id -g)" -v "$(CURDIR):/src" -w /src --entrypoint Rscript $(IMAGE) tests/make_fixture.R data/fixture both
endif

image:
	docker build --platform linux/amd64 -t $(IMAGE) .

run: clean
	@if [ -z "$(INPUT)" ]; then $(MAKE) fixture; fi
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	INPUT_HOST="$(or $(INPUT),./data/fixture/linked)" docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
