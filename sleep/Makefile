CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
T ?= 0.045

REPO_ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ART := $(REPO_ROOT)/../artifacts

.PHONY: all test sweep grades bench knee derive experiments plot clean

all: sleep

sleep: src/sleep.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sleep
	./sleep test

derive: sleep
	./sleep derive

knee: sleep
	./sleep knee

bench: sleep
	./sleep bench

# threshold sweep CSV + official three-grade JSON run, into ../artifacts
sweep: sleep
	mkdir -p $(ART)
	./sleep sweep $(ART)/sleep_sweep.csv
	./sleep grades $(T) 2>/dev/null > $(ART)/sleep_grade_runs.json
	./sleep knee > $(ART)/sleep_knee.csv
	./sleep bench | tee $(ART)/sleep_bench.txt

plot: sweep
	python3 tools/plot_sweep.py $(ART)/sleep_sweep.csv $(ART)/sleep_grade_runs.json $(ART)

# full reproducible evidence bundle: tests, data, curves, benchmark
experiments: test plot
	@echo "experiments complete -> $(ART)"

clean:
	rm -f sleep
