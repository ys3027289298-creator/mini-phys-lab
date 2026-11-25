CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
PYTHON ?= python3

all: integrate

integrate: src/integrate.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: integrate
	./integrate

# penetration-vs-truncation curves; writes artifacts/substep_curves.png + .csv
plot:
	$(PYTHON) -I tools/substep_study.py

.PHONY: all test plot clean

clean:
	rm -f integrate
