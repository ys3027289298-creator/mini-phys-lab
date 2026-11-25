CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
PYTHON ?= python3

all: solve

solve: src/solve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: solve
	./solve

# friction-cone design curves (writes artifacts/friction_cone_ratios.png)
plot:
	$(PYTHON) -I tools/friction_cone.py

clean:
	rm -f solve
