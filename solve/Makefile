CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: solve bench

solve: src/solve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

bench: src/bench_solve.c
	$(CC) $(CFLAGS) $< -o bench_solve -lm

test: solve
	./solve

bench-test: bench
	./bench_solve selftest

bench-run: bench
	./bench_solve

clean:
	rm -f solve bench_solve
