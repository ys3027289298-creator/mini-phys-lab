CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: solve

solve: src/solve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: solve
	./solve

clean:
	rm -f solve
