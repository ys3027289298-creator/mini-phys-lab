CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: solve

solve: src/solve.c src/solve.h
	$(CC) $(CFLAGS) src/solve.c -o $@ -lm

test: solve
	./solve

contract: test/contract.c
	$(CC) $(CFLAGS) test/contract.c -o contract
	./contract src ../integrate/src

clean:
	rm -f solve contract

.PHONY: all test contract clean
