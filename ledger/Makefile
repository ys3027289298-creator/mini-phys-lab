CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ledger

ledger: src/ledger.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ledger
	./ledger

clean:
	rm -f ledger
