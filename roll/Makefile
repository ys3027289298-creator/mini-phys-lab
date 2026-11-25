CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: roll

roll: src/roll.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: roll
	./roll

clean:
	rm -f roll

.PHONY: all test clean
