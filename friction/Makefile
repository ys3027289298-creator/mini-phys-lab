CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: friction

friction: src/friction.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: friction
	./friction selftest

clean:
	rm -f friction
