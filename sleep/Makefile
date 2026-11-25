CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Isrc

.PHONY: all test sweep clean

all: sleep

sleep: src/sleep.c src/sleep.h
	$(CC) $(CFLAGS) $(CPPFLAGS) src/sleep.c -o $@

test: sleep
	./sleep

sweep: sleep
	mkdir -p ../artifacts
	./sleep --sweep ../artifacts/sleep_tradeoff.csv

clean:
	rm -f sleep
