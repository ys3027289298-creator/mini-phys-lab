CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: rollfeel

rollfeel: src/rollfeel.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: rollfeel
	./rollfeel

clean:
	rm -f rollfeel
