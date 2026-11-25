CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: chain

chain: src/chain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: chain
	./chain

clean:
	rm -f chain
