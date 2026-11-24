CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: collide

collide: src/collide.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: collide
	./collide

clean:
	rm -f collide
