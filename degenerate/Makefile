CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: selftest libdegenerate.so

selftest: src/selftest.c src/degenerate.c src/degenerate.h
	$(CC) $(CFLAGS) src/selftest.c src/degenerate.c -o $@ -lm

libdegenerate.so: src/degenerate.c src/degenerate.h
	$(CC) $(CFLAGS) -fPIC -shared src/degenerate.c -o $@ -lm

test: all
	./selftest

clean:
	rm -f selftest libdegenerate.so
