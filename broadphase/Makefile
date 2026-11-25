CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: broadphase

broadphase: src/broadphase.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: broadphase
	./broadphase test

sweep: broadphase
	./broadphase sweep results.csv

clean:
	rm -f broadphase results.csv
