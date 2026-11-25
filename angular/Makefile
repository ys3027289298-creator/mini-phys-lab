CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Iinclude
LDFLAGS ?= -lm

ART ?= ../artifacts

.PHONY: all test bench bench-run plot clean

all: angular

angular: src/angular.c src/am.c include/am.h
	$(CC) $(CFLAGS) src/angular.c src/am.c -o $@ $(LDFLAGS)

test: angular
	./angular

bench: bench/am_bench.c src/am.c include/am.h
	$(CC) $(CFLAGS) -O2 bench/am_bench.c src/am.c -o am_bench $(LDFLAGS)

bench-run: bench
	mkdir -p $(ART)
	./am_bench $(ART)

plot:
	python3 tools/plot_tradeoff.py $(ART)

clean:
	rm -f angular am_bench
