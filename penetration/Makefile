CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS ?=

SRC = src/penetration.c
HDR = src/pen_depth.h src/pen_model.h

all: penetration

penetration: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -Isrc $(SRC) -o $@ -lm

sweep: tools/sweep.c $(HDR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -Isrc tools/sweep.c -o $@ -lm

test: penetration
	./penetration

curves: sweep
	./sweep > ../docs/penetration-curves.csv

clean:
	rm -f penetration sweep

.PHONY: all test curves clean
