CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: deform

deform: src/deform.c src/deform.h src/video_measurements.h
	$(CC) $(CFLAGS) src/deform.c -o $@ -lm

test: deform
	./deform

sweep: deform
	./deform --sweep ../artifacts/stiffness_sweep.csv

clean:
	rm -f deform
