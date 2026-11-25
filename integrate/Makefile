CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: integrate

integrate: src/integrate.c src/integrate.h
	$(CC) $(CFLAGS) src/integrate.c -o $@ -lm

test: integrate
	./integrate

clean:
	rm -f integrate
