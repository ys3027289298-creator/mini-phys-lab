CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: integrate

integrate: src/integrate.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: integrate
	./integrate

clean:
	rm -f integrate
