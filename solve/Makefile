CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: solve order_compare

solve: src/solve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

order_compare: src/order_compare.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: solve order_compare
	./solve
	mkdir -p ../artifacts/friction_order
	./order_compare > ../artifacts/friction_order/c_five_cases.csv

plot: order_compare
	mkdir -p ../artifacts/friction_order
	./order_compare > ../artifacts/friction_order/c_five_cases.csv
	python3 order_experiment.py

clean:
	rm -f solve order_compare
