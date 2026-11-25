CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

BIN := out/inclined
DATA := out/curves.csv out/tiers.csv out/timeseries.csv out/sweep.csv

all: $(BIN)

$(BIN): src/inclined.c
	@mkdir -p out
	$(CC) $(CFLAGS) $< -o $@ -lm

test: $(BIN)
	./$(BIN) self

# regenerate every CSV the report scripts read
data: $(BIN)
	./$(BIN) curves      > out/curves.csv
	./$(BIN) tiers       > out/tiers.csv
	./$(BIN) timeseries  > out/timeseries.csv
	./$(BIN) sweep       > out/sweep.csv

figures: data
	python3 report/plot.py

mc: data
	python3 report/mc_check.py

report: figures mc

clean:
	rm -f $(BIN) $(DATA)

.PHONY: all test data figures mc report clean
