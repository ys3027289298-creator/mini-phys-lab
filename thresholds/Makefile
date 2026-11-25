PY ?= python3

.PHONY: all derive test plots clean

all: derive test

derive:
	$(PY) -I src/derive.py

test:
	$(PY) -I src/experiments.py

plots:
	$(PY) -I src/plots.py

clean:
	rm -f results/report.json results/calibration.csv results/*.png
