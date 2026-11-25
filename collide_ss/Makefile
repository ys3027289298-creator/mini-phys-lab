CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -Isrc

.PHONY: all test repro sweep plot clean

all: test repro

collide_ss.o: src/collide_ss.c src/collide_ss.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/collide_ss.c -o $@

test: collide_ss.o
	$(CC) $(CFLAGS) $(CPPFLAGS) src/test_collide_ss.c collide_ss.o -o test_collide_ss -lm
	./test_collide_ss

# 复现制品：朴素实现的 SIGFPE 崩溃 + 默认 IEEE 语义下的 NaN 污染。
# 只负责构建；崩溃复现请手动执行 ./naive_repro（预期 exit=136）。
repro:
	$(CC) $(CFLAGS) src/naive_repro.c -o naive_repro -lm
	$(CC) $(CFLAGS) src/nan_demo.c -o nan_demo -lm

# 阈值扫描：SWEEP_ALLOW_ZERO_EPS 让 eps=0 表达“无守卫”基线
sweep: collide_ss.o
	mkdir -p results
	$(CC) $(CFLAGS) $(CPPFLAGS) -DSWEEP_ALLOW_ZERO_EPS tools/sweep.c collide_ss.o -o sweep -lm
	./sweep results/sweep.csv

plot: sweep
	python3 -I tools/plot_sweep.py results/sweep.csv results/threshold_curves.png

clean:
	rm -f collide_ss.o test_collide_ss naive_repro nan_demo sweep
	rm -rf results
