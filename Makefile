.PHONY: all test clean

DIRS = collide solve integrate diagnose

all test:
	@for d in $(DIRS); do $(MAKE) -C $$d $@ || exit $$?; done

clean:
	@for d in $(DIRS); do $(MAKE) -C $$d clean; done
