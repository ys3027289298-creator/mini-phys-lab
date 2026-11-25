MODULES = collide solve integrate sleep

.PHONY: all test clean

all test clean:
	@for d in $(MODULES); do $(MAKE) -C $$d $@; done
