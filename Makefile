all:
	$(MAKE) -C src
test:
	$(MAKE) -C examples tests.opt.log
