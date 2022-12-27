ifeq ($(wildcard config),)
$(error Please run ./configure first)
endif

include ./config

all:
	$(MAKE) -C src
test:
	$(MAKE) -C examples tests.opt.log
doc:
	$(MAKE) -C src
	$(MAKE) -C doc

loc:
	@printf "OCaml (ocamlwc: lines of code, lines of comments)\n"
	@ocamlwc $(shell find src -iregex '.*ml$$' -o -iregex '.*mli$$') \
	    | tail -1
	@printf "Running cloc:\n"
	@cloc --quiet src

### Preparing Releases

# Sets up a sandbox and makes sure everything compiles.  You need
# everything to do this: gcc, git, OCaml >= 4, OCamlMPI, lapack, etc.
# The sandbox is configured with $(CONFIG_COMMAND), which defaults to
# the command you used to configure the current repository.
distcheck:
	git clone . sandbox
	cd sandbox && $(CONFIG_COMMAND) > config.log
	@ # Make sure everything is enabled.
	! grep "NOT FOUND" sandbox/config.log
	! grep -i "DISABLED" sandbox/config.log
	! grep "without lapack" sandbox/config.log
	@ # Tests.  Important ones and short ones first.
	$(MAKE) -C sandbox all
	$(MAKE) -C sandbox/examples tests.opt.log
	$(MAKE) -C sandbox/examples/ocaml
	$(MAKE) -C sandbox doc
	@ # Because perf.opt.log repeats tests, it sometimes uncovers
	@ # GC bugs that tests.opt.log doesn't.
	$(MAKE) -C sandbox/examples PERF_DATA_POINTS=1 perf.opt.log
	$(MAKE) -C sandbox/examples tests.byte.log
	rm -rf sandbox

# install-sys installs to OCaml's system directory -- /usr/lib/ocaml
# on Debian derivatives.

# install-doc installs the doc in $(DOCDIR), which is set by
# configure.

# install-ocamlfind installs via ocamlfind.  install-findlib is an
# alias for install-ocamlfind.

# install performs install-sys and install-doc.
install install-sys install-doc install-ocamlfind install-findlib:
	$(MAKE) -C src $@

uninstall uninstall-sys uninstall-doc uninstall-ocamlfind uninstall-findlib:
	$(MAKE) -C src $@

clean:
	$(MAKE) -C src $@
	$(MAKE) -C examples $@

cleanall distclean:
	$(MAKE) -C doc $@
	$(MAKE) -C src $@
	$(MAKE) -C examples $@
	-@$(RM) -f config config.log

cleandoc:
	$(MAKE) -C doc clean

.PHONY: all test doc clean cleanall distclean cleandoc \
	install install-sys install-doc \
	install-ocamlfind install-findlib
