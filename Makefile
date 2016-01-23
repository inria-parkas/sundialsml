ifeq ($(wildcard config),)
$(error Please run ./configure first)
endif


all:
	$(MAKE) -C src
test:
	$(MAKE) -C examples tests.opt.log
doc:
	$(MAKE) -C doc

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

distclean:
	$(MAKE) -C doc $@
	$(MAKE) -C src $@
	$(MAKE) -C examples $@
	-@$(RM) -f config config.log

cleandoc:
	$(MAKE) -C doc clean

.PHONY: all test doc clean distclean cleandoc install install-sys install-doc \
	install-ocamlfind install-findlib
