all:
	$(MAKE) -C src
test:
	$(MAKE) -C examples tests.opt.log

install: install-sys $(if $(INSTALL_DOCS),install-doc)

# install-sys installs to OCaml's system directory -- /usr/lib/ocaml
# on Debian derivatives.

# install-doc installs the doc in $(DOCDIR), which is set by
# configure.

# install-ocamlfind installs via ocamlfind.  install-findlib is an
# alias for install-ocamlfind.
install-sys install-doc install-ocamlfind install-findlib:
	$(MAKE) -C src $@

uninstall uninstall-sys uninstall-doc uninstall-ocamlfind uninstall-findlib:
	$(MAKE) -C src $@
