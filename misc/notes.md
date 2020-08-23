Release procedure
-----------------

1. Update and commit performance graph if necessary:
```
make doc/html/perf.opt.png
```

2. Test, test, test!
```
make distcheck
```

3. Branch from trunk if necessary for major version.
```
git checkout -b 2.5.0-release
```

4. Check download links (to Sundials page, etc.)
   README.md, intro.doc, config.in

5. Check version info in configure
```
./configure --help | head -n 1
```

6. Tag release in release branch with release notes in commit message:
```
git tag -a v2.5.0p0 --sign
git push origin v2.5.0p0
```

7. Update gh-pages:
   a. ensure MathJax is http://
   b. current docs (root level)
   c. add a subdirectory for the release

8. OPAM:
   - update package files in forked opam repository
   - test
   - submit pull request

9. Github: check release text

10. Announce on: `caml-list@inria.fr`

Debugging Tips
--------------

* Configure with `--enable-debug`.
* Compile OCaml programs with `-runtime-variant d` to enable debugging 
  assertions in the OCaml runtime.
* Use `valgrind --track-origins=yes` to run OCaml programs.
* Use `ltrace -l 'libsundials_*'
* Use the `misc/ltrace_mpi_wrapper` to debug mpi problems:
  `mpirun -np 4 ../../../misc/ltrace_mpi_wrapper -l 'libsundials_*' ./prog.opt`
* Debug with both `valgrind` and `gdb`:
```
valgrind --track-origins=yes --vgdb=yes --vgdb-error=0 ./blah.opt
gdb ./blah.opt
```

