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

Performance Graphs
------------------
Requires _octave_ with statistics package
(`apt install octave octave-statistics`)

1. Compile Sundials with `-DCMAKE_BUILD_TYPE=Release`
2. Compile Sundials/ML with `--unsafe`
```
make distclean
./configure --unsafe --disable-openmp \
  SUNDIALS_DIR=<sundials-5.8.0-install> \
  SUPERLUMT_DIR=<SuperLU_MT_3.1> \
  KLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu \
  KLU_INCLUDE_DIR=/usr/include/suitesparse
make -j
cd examples
make perf-intv.opt.pdf GC_AT_END=1 PERF_DATA_POINTS=40
make SIZE=2000,1200 FONT=Arial,14 DOTSIZE=1 perf-intv.opt.pngcairo
cp perf-intv.opt.pngcairo ../doc/html/perf-intv.opt.png
cp perf-intv.opt.log ../doc/html/perf-intv.opt.log
cp perf-intv.opt.pdf ../doc/html/perf-intv.opt.pdf
```
3. Compile Sundials/ML normally
```
make distclean
./configure --disable-openmp \
  SUNDIALS_DIR=<sundials-5.8.0-install> \
  SUPERLUMT_DIR=<SuperLU_MT_3.1> \
  KLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu \
  KLU_INCLUDE_DIR=/usr/include/suitesparse
make -j
cd examples
make perf-intv.opt.pdf GC_AT_END=1 PERF_DATA_POINTS=40
make SIZE=2000,1200 FONT=Arial,14 DOTSIZE=1 perf-intv.opt.pngcairo
cp perf-intv.opt.pngcairo ../doc/html/perf-intv-safe.opt.png
cp perf-intv.opt.log ../doc/html/perf-intv-safe.opt.log
cp perf-intv.opt.pdf ../doc/html/perf-intv-safe.opt.pdf
```

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

