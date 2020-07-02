Sundials/ML Development Notes
=============================

Debian
------

Packages
- `cmake`
- `libopenmpi-dev` (for parallel nvectors)
- `libsuitesparse-dev` (for KLU)

### SuperLU/MT

Download the source from 
[https://portal.nersc.gov/project/sparse/superlu/#superlu_mt](SuperLU_MT 
Version 3.1). Edit `make.inc` by commenting out `BLASDEF=...` and setting
`BLASLIB=../lib/libblas$(PLAT).a`. Run the following commands.

```
export SUPERLUMT_DIR=/path/to/superlumt/download
make blaslib
make # ignore the implicit function declarations
```

### Building Sundials

Create directories `sundials-build` and `sundials-install`. Change to the 
former. Run the following commands.

```
SUNDIALS_SOURCE=/path/to/sundials/source
export SUNDIALS_DIR=/path/to/sundials-install
cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="${SUNDIALS_DIR}" \
    -Wno-dev "$SUNDIALS_SOURCE" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DOPENMP_ENABLE=1 \
    -DPTHREAD_ENABLE=1 \
    -DMPI_ENABLE=1 \
    -DKLU_ENABLE=1 \
    -DKLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu \
    -DKLU_INCLUDE_DIR=/usr/include/suitesparse \
    -DSUPERLUMT_ENABLE=1 \
    -DSUPERLUMT_LIBRARY_DIR="$SUPERLUMT_DIR/lib" \
    -DSUPERLUMT_INCLUDE_DIR="$SUPERLUMT_DIR/SRC" \
    -DSUPERLUMT_LIBRARIES=-lblas
make -j install
```

It should normally be possible to add the following lines to build with 
lapack.

```
-DLAPACK_ENABLE=1 -DLAPACK_LIBRARIES='-llapack -lblas'`
```

But this currently (buster) results in an error message during Sundials 
build about 64-bit integer indexes.

### Building Sundials/ML

Run the following command. This command assumes that `SUPERLUMT_DIR` and 
`SUNDIALS_DIR` have been exported as specified above (otherwise they must be 
specified on the command-line).

```
export LD_LIBRARY_PATH="$SUNDIALS_DIR/lib:$LD_LIBRARY_PATH"
./configure --enable-debug
make -j
```

