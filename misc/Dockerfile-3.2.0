# docker build -t sundialsml-env:3.2.0 - <Dockerfile-3.2.0
# docker run -v $HOME/Projects/sundialsml:/home/opam/sundialsml -it sundialsml-env:3.2.0 bash
FROM ocaml/opam2
ENV SUNDIALSML_SRC=https://raw.githubusercontent.com/inria-parkas/sundialsml/master\
 SUPERLUMT_TARGZ=http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_mt_3.1.tar.gz\
 SUNDIALS_GIT=https://github.com/LLNL/sundials.git\
 SUNDIALS_TAG=v3.2.0\
 SUNDIALSML_GIT=https://github.com/inria-parkas/sundialsml\
 BUILDDIR=$HOME
USER root
RUN apt-get update\
 && apt-get install -y\
      cmake liblapack-dev libopenmpi-dev openmpi-bin libsuitesparse-dev python\
      wget csh\
      pkg-config m4\
 && rm -rf /var/lib/apt/lists/*
USER opam
RUN opam update && opam install mpi
USER root
RUN wget -qO- $SUPERLUMT_TARGZ | tar xz -C $BUILDDIR\
 && wget -qO- $SUNDIALSML_SRC/misc/superlu_mt_3.1.patch\
	| patch -d $BUILDDIR/SuperLU_MT_3.1 -p1\
 && make -C $BUILDDIR/SuperLU_MT_3.1 install lib\
 && mkdir -p /usr/local/include/superlu_mt\
 && cp $BUILDDIR/SuperLU_MT_3.1/SRC/*.h /usr/local/include/superlu_mt/\
 && cp $BUILDDIR/SuperLU_MT_3.1/lib/* /usr/local/lib/\
 && rm -rf $BUILDDIR/SuperLU_MT_3.1
#-DLAPACK_ENABLE=1 -DLAPACK_LIBRARIES='-llapack -lblas'
RUN git clone $SUNDIALS_GIT --branch $SUNDIALS_TAG --depth 1 $BUILDDIR/sundials\
 && wget -qO- $SUNDIALSML_SRC/misc/sundials-3.2.0.patch\
	| patch -d $BUILDDIR/sundials -p1\
 && mkdir $BUILDDIR/sundials-build\
 && cd $BUILDDIR/sundials-build\
 && cmake -Wno-dev $BUILDDIR/sundials\
 -DCMAKE_BUILD_TYPE=Release\
 -DOPENMP_ENABLE=1\
 -DMPI_ENABLE=1\
 -DMPI_C_COMPILER=/usr/bin/mpicc\
 -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec\
 -DBLAS_ENABLE=1\
 -DPTHREAD_ENABLE=1\
 -DEXAMPLES_ENABLE_C=1\
 -DKLU_ENABLE=1\
 -DKLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu\
 -DKLU_INCLUDE_DIR=/usr/include/suitesparse\
 -DSUPERLUMT_ENABLE=1\
 -DSUPERLUMT_THREAD_TYPE=PTHREAD\
 -DSUPERLUMT_LIBRARY_DIR=/usr/local/lib\
 -DSUPERLUMT_INCLUDE_DIR=/usr/local/include/superlu_mt\
 -DSUPERLUMT_LIBRARIES=-lblas\
 && make -j\
 && make install\
 && printf "#!/bin/sh\nexec ./configure SUPERLUMT_INCLUDE_DIR=/usr/local/include/superlu_mt\n"\
	> $HOME/run-sundialsml-configure\
 && chmod ugo+x $HOME/run-sundialsml-configure
USER opam
ENV LD_LIBRARY_PATH=/usr/local/lib
#RUN git clone $SUNDIALSML_GIT $HOME/sundialsml
LABEL sundials-tag=$SUNDIALS_TAG
WORKDIR $HOME/sundialsml
