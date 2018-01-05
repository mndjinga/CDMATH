#!/bin/bash

export CDMATH_DIR=@CMAKE_INSTALL_PREFIX@
export PETSC_DIR=@PETSC_DIR@
export PETSC_ARCH="@PETSC_ARCH@"

export LD_LIBRARY_PATH=$CDMATH_DIR/lib/:$CDMATH_DIR/lib/medcoupling:$CDMATH_DIR/lib/med:$PETSC_DIR/$PETSC_ARCH/lib:${PETSC_DIR}/lib:/usr/lib64/:${LD_LIBRARY_PATH}
export PYTHONPATH=$CDMATH_DIR/lib/:$CDMATH_DIR/lib/cdmath:$CDMATH_DIR/lib/medcoupling:$CDMATH_DIR/lib/med:$CDMATH_DIR/bin/cdmath:${PYTHONPATH}

