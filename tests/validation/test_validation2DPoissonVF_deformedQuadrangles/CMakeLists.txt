
SET(SCRIPT
  ../scripts/Poisson2DVF/FiniteVolumes2DPoisson_SQUARE.py
  test_validation2DPoissonVF_deformedQuadrangles.py
  )
FILE(GLOB MY_MESHES_AND_PICTURES
  "../../ressources/2DDeformedQuadrangles/squareWithDeformedQuadrangles_*"
)
SET(NOTEBOOK ../../doc/2DPoissonVF_deformedQuadrangles/Convergence_Poisson_FV5_SQUARE_deformedQuadrangles.ipynb )


file(COPY ${SCRIPT} ${NOTEBOOK} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)

    ADD_TEST(validationPoisson_2DFV_deformedQuadrangles ${PYTHON_EXECUTABLE} test_validation2DPoissonVF_deformedQuadrangles.py)

    SET_PROPERTY(TEST validationPoisson_2DFV_deformedQuadrangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")

    SET_PROPERTY(TEST validationPoisson_2DFV_deformedQuadrangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


