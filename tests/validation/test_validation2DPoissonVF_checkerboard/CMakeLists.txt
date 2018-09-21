
SET(SCRIPT
  ../scripts/Poisson2DVF/FiniteVolumes2DWithCDMATH.py
  )
SET(PICTURES
  ../../ressources/2DCheckerboard/checkerboard_4x4.png
  ../../ressources/2DCheckerboard/checkerboard_8x8.png
  ../../ressources/2DCheckerboard/checkerboard_16x16.png
  ../../ressources/2DCheckerboard/checkerboard_32x32.png
  ../../ressources/2DCheckerboard/checkerboard_64x64.png
)
SET(NOTEBOOK
  Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb
)

file(COPY ${SCRIPT} ${PICTURES} ${NOTEBOOK} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} ${PICTURES} ${NOTEBOOK} DESTINATION share/validation/test_validation2DPoissonVF_checkerboard)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_TEST(validationPoisson_2DFV_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DPoissonVF_checkerboard.py)
    SET_PROPERTY(TEST validationPoisson_2DFV_checkerboard PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validationPoisson_2DFV_checkerboard APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")
endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


