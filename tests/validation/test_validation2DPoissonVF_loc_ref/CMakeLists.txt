
SET(SCRIPT
  ../scripts/Poisson2DVF/FiniteVolumes2DWithCDMATH.py
  )
SET(PICTURES
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_1.png
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_2.png
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_3.png
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_4.png
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_5.png
  ../../ressources/2DLocRefinedSquares/squareWithLocRefSquares_6.png
)
SET(NOTEBOOK
  Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb
)

file(COPY ${SCRIPT} ${PICTURES} ${NOTEBOOK} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} ${PICTURES} ${NOTEBOOK} DESTINATION share/validation/test_validation2DPoissonVF_loc_ref)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_TEST(validationPoisson_2DFV_loc_ref ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation2DPoissonVF_loc_ref.py)
    SET_PROPERTY(TEST validationPoisson_2DFV_loc_ref PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validationPoisson_2DFV_loc_ref APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")
endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


