
SET(SCRIPT
  ../scripts/WaveSystemUpwind/WaveSystemUpwind.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation3DWaveSystemUpwindTetrahedra)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_TEST(test_validation3DWaveSystemUpwindTetrahedra ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindTetrahedra.py)
    SET_PROPERTY(TEST test_validation3DWaveSystemUpwindTetrahedra PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST test_validation3DWaveSystemUpwindTetrahedra APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")
endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


