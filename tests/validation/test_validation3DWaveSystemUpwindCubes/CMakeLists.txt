
SET(SCRIPT
  ../scripts/WaveSystemUpwind/WaveSystemUpwind.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_BINARY_DIR}/tests/validation/test_validation3DWaveSystemUpwindCubes)
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation3DWaveSystemUpwindCubes)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_TEST(test_validation3DWaveSystemUpwindCubes ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation3DWaveSystemUpwindCubes.py)
    SET_PROPERTY(TEST test_validation3DWaveSystemUpwindCubes PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib")
    SET_PROPERTY(TEST test_validation3DWaveSystemUpwindCubes APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}")
endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


