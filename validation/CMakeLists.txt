SET(MESH_MED
  ../tests/ressources/meshSquare.med
  ../tests/ressources/meshCube.med
  ../tests/ressources/meshSphere.med
  )
SET(MESH_MED2
  ../tests/ressources/2DTriangles
  ../tests/ressources/2DLocRefinedSquares
  ../tests/ressources/3DTetrahedron
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/validation)

file(GLOB NICE_EXAMPLES_TO_INSTALL thermique1d transport1d transport2d_ns transport2d_s Poisson2DEF_ns Poisson3DEF_ns Poisson2DVF_ns Poisson2DVF_s PoissonCurvedSurfaceEF_ns)

install(DIRECTORY ${NICE_EXAMPLES_TO_INSTALL} DESTINATION share/Validation)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    #ADD_SUBDIRECTORY(PoissonCurvedSurfaceEF_ns)
    #ADD_SUBDIRECTORY(Poisson2DVF_s)
    #ADD_SUBDIRECTORY(Poisson2DEF_ns)
    #ADD_SUBDIRECTORY(Poisson3DEF_ns)
ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_TEST(validation ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validation.py)
    SET_PROPERTY(TEST validation PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib")
    SET_PROPERTY(TEST validation APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/Poisson2DEF_ns:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}")
endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


