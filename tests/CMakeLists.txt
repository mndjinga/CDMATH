# nom du projet tests
PROJECT(tests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

ADD_SUBDIRECTORY(cdmath)
ADD_SUBDIRECTORY(ressources)

IF (CDMATH_WITH_PYTHON)
    ADD_SUBDIRECTORY(swig)
ENDIF (CDMATH_WITH_PYTHON)

FILE(COPY doc/ DESTINATION convergence)
INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/convergence DESTINATION share/doc/)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
  # Examples tests
  set (examples_dir ${CDMATH_SOURCE_DIR}/tests/examples)
  add_subdirectory (${examples_dir})
  add_custom_target (examples COMMAND ctest -R Example)# may be replace ctest -R with ctest -L

  #validation tests
  add_subdirectory (validation)
  add_custom_target (validation COMMAND ctest -R validation)# may be replace ctest -R with ctest -L
  add_custom_target (check COMMAND ctest -E validation)# may be replace ctest -E with ctest -LE
endif ()


INCLUDE_DIRECTORIES(
  ${BASE_DIR}/inc
  ${MESH_DIR}/inc
  ${MED_COUPLING_DIR}
  ${PRE_REQUIS_DIR}/INTERP_KERNEL
  ${PRE_REQUIS_DIR}/INTERP_KERNEL/Bases
  cdmath
  )

SET(_extra_lib_TESTS)
if (CDMATH_WITH_PETSC)
    INCLUDE_DIRECTORIES(
      ${LINEARSOLVER_DIR}/inc
      ${PETSC_INCLUDES}
      )
    SET(_extra_lib_TESTS linearsolver ${PETSC_LIBRARIES})
endif (CDMATH_WITH_PETSC)

SET(MAIN_SOURCES
  MainTests.cxx
  )

ADD_EXECUTABLE(TestCDMATH ${MAIN_SOURCES})

if (CDMATH_WITH_PETSC)
   SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS -DCDMATH_WITH_PETSC)
else(CDMATH_WITH_PETSC)
   SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS "")
endif (CDMATH_WITH_PETSC)

TARGET_LINK_LIBRARIES(TestCDMATH TestMesh base mesh medcoupling medloader ${CPPUNIT_LIBS} ${_extra_lib_TESTS})

ADD_TEST(UnitTestCDMATH TestCDMATH)
SET_PROPERTY(TEST UnitTestCDMATH PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib")

INSTALL(TARGETS TestCDMATH DESTINATION bin/tests)
