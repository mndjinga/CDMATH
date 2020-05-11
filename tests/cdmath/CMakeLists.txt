INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${PRE_REQUIS_DIR}/MEDLoader
  ${PRE_REQUIS_DIR}/MEDCoupling
  ${PRE_REQUIS_DIR}/INTERP_KERNEL
  ${PRE_REQUIS_DIR}/INTERP_KERNEL/Bases
  ${MESH_DIR}/inc 
  ${BASE_DIR}/inc
  )

if (CDMATH_WITH_PETSC)
    INCLUDE_DIRECTORIES(
      ${PETSC_INCLUDES}
      ${LINEARSOLVER_DIR}/inc
      )
else (CDMATH_WITH_PETSC)
endif (CDMATH_WITH_PETSC)


SET(MESH_MED
  ../ressources/meshSquare.med
  ../ressources/meshSphere.med
  ../ressources/meshCube.med
	)

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/tests) 

if (CDMATH_WITH_PETSC)
    FILE(GLOB TestMesh_SOURCES . *.cxx)
    ADD_LIBRARY(TestMesh SHARED ${TestMesh_SOURCES})
    TARGET_LINK_LIBRARIES(TestMesh base mesh linearsolver ${HDF5_LIBRARIES} medloader ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
else (CDMATH_WITH_PETSC)  
    FILE(GLOB TestMesh_SOURCES . CellTests.cxx DoubleTabTests.cxx FaceTests.cxx FieldTests.cxx IntTabTests.cxx MatrixTests.cxx MeshTests.cxx NodeTests.cxx PointTests.cxx VectorTests.cxx)
    ADD_LIBRARY(TestMesh SHARED ${TestMesh_SOURCES})
    TARGET_LINK_LIBRARIES(TestMesh base mesh ${HDF5_LIBRARIES} medloader)
endif (CDMATH_WITH_PETSC) 
  

INSTALL(TARGETS TestMesh DESTINATION lib)
