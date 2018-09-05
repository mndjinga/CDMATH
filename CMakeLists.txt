cmake_minimum_required (VERSION 2.8)

# Project name
project (CDMATH)
set (CDMATH_VERSION_MAJOR 1)
set (CDMATH_VERSION_MINOR 0)

# Project options
option (CDMATH_WITH_PETSC "Compile CDMATH with PETSc linking." OFF)
option (CDMATH_WITH_DOCUMENTATION "Generate documentation with doxygen." OFF)
option (CDMATH_WITH_PYTHON "Compile Python interface for CDMATH." OFF)
option (CDMATH_WITH_POSTPRO "Install postprocessing Python modules." OFF)
option (CDMATH_WITH_PACKAGE "Generate RPM, Debian and tarball packages." OFF)
option (CDMATH_WITH_TESTS "Compile unit testing." OFF)

# Base directories
set (BASE_DIR ${CDMATH_SOURCE_DIR}/base)
set (MESH_DIR ${CDMATH_SOURCE_DIR}/mesh)
set (LINEARSOLVER_DIR ${CDMATH_SOURCE_DIR}/linearsolver)
set (CDMATH_SWIG_DIR ${CDMATH_SOURCE_DIR}/swig)
set (CDMATH_POSTPRO_DIR ${CDMATH_SOURCE_DIR}/postprocessing)
set (PRE_REQUIS_DIR ${CDMATH_SOURCE_DIR}/pre_requis)
set (INTERP_KERNEL_DIR ${PRE_REQUIS_DIR}/INTERP_KERNEL)
set (MED_COUPLING_DIR ${PRE_REQUIS_DIR}/MEDCoupling)
set (MED_COUPLING_SWIG_DIR ${PRE_REQUIS_DIR}/MEDCoupling_Swig)
set (MED_LOADER_DIR ${PRE_REQUIS_DIR}/MEDLoader)
set (MED3_DIR ${PRE_REQUIS_DIR}/MEDFile)
set (MEDFILE_INCLUDE_DIRS ${PRE_REQUIS_DIR}/MEDIncludes)

set (TESTS_DIR ${CDMATH_SOURCE_DIR}/tests)
list (APPEND CMAKE_MODULE_PATH "${CDMATH_SOURCE_DIR}/cmake_files")

set (MACHINE PCLINUX)
set (MED3_DEFINITIONS "-D${MACHINE}")
set (MED_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/MEDFile)
set (MED3_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/pre_requis/MEDFile/include ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/MEDFile/include)
set (MEDFILE_C_LIBRARIES med medC)

# HDF5
find_package (HDF5 REQUIRED)
set (CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-dev")
set (CPACK_RPM_PACKAGE_REQUIRES "hdf5-devel")

# PETSc
if (CDMATH_WITH_PETSC)
    message ( STATUS "Checking variable PETSC_DIR : " $ENV{PETSC_DIR} )

  if ( DEFINED ENV{PETSC_DIR})

    find_package (PETSc 3.4 REQUIRED)
    petsc_get_version ()

  else ()
    message ( STATUS "PETSC_DIR not set, searching petsc in the system" )

    if ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so")  #Case fedora/redhat system install
      message ( STATUS "PETSC includes found in /usr/include/petsc/" )
      message ( STATUS "PETSC library found in /usr/lib64" )
      set(PETSC_DIR /usr/)
      set(PETSC_INCLUDES  /usr/include/petsc /usr/include/petsc/petsc/mpiuni)
      set(PETSC_LIBRARIES /usr/lib64/libpetsc.so)

      set(PETSC_VERSION "3.7") #value for fedora 25/26 
      set (CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, petsc-devel (>= 3.4)") # This is not fully working yet. 

    #elseif ( IS_DIRECTORY "/usr/lib/petsc/") #Case ubuntu/debian system install
    #  message ( STATUS "PETSC found in /usr/lib/petsc/" )
    #  set(PETSC_DIR /usr/lib/petsc/)

    #  find_package (PETSc 3.4 REQUIRED)
    #  petsc_get_version ()
    #  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev   (>= 3.4)") # This is not fully working yet. 

    #elseif ( IS_DIRECTORY "/usr/local/lib/python2.7/dist-packages/petsc") #Case ubuntu/debian system pip install
    #  message ( STATUS "PETSC found in /usr/local/lib/python2.7/dist-packages/petsc" )
    #  set(PETSC_DIR /usr/local/lib/python2.7/dist-packages/petsc)
    #  set(PETSC_INCLUDES  /usr/local/lib/python2.7/dist-packages/petsc/include /usr/include/openmpi)
    #  set(PETSC_LIBRARIES /usr/local/lib/python2.7/dist-packages/petsc/lib/libpetsc.so)

    #  set(PETSC_VERSION "3.8") #value for Ubuntu 16.04 

    else ()  # No petsc in system, do compile petsc
      message (STATUS "PETSC not found in the system, petsc 3.8.0 will be installed" )
      set(PETSC_SRC ${CDMATH_SOURCE_DIR}/pre_requis/PETSc)
      set(PETSC_DIR ${CMAKE_INSTALL_PREFIX}/share/petsc-3.8.0) # Location of the final install 
      set(PETSC_INCLUDES  ${PETSC_DIR}/include/  ${PETSC_DIR}/include/petsc/mpiuni)
      set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc.so)
      set(PETSC_VERSION "3.8.0")  

      include(ExternalProject)

      ExternalProject_Add (PETSc
      URL               ${PETSC_SRC}/petsc-lite-3.8.0.tar.gz   
      SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/petsc-3.8.0
      BUILD_IN_SOURCE   TRUE
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${PETSC_DIR} --with-debugging=0 --with-mpi=0 --download-fblaslapack=${PETSC_SRC}/fblaslapack-3.4.2.tar.gz
      BUILD_COMMAND     make all
      TEST_COMMAND      make test
      INSTALL_COMMAND   make install
      INSTALL_DIR       ${PETSC_DIR}
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
       )

    endif ()
  endif ()

  if (${PETSC_VERSION} VERSION_GREATER 3.5)
    add_definitions(-DPETSC_VERSION_GREATER_3_5)
  endif ()

endif ()

# Swig interface
if (CDMATH_WITH_PYTHON)
   find_package (PYTHON REQUIRED)
   find_package (SWIG REQUIRED)
   set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")
   set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python-devel, numpy, swig")
   add_subdirectory (${CDMATH_SWIG_DIR})
endif ()

add_subdirectory (${PRE_REQUIS_DIR})
add_subdirectory (${BASE_DIR})
add_subdirectory (${MESH_DIR})

if (CDMATH_WITH_PETSC)
    add_subdirectory (${LINEARSOLVER_DIR})
endif ()

# Documentation
if (CDMATH_WITH_DOCUMENTATION)
  find_package (Doxygen)
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, doxygen, graphviz, mscgen")
  set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, doxygen, graphviz, mscgen")
  if (NOT DOXYGEN_FOUND)
    message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
  endif (NOT DOXYGEN_FOUND)
  configure_file (Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
  							  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  if ( IS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html)
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/cdmath)
  endif()
endif ()

# Tests
if (CDMATH_WITH_TESTS)
  set (CMAKE_EXE_LINKER_FLAGS "-pg -fprofile-arcs -ftest-coverage -lgcov")
  enable_testing ()
  find_package (CPPUNIT REQUIRED)
  add_subdirectory (${TESTS_DIR})
endif ()

# Packaging
if (CDMATH_WITH_PACKAGE)
  include (CPackLists.txt)
endif ()

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)
  add_subdirectory (postprocessing)
  set (PV_LIB_DIR /usr/lib/paraview/:/usr/lib64/paraview/)
  set (PV_PYTHON_DIR /usr/lib/paraview/site-packages/paraview/:/usr/lib/paraview/site-packages/paraview/vtk:/usr/lib64/paraview/site-packages/paraview/:/usr/lib64/paraview/site-packages/vtk)

  add_custom_target (jupyter COMMAND jupyter-nbconvert --to html ${PROJECT_BINARY_DIR}/tests/validation/test_synthesis/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to pdf ${PROJECT_BINARY_DIR}/tests/validation/test_synthesis/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to latex ${PROJECT_BINARY_DIR}/tests/validation/test_synthesis/test_synthesis_nbook.ipynb)
endif ()

# Configuration file
configure_file(
    ${PROJECT_SOURCE_DIR}/env_CDMATH.sh
    ${CMAKE_INSTALL_PREFIX}/env_CDMATH.sh
    @ONLY
)
													    

