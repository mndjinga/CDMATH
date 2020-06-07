cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

set (TESTS_DIR ${CDMATH_SOURCE_DIR}/tests)
list (APPEND CMAKE_MODULE_PATH "${CDMATH_SOURCE_DIR}/cmake_files")

# PETSc and HDF5
if (CDMATH_WITH_PETSC)
    message ( STATUS "Checking variable PETSC_DIR : " $ENV{PETSC_DIR} )

  if ( DEFINED ENV{PETSC_DIR} )

    find_package (PETSc 3.4 REQUIRED)
    petsc_get_version ()
    string(REPLACE ";" ":"  PETSC_INCLUDES_PATH "${PETSC_INCLUDES}")# use colon instead of semicolon in environment file env_CDMATH.sh
    message ( STATUS "PETSC found. Version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" )

    #define slepc variables
    set(SLEPC_DIR ${PETSC_DIR}/${PETSC_ARCH})
    message ( STATUS "Checking variable SLEPC_DIR" )
    if ( IS_DIRECTORY ${SLEPC_DIR}/include AND EXISTS ${SLEPC_DIR}/lib/libslepc.so)
      set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
      set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
      message( STATUS "SLEPC found at ${SLEPC_DIR}" )
    else()
      message( WARNING "SLEPC not found at ${SLEPC_DIR}" )
    endif()

    #HDF5 to be found in petsc external packages
    set(HDF5_ROOT $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})# define hint for hdf5 installation in petsc external packages

  else ( DEFINED ENV{PETSC_DIR} )
    message ( STATUS "PETSC_DIR not set, searching petsc in the system" )

    if ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" AND DEFINED $ENV{SLEPC_DIR} )  #Case fedora/redhat system install
      message ( STATUS "PETSC includes found in /usr/include/petsc/" )
      message ( STATUS "PETSC library found in /usr/lib64" )
      set(PETSC_DIR /usr/)
      set(PETSC_INCLUDES  /usr/include/petsc /usr/include/petsc/petsc/mpiuni)
      set(PETSC_LIBRARIES /usr/lib64/libpetsc.so)

      #set(PETSC_VERSION "3.8") #3.8 for fedora 26 , 3.9 for fedora 29 , 3.10 for fedora 30, , 3.12 for fedora 32
      petsc_get_version ()
      set (CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, petsc-devel (>= 3.4)") # This is not fully working yet. 

      message ( STATUS "Checking variable SLEPC_DIR : " $ENV{SLEPC_DIR} )
      if ( IS_DIRECTORY ${SLEPC_DIR}/include AND EXISTS ${SLEPC_DIR}/lib/libslepc.so)
        set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
        set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
        message( STATUS "SLEPC found at $ENV{SLEPC_DIR}" )
      else()
        message( WARNING "SLEPC not found" )
      endif()

      #HDF5 to be found in the system
      set (CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-dev")
      set (CPACK_RPM_PACKAGE_REQUIRES "hdf5-devel")

    #elseif ( IS_DIRECTORY "/usr/lib/petsc/" AND EXISTS "/usr/lib/petsc/lib/libpetsc_real.so") #Case ubuntu/debian system install
    #  message ( STATUS "PETSC found in /usr/lib/petsc/ abnd exists /usr/lib/petsc/lib/libpetsc_real.so")
    #  set(PETSC_DIR /usr/lib/petsc/)
    #  set(PETSC_INCLUDES  ${PETSC_DIR}/include ${PETSC_DIR}/include/petsc/mpiuni)
    #  set(PETSC_LIBRARIES ${PETSC_DIR}/libpetsc_real.so)

    #  find_package (PETSc 3.4 REQUIRED)
    #  petsc_get_version ()
    #  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev   (>= 3.4)") # This is not fully working yet. 

    #  set(SLEPC_DIR /usr/lib/slepc/)
    #  set(SLEPC_INCLUDES  ${SLEPC_DIR}/include )
    #  set(SLEPC_LIBRARIES ${SLEPC_DIR}/libslepc_real.so)

    #elseif ( IS_DIRECTORY "/usr/local/lib/python2.7/dist-packages/petsc") #Case ubuntu/debian system pip install
    #  message ( STATUS "PETSC found in /usr/local/lib/python2.7/dist-packages/petsc" )
    #  set(PETSC_DIR /usr/local/lib/python2.7/dist-packages/petsc)
    #  set(PETSC_INCLUDES  ${PETSC_DIR}/include /usr/include/openmpi)
    #  set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc.so)
    #  petsc_get_version ()
    #  set(PETSC_VERSION "3.4") #value for Ubuntu 14.04 
    #  set(PETSC_VERSION "3.6") #value for Ubuntu 16.04 
    #  set(PETSC_VERSION "3.7") #value for Ubuntu 18.04 

    else ()  # No petsc in system, do compile petsc along with slepc as an external package
      message (STATUS "PETSC not found in the system, petsc 3.13.1 will be installed" )
      set(PETSC_SRC ${CDMATH_SOURCE_DIR}/pre_requis/PETSc)
      set(PETSC_DIR ${CMAKE_INSTALL_PREFIX}/share/petsc-3.13.1) # Location of the final install 
      set(PETSC_INCLUDES  ${PETSC_DIR}/include/  ${PETSC_DIR}/include/petsc/mpiuni)
      set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc.so)
      set(PETSC_VERSION "3.13.1")  

      include(ExternalProject)

      ExternalProject_Add (PETSc
      URL               ${PETSC_SRC}/petsc-lite-3.13.1.tar.gz   
      SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/petsc-3.13.1
      BUILD_IN_SOURCE   TRUE
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${PETSC_DIR} --with-debugging=0 --with-mpi=0 --download-f2cblaslapack=${PETSC_SRC}/f2cblaslapack-3.4.2.q4.tar.gz --with-fc=0 --download-slepc=${PETSC_SRC}/slepc-3.13.1.tar.gz --download-hdf5=${PETSC_SRC}/hdf5-1.10.3.tar.gz
      BUILD_COMMAND     make all
      TEST_COMMAND      make check
      INSTALL_COMMAND   make install
      INSTALL_DIR       ${PETSC_DIR}
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
       )

      #define slepc variables
      set(SLEPC_DIR ${PETSC_DIR}/${PETSC_ARCH})
      set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
      set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
      string(REPLACE ";" ":"  SLEPC_INCLUDES_PATH "${SLEPC_INCLUDE_DIR}")# use colon instead of semicolon in environment file env_CDMATH.sh

      #HDF5 to be found in petsc external packages
      set(HDF5_ROOT $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})# define hint for hdf5 installation in petsc external packages

    endif ()
  endif ( DEFINED ENV{PETSC_DIR} )

  if( ${PETSC_VERSION} VERSION_GREATER 3.5 )
    add_definitions(-DPETSC_VERSION_GREATER_3_5)
  endif( ${PETSC_VERSION} VERSION_GREATER 3.5 )

endif( CDMATH_WITH_PETSC )

#HDF5
find_package (HDF5 REQUIRED)

#MED
message (STATUS "MED installation" )
message (STATUS "HDF5 home is ${HDF5_ROOT}")
string(FIND "${HDF5_LIBRARIES}" "libhdf5.so" pos)
string(SUBSTRING "${HDF5_LIBRARIES}" 0 ${pos} HDF5_LIBRARY_DIR)
message (STATUS "HDF5 library dir is ${HDF5_LIBRARY_DIR}")
message (STATUS "HDF5 include dir is ${HDF5_INCLUDE_DIRS}")

set(MED_SRC ${CDMATH_SOURCE_DIR}/pre_requis/)
set(MACHINE PCLINUX)
set(MED3_DEFINITIONS "-D${MACHINE} -DMEDFILE_INSTALL_DOC=OFF")#Comment ne pas compiler la doc ?
set(MED_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/med-4.0.0_install) # Location of the final install

include(ExternalProject)

ExternalProject_Add (MED
      URL               ${MED_SRC}/med-4.0.0.tar.gz  #location of med tarball
      SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/med-4.0.0
      BUILD_IN_SOURCE   TRUE
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${MED_INSTALL_DIR}  --with-hdf5=${HDF5_ROOT} --with-hdf5-include=${HDF5_INCLUDE_DIRS} --with-hdf5-lib=${HDF5_LIBRARY_DIR} --with-hdf5-bin=${HDF5_ROOT}/bin  --with-swig=yes 
      BUILD_COMMAND     make
      INSTALL_COMMAND   make install
      INSTALL_DIR       ${MED_INSTALL_DIR}
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
)

set (MED3_DIR ${MED_INSTALL_DIR} )
set (MEDFILE_INCLUDE_DIRS ${MED3_DIR}/include)# A mettre à jour
set (MED3_INCLUDE_DIRS ${MED_INSTALL_DIR}/include)
set (MEDFILE_C_LIBRARIES med medC)

#MEDCoupling
message (STATUS "MEDCoupling installation" )

set(MEDCOUPLING_SRC ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/medCoupling-9.4.0/) # location of sources after unzip
set(MEDCOUPLING_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/medCoupling-9.4.0_install) # Location of the final install

include(ExternalProject)

ExternalProject_Add (MEDCoupling
      URL               ${CDMATH_SOURCE_DIR}/pre_requis/medCoupling-9.4.0.tar.gz  #location of medcoupling tarball
      SOURCE_DIR        ${MEDCOUPLING_SRC}
      BUILD_IN_SOURCE   TRUE
      CMAKE_COMMAND cmake <SOURCE_DIR>/medCoupling-9.4.0 -DCMAKE_INSTALL_PREFIX=${MEDCOUPLING_INSTALL_DIR} -DCONFIGURATION_ROOT_DIR=<SOURCE_DIR>/configuration-9.4.0/ -DMEDFILE_ROOT_DIR=${MED_INSTALL_DIR} -DMEDCOUPLING_ENABLE_PARTITIONER=OFF -DMEDCOUPLING_PARTITIONER_METIS=OFF -DMEDCOUPLING_PARTITIONER_SCOTCH=OFF -DMEDCOUPLING_ENABLE_RENUMBER=OFF -DMEDCOUPLING_BUILD_DOC=OFF
      BUILD_COMMAND     make
      INSTALL_COMMAND   make install
      INSTALL_DIR       ${MEDCOUPLING_INSTALL_DIR}
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
)

set (INTERP_KERNEL_DIR     ${MEDCOUPLING_SRC}/medCoupling-9.4.0/src/INTERP_KERNEL)
set (MED_COUPLING_DIR      ${MEDCOUPLING_SRC}/medCoupling-9.4.0/src/MEDCoupling)
set (MED_COUPLING_SWIG_DIR ${MEDCOUPLING_SRC}/medCoupling-9.4.0/src/MEDCoupling_Swig)
set (MED_LOADER_DIR        ${MEDCOUPLING_SRC}/medCoupling-9.4.0/src/MEDLoader)
set (MED_COUPLING_INCLUDE_DIR      ${MEDCOUPLING_INSTALL_DIR}/include)

# Paraview
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)
  set (PV_LIB_DIR /usr/lib/python2.7/dist-packages/:/usr/lib64/paraview/)
  set (PV_PYTHON_DIR /usr/lib/python2.7/dist-packages:/usr/lib64/paraview/site-packages/:/usr/lib64/paraview/site-packages/paraview/:/usr/lib64/paraview/site-packages/vtk/:/usr/lib64/paraview/python2.7/site-packages/:/usr/lib64/paraview/python2.7/site-packages/paraview/:/usr/lib64/paraview/python2.7/site-packages/vtkmodules:/usr/lib64/paraview/python3.7/site-packages/:/usr/lib64/paraview/python3.7/site-packages/paraview/:/usr/lib64/paraview/python3.7/site-packages/vtkmodules)
endif ()

# Swig interface
if (CDMATH_WITH_PYTHON)
   find_package (PYTHON REQUIRED)
   find_package (SWIG 3.0 REQUIRED)
   set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")
   set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python-devel, numpy, swig")
   add_subdirectory (${CDMATH_SWIG_DIR})
endif ()

find_package (XDR REQUIRED)
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
endif ()
if (CDMATH_WITH_TESTS AND CDMATH_WITH_PYTHON)
  enable_testing ()
  add_subdirectory (${TESTS_DIR})#contains c++ and python tests
endif()

# Packaging
if (CDMATH_WITH_PACKAGE)
  include (CPackLists.txt)
endif ()

# Postprocessing
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)
  add_subdirectory (postprocessing)
endif ()

# Configuration file
configure_file(
    ${PROJECT_SOURCE_DIR}/env_CDMATH.sh
    ${CMAKE_INSTALL_PREFIX}/env_CDMATH.sh
    @ONLY
)
													    

