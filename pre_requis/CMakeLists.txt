# Copyright (C) 2015-2016  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

# Directories
#
# Directories have to be given after prerequisites (to be able to use
# Python version string for example).
# ===========
SET(MEDCOUPLING_INSTALL_BINS bin/cdmath CACHE PATH "Install path: MEDCoupling binaries")
SET(MEDCOUPLING_INSTALL_LIBS lib CACHE PATH "Install path: MEDCoupling libs")
SET(MEDCOUPLING_INSTALL_HEADERS include CACHE PATH "Install path: MEDCoupling headers")
SET(MEDCOUPLING_INSTALL_SCRIPT_SCRIPTS ${MEDCOUPLING_INSTALL_BINS} CACHE PATH "Install path: MEDCoupling scripts")
SET(MEDCOUPLING_INSTALL_TESTS tests CACHE PATH "Install path: MEDCoupling tests")
# SET(MEDCOUPLING_INSTALL_SCRIPT_DATA ${MEDCOUPLING_INSTALL_BINS} CACHE PATH "Install path: MEDCoupling script data")
SET(MEDCOUPLING_INSTALL_SCRIPT_PYTHON ${MEDCOUPLING_INSTALL_BINS} CACHE PATH "Install path: MEDCoupling Python scripts")
SET(MEDCOUPLING_INSTALL_CMAKE_LOCAL cmake_files CACHE PATH "Install path: local MEDCoupling CMake files")

IF (CMAKE_CDMATH_SWIG)
  SET(MEDCOUPLING_ENABLE_PYTHON ON)
  SET(_pydir lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
  SET(MEDCOUPLING_INSTALL_PYTHON ${_pydir} CACHE PATH "Install path: MEDCoupling Python stuff")
  SET(MEDCOUPLING_INSTALL_PYTHON_SHARED ${MEDCOUPLING_INSTALL_PYTHON}/shared_modules CACHE PATH
    "Install path: MEDCoupling Python shared modules")
ENDIF (CMAKE_CDMATH_SWIG)

SET(MEDCOUPLING_INSTALL_RES share/resources CACHE PATH "Install path: MEDCoupling resources")
SET(MEDCOUPLING_INSTALL_DOC share/doc CACHE PATH "Install path: MEDCoupling documentation")

# Med specific:

SET(MEDCOUPLING_INSTALL_RES_DATA "${MEDCOUPLING_INSTALL_RES}/med" CACHE PATH "Install path: MEDCoupling specific data")
#SET(MEDCOUPLING_INSTALL_RES_SCRIPTS "${MEDCOUPLING_INSTALL_RES}/med" CACHE PATH "Install path: MEDCouplng specific scripts")

MARK_AS_ADVANCED(MEDCOUPLING_INSTALL_BINS MEDCOUPLING_INSTALL_LIBS MEDCOUPLING_INSTALL_IDLS MEDCOUPLING_INSTALL_HEADERS)
MARK_AS_ADVANCED(MEDCOUPLING_INSTALL_SCRIPT_SCRIPTS MEDCOUPLING_INSTALL_SCRIPT_DATA MEDCOUPLING_INSTALL_SCRIPT_PYTHON)
MARK_AS_ADVANCED(MEDCOUPLING_INSTALL_APPLISKEL_SCRIPTS  MEDCOUPLING_INSTALL_APPLISKEL_PYTHON MEDCOUPLING_INSTALL_CMAKE MEDCOUPLING_INSTALL_CMAKE_LOCAL MEDCOUPLING_INSTALL_RES)
MARK_AS_ADVANCED(MEDCOUPLING_INSTALL_PYTHON MEDCOUPLING_INSTALL_PYTHON_SHARED MEDCOUPLING_INSTALL_RES_DATA MEDCOUPLING_INSTALL_DOC)


# InterpKERNEL
ADD_SUBDIRECTORY(INTERP_KERNEL)

# MEDCoupling
ADD_SUBDIRECTORY(MEDCoupling)
IF(MEDCOUPLING_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(MEDCoupling_Swig)
ENDIF(MEDCOUPLING_ENABLE_PYTHON)

# MEDLoader
ADD_SUBDIRECTORY(MEDLoader)

# MEDFile
ADD_SUBDIRECTORY(MEDFile)

