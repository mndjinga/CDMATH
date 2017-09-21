// Copyright (C) 2015-2016  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#if !defined(__MEDCOUPLING_VERSION_H__)
#define __MEDCOUPLING_VERSION_H__

/*!
  Specify version of MEDCoupling, as follows

  MEDCOUPLING_VERSION_MAJOR       : (integer) number identifying major version
  MEDCOUPLING_VERSION_MINOR       : (integer) number identifying minor version
  MEDCOUPLING_VERSION_MAINTENANCE : (integer) number identifying maintenance version
  MEDCOUPLING_VERSION_STR         : (string)  complete version number "major.minor.maintenance"
  MEDCOUPLING_VERSION             : (hex)     complete version number (major << 16) + (minor << 8) + maintenance
*/

#define MEDCOUPLING_VERSION_MAJOR       8
#define MEDCOUPLING_VERSION_MINOR       3
#define MEDCOUPLING_VERSION_MAINTENANCE 0
#define MEDCOUPLING_VERSION_STR         "8.3.0"
#define MEDCOUPLING_VERSION             0x080300

#endif // __MEDCOUPLING_VERSION_H__
