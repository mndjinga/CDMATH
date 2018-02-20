cmake_minimum_required (VERSION 2.8)
PROJECT(postprocessing)

INSTALL(FILES analyticalFunction_discretizer.py VTK_routines.py MED_routines.py DESTINATION bin/cdmath/postprocessing)

