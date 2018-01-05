
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/analyticalFunction_discretizer.py DESTINATION bin/postprocessing)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/VTU_routines.py DESTINATION bin/postprocessing)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/MED_routines.py DESTINATION bin/postprocessing)
INSTALL(FILES ${MED_LOADER_DIR}/swig/VTKReader.py DESTINATION bin/postprocessing)
