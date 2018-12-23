
SET(SCRIPT
  ./TransportEquation1DUpwindExplicit.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DTransportEquationUpwindExplicit)

if (CDMATH_WITH_PYTHON )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL1_NonSmooth ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )
    SET_PROPERTY(TEST validation1DTransportEquationUpwindExplicit_CFL1_NonSmooth PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST validation1DTransportEquationUpwindExplicit_CFL1_NonSmooth APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL1_Smooth ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    SET(SMOOTHNESS  0 )

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL0.5_NonSmooth ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )


    SET(SMOOTHNESS  1 )

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL0.5_Smooth ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL} ${SMOOTHNESS} )

endif (CDMATH_WITH_PYTHON )


