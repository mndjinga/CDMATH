
SET(SCRIPT
  ./1DTransportEquationUpwindExplicit.py
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/test_validation1DTransportEquationUpwindExplicit)

if (CDMATH_WITH_PYTHON )

    SET(CFL  0.99  )#Courant Friedrichs Lewy number

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL})

    SET(CFL  0.5  )#Courant Friedrichs Lewy number

    ADD_TEST(validation1DTransportEquationUpwindExplicit_CFL0.5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_validation1DTransportEquationUpwindExplicit.py  ${CFL})

endif (CDMATH_WITH_PYTHON )


