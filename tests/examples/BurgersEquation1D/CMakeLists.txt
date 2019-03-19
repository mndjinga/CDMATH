
if (CDMATH_WITH_PYTHON )

    ADD_TEST(ExampleBurgersEquation_1DFV_shock_formation ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/1DBurgersEquation_shock_formation.py )

endif (CDMATH_WITH_PYTHON )


