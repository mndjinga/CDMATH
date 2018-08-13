
SET(SCRIPT
  ./test_synthesis_nbook.ipynb
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/)

