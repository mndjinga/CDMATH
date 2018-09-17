
SET(SCRIPT
  ./test_synthesis_nbook.ipynb
  ./test_synthesis_all.csv
  ./test_synthesis_short.csv
  ./test_synthesis_all.xlsx
  ./test_synthesis_short.xlsx
  )

file(COPY ${SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${SCRIPT} DESTINATION share/validation/)

