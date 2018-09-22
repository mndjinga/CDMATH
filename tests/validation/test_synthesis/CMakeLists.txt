
SET(NOTEBOOK
  ./test_synthesis_nbook.ipynb
  )
SET(TABLES
  ./test_synthesis_all.csv
  ./test_synthesis_short.csv
  ./test_synthesis_all.xlsx
  ./test_synthesis_short.xlsx
  )

file(COPY ${NOTEBOOK} ${TABLES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${NOTEBOOK} ${TABLES} DESTINATION share/validation/)

add_custom_target (test_synthesis COMMAND jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}s/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis_nbook.ipynb
				&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis_nbook.ipynb)


