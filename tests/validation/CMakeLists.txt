
SET( VALIDATION_FOLDERS
2DPoissonEF 
2DPoissonVF_checkerboard 
2DPoissonVF_loc_ref 
2DPoissonVF_triangles 
2DPoissonVF_squares 
2DPoissonVF_deformedQuadrangles 
2DPoissonVF_hexagons 
2DPoissonVF_brickwall
2DPoissonVF_cross_triangles
2DPoissonVF_long_rectangles
2DPoissonVF_long_triangles
2DPoissonVF_flat_triangles
3DSpherePoissonEF 
3DPoissonEF 
3DPoissonVF_checkerboard 
3DPoissonVF_tetrahedra 
3DPoissonVF_cubes 
2DWaveSystemPStagSquares 
2DWaveSystemPStagTriangles 
2DWaveSystemPStagCrossTriangles 
2DWaveSystemPStagFlatCrossTriangles 
2DWaveSystemPStagCheckerboard
2DWaveSystemPStagBrickWall
2DWaveSystemStaggeredSquares 
2DWaveSystemCenteredSquares 
2DWaveSystemCenteredTriangles 
2DWaveSystemCenteredCrossTriangles 
2DWaveSystemCenteredFlatCrossTriangles 
2DWaveSystemCenteredCheckerboard
2DWaveSystemCenteredBrickWall
2DWaveSystemUpwindSquares 
2DWaveSystemUpwindTriangles 
2DWaveSystemUpwindCrossTriangles 
2DWaveSystemUpwindFlatCrossTriangles 
2DWaveSystemUpwindCheckerboard
2DWaveSystemUpwindBrickWall
2DWaveSystemUpwindHexagons
2DWaveSystemUpwindDeformedQuadrangles
3DWaveSystemUpwindCubes 
3DWaveSystemUpwindTetrahedra 
test_synthesis 
)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

  foreach(dir ${VALIDATION_FOLDERS} )
    ADD_SUBDIRECTORY(${dir})

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/validation/${dir}
      FILES_MATCHING PATTERN "*.py"
                     PATTERN "*.med"
                     PATTERN "*_0.vtu"
                     PATTERN "*.json"
    )

    string(FIND ${dir} "2DPoissonVF" is_2DPoissonVF)
    if(NOT ${is_2DPoissonVF} EQUAL -1)
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/2DPoissonVF
      FILES_MATCHING PATTERN "*_SQUARE.ipynb"
                     PATTERN "*_SQUARE.pdf"
    )
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/2DPoissonVF/${dir}
      FILES_MATCHING PATTERN "*.png"
                     PATTERN "*.ipynb"
                     PATTERN "*.py"
                     PATTERN "*.html"
                     PATTERN "*.pdf"
                     PATTERN "*.tex"
                     PATTERN "*.md"
                     PATTERN "*.xlsx"
                     PATTERN "*.csv"
                     PATTERN "*.json"
    )
    endif(NOT ${is_2DPoissonVF} EQUAL -1)
    string(FIND ${dir} "2DWaveSystem" is_2DWaveSystem)
    if(NOT ${is_2DWaveSystem} EQUAL -1)
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/2DWaveSystemVF_stationary
      FILES_MATCHING PATTERN "*_SQUARE.ipynb"
                     PATTERN "*_SQUARE.pdf"
    )
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/2DWaveSystemVF_stationary/${dir}
      FILES_MATCHING PATTERN "*.png"
                     PATTERN "*.ipynb"
                     PATTERN "*.py"
                     PATTERN "*.html"
                     PATTERN "*.pdf"
                     PATTERN "*.tex"
                     PATTERN "*.md"
                     PATTERN "*.xlsx"
                     PATTERN "*.csv"
                     PATTERN "*.json"
    )
    endif(NOT ${is_2DWaveSystem} EQUAL -1)
  endforeach(dir)

add_custom_target (
notebook COMMAND 
     jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis/test_synthesis_nbook.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis/test_synthesis_nbook.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis/test_synthesis_nbook.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis/test_synthesis_nbook.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/test_synthesis/test_synthesis_nbook.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_brickwall/Convergence_Poisson_FV5_SQUARE_brickwall.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_brickwall/Convergence_Poisson_FV5_SQUARE_brickwall.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_brickwall/Convergence_Poisson_FV5_SQUARE_brickwall.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_brickwall/Convergence_Poisson_FV5_SQUARE_brickwall.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_brickwall/Convergence_Poisson_FV5_SQUARE_brickwall.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_hexagons/Convergence_Poisson_FV5_SQUARE_hexagons.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_squares/Convergence_Poisson_FV5_SQUARE_squares.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_squares/Convergence_Poisson_FV5_SQUARE_squares.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_squares/Convergence_Poisson_FV5_SQUARE_squares.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_squares/Convergence_Poisson_FV5_SQUARE_squares.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_squares/Convergence_Poisson_FV5_SQUARE_squares.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_triangles/Convergence_Poisson_FV5_SQUARE_triangles.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_triangles/Convergence_Poisson_FV5_SQUARE_triangles.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_triangles/Convergence_Poisson_FV5_SQUARE_triangles.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_triangles/Convergence_Poisson_FV5_SQUARE_triangles.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_triangles/Convergence_Poisson_FV5_SQUARE_triangles.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_checkerboard/Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_checkerboard/Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_checkerboard/Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_checkerboard/Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_checkerboard/Convergence_Poisson_FV5_SQUARE_checkerboard.ipynb

&&   jupyter-nbconvert --to notebook --execute ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb
&&   jupyter-nbconvert --to html     ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb
&&   jupyter-nbconvert --to pdf      ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb
&&   jupyter-nbconvert --to latex    ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb
&&   jupyter-nbconvert --to markdown ${CMAKE_CURRENT_BINARY_DIR}/2DPoissonVF_loc_ref/Convergence_Poisson_FV5_SQUARE_loc_ref.ipynb

)


ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

