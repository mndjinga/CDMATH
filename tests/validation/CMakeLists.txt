
SET( VALIDATION_FOLDERS
2DPoissonEF 
2DPoissonVF_checkerboard 
2DPoissonVF_loc_ref 
2DPoissonVF_triangles 
2DPoissonVF_squares 
2DPoissonVF_deformedQuadrangles 
2DPoissonVF_hexagons 
2DPoissonVF_brickwall
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
2DWaveSystemStaggeredSquares 
2DWaveSystemCenteredSquares 
2DWaveSystemCenteredTriangles 
2DWaveSystemUpwindSquares 
2DWaveSystemUpwindTriangles 
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

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/${dir}
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

