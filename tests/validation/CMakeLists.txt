
SET( VALIDATION_FOLDERS
test_validation2DPoissonEF 
test_validation2DPoissonVF_checkerboard 
test_validation2DPoissonVF_loc_ref 
test_validation2DPoissonVF_triangles 
test_validation2DPoissonVF_squares 
test_validation2DPoissonVF_deformedQuadrangles 
test_validation2DPoissonVF_hexagons 
test_validation2DPoissonVF_brickwall
test_validation3DSpherePoissonEF 
test_validation3DPoissonEF 
test_validation3DPoissonVF_checkerboard 
test_validation3DPoissonVF_tetrahedra 
test_validation3DPoissonVF_cubes 
test_validation2DWaveSystemPStagSquares 
test_validation2DWaveSystemPStagTriangles 
test_validation2DWaveSystemUpwindSquares 
test_validation2DWaveSystemUpwindTriangles 
test_validation2DWaveSystemUpwindCheckerboard
test_validation2DWaveSystemUpwindBrickWall
test_validation2DWaveSystemUpwindHexagons
test_validation2DWaveSystemUpwindDeformedQuadrangles
test_validation3DWaveSystemUpwindCubes 
test_validation3DWaveSystemUpwindTetrahedra 
test_synthesis 
)

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)

  foreach(dir ${VALIDATION_FOLDERS} )

    ADD_SUBDIRECTORY(${dir})

    string(REPLACE "test_validation" "" short_name_dir ${dir})

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/validation/${short_name_dir}
      FILES_MATCHING PATTERN "*.py"
                     PATTERN "*.med"
                     PATTERN "*_0.vtu"
                     PATTERN "*.json"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      DESTINATION share/doc/convergence/${short_name_dir}
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

ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)

