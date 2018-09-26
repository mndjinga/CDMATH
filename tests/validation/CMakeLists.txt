
IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)
    ADD_SUBDIRECTORY(test_validation2DPoissonEF)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_checkerboard)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_loc_ref)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_triangles)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_squares)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_deformedQuadrangles)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_hexagons)
    ADD_SUBDIRECTORY(test_validation2DPoissonVF_brickwall)
    ADD_SUBDIRECTORY(test_validation3DSpherePoissonEF)
    ADD_SUBDIRECTORY(test_validation3DPoissonEF)
    ADD_SUBDIRECTORY(test_validation3DPoissonVF_checkerboard)
    ADD_SUBDIRECTORY(test_validation3DPoissonVF_tetrahedra)
    ADD_SUBDIRECTORY(test_validation3DPoissonVF_cubes)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemPStagSquares)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemPStagTriangles)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindSquares)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindTriangles)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindCheckerboard)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindBrickWall)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindHexagons)
    ADD_SUBDIRECTORY(test_validation2DWaveSystemUpwindDeformedQuadrangles)
    ADD_SUBDIRECTORY(test_validation3DWaveSystemUpwindCubes)
    ADD_SUBDIRECTORY(test_validation3DWaveSystemUpwindTetrahedra)
    ADD_SUBDIRECTORY(test_synthesis)
ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC)


file(GLOB VALIDATION_FOLDERS_TO_INSTALL 
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
scripts)

install(DIRECTORY ${VALIDATION_FOLDERS_TO_INSTALL} DESTINATION share/validation)

install(DIRECTORY ${VALIDATION_FOLDERS_TO_INSTALL} DESTINATION share/doc/validation
FILES_MATCHING PATTERN "*.png"
               PATTERN "*.ipynb"
               PATTERN "*.py"
               PATTERN "*.html"
               PATTERN "*.pdf"
               PATTERN "*.tex"
               PATTERN "*.md"
               PATTERN "*.xlsx"
               PATTERN "*.csv"
)

