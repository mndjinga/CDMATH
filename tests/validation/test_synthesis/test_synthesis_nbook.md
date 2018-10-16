

```python
import glob
import json
import pandas as pd
```

Let us gather the name of all json files into a list


```python
description_files = glob.glob('../*/test_*.json')

print("The following test description files : ")
print(description_files)
print("will be imported")
```

    The following test description files : 
    ['../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles16Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles256Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles64Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles1024Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall2500Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall900Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall225Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall25Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_25Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_11500Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_1020Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_2850Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_255Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles6422Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles934Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles40Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles224Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_21542Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_572Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_9020Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_2244Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_5272Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares16Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares4096Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares64Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares1024Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares256Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra55566Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra7986Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra750Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra105456Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_2560Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_160Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_640Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_10240Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_163840Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_40Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_40960Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_40Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_224Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_6422Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_25872Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_934Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_9261Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_68921Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_1331Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_256Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_1024Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_4096Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_16Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_64Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares64Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares16Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares256Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares1024Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares4096Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes216Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes9261Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes1331Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_6422Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_224Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_25872Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_40Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_934Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles6422Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles224Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles40Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles934Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_16834Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_2081Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_7629Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_28561Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_63249Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_270Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_4077Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons255Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons25Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons1020Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons2850Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_7629Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_16834Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_63249Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_4077Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_28561Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_270Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_2081Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_36Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_18432Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_2304Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_288Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard40Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard2560Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard160Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard640Cells.json']
    will be imported


Each json file content will be imported into a python dict, all these dict will be gathered into a list called `all_descriptions`


```python
# Let's import all json files into a list of dictionaries
all_descriptions = []

for file_name in description_files:
    with open(file_name, 'r') as fd:
        all_descriptions.append(json.load(fd))

print("json files have been imported")
```

    json files have been imported


In order to print a list (or a sublist), we need to import the pprint python package


```python
import pprint as pp
#pp.pprint(all_descriptions)
```


```python
# Let's create a pandas dataframe out of our dict list
df = pd.DataFrame(all_descriptions)
print("The pandas dataframe has been created")

list_of_all_columns = df.columns
print("Printing the columns of the dataframe : these are the parameters of the database")
pp.pprint(list_of_all_columns)
```

    The pandas dataframe has been created
    Printing the columns of the dataframe : these are the parameters of the database
    Index([u'Absolute_error', u'Boundary_conditions',
           u'Computational_time_taken_by_run', u'Geometry', u'Global_comment',
           u'Global_name', u'Initial_data', u'Linear_solver_algorithm',
           u'Linear_solver_maximum_iterations', u'Linear_solver_precision',
           u'Linear_solver_preconditioner', u'Linear_solver_with_scaling',
           u'Linear_system_max_actual_condition number',
           u'Linear_system_max_actual_error',
           u'Linear_system_max_actual_iterations_number', u'Mesh_cell_type',
           u'Mesh_dimension', u'Mesh_is_unstructured',
           u'Mesh_max_number_of_neighbours', u'Mesh_number_of_elements',
           u'Mesh_type', u'Numerical_method_name',
           u'Numerical_method_space_discretization',
           u'Numerical_method_time_discretization', u'Numerical_parameter_cfl',
           u'Numerical_parameter_space_step', u'Numerical_parameter_time_step',
           u'PDE_is_stationary', u'PDE_model',
           u'PDE_search_for_stationary_solution',
           u'Part_of_mesh_convergence_analysis', u'Relative_error',
           u'Simulation_final_number_of_time_steps_after_run',
           u'Simulation_final_time_after_run', u'Simulation_output_frequency',
           u'Simulation_parameter_maximum_time',
           u'Simulation_parameter_maximum_time_step', u'Space_dimension',
           u'Test_color'],
          dtype='object')



```python
#print("Printing the dataframe")
#pp.pprint(df)

# print values of columns: the name of the column can be used as attributes
#print("Values of the Name column")
#print(df.Global_name)
```


```python
# a new dataframe with a few columns only
#Tous les résultats avec cfl >1
column_list = ['Geometry', 'Numerical_parameter_cfl']
sub_df1 = df[column_list]
#print("sub_df1")
#pp.pprint(sub_df1)
```


```python
# a new dataframe according to the CFL value
sub_df2 = df[df.Numerical_parameter_cfl > 0.1]
#print("sub_df2")
#pp.pprint(sub_df2)
```


```python
# sorting a dataframe
df.sort_values(by=['Global_name', 'Numerical_parameter_cfl'], ascending=True)

sub_df3 = df[df['Boundary_conditions'].isin(['Dirichlet'])]
#print("sub_df3")
#pp.pprint(sub_df3)
```

# Displaying validation test tables with qgrid

Let's play with `qgrid now`.
First extract the most interesting columns and visualise them in a widget.


```python
import qgrid

# here's a cool dictionnary of options for displaying data
gopt={
    'fullWidthRows': True,
    'syncColumnCellResize': True,
    'forceFitColumns': True,
    'defaultColumnWidth': 150,
    'rowHeight': 28,
    'enableColumnReorder': True,
    'enableTextSelectionOnCells': True,
    'editable': False,
    'autoEdit': False,
    'explicitInitialization': True,
    'maxVisibleRows': 40,
    'minVisibleRows': 8,
    'sortable': True,
    'filterable': True,
    'highlightSelectedCell': False,
    'highlightSelectedRow': True
}

# Extract the most interesting column from df into a second dataframe df2
df2=df[['PDE_model','Numerical_method_name','Mesh_dimension','Mesh_type','Mesh_cell_type','Test_color','Mesh_number_of_elements','Computational_time_taken_by_run']]

# Let's create a jupyter table widget from the dataframe df2
qgrid_widget=qgrid.show_grid(df2, grid_options=gopt, show_toolbar=False)

# let's output this widget
qgrid_widget
```


    QgridWidget(grid_options={'defaultColumnWidth': 150, 'highlightSelectedRow': True, 'enableTextSelectionOnCells…


# Exporting validation test table to CSV and Excel format

`pandas` can be used to export to csv and excel, this is useful!

Let us first export the large database df.


```python
df.to_csv('test_synthesis_all.csv')#Saving using csv format
output_file_name='test_synthesis_all.xlsx'
writer = pd.ExcelWriter(output_file_name)
df.to_excel(writer,'Sheet1')
writer.save()
print("Done writing file "+output_file_name)
```

    Done writing file test_synthesis_all.xlsx


Let us now export the short database df2.


```python
df2.to_csv('test_synthesis_short.csv')#Saving using csv format
output_file_name='test_synthesis_short.xlsx'
writer = pd.ExcelWriter(output_file_name)
df2.to_excel(writer,'Sheet1')
writer.save()
print("Done writing file "+output_file_name)
```

    Done writing file test_synthesis_short.xlsx



```python
ls
```

    [0m[38;5;33mCMakeFiles[0m/          Makefile                 test_synthesis_nbook.ipynb
    cmake_install.cmake  test_synthesis_all.csv   test_synthesis_short.csv
    CTestTestfile.cmake  test_synthesis_all.xlsx  test_synthesis_short.xlsx


# Convergence study table


```python
convergence_files = glob.glob('../*/Convergence_*.json')

print("The following convergence description files : ")
print(description_files)
print("will be imported")
```

    The following convergence description files : 
    ['../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles16Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles256Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles64Cells.json', '../test_validation2DWaveSystemUpwindDeformedQuadrangles/test_WaveSystem2DUpwind_SquareWithDeformedQuadrangles1024Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall2500Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall900Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall225Cells.json', '../test_validation2DWaveSystemUpwindBrickWall/test_WaveSystem2DUpwind_squareWithBrickWall25Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_25Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_11500Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_1020Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_2850Cells.json', '../test_validation2DPoissonVF_hexagons/test_Poisson2D_VF_255Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles6422Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles934Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles40Cells.json', '../test_validation2DWaveSystemUpwindTriangles/test_WaveSystem2DUpwind_squareWithTriangles224Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_21542Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_572Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_9020Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_2244Cells.json', '../test_validation3DSpherePoissonEF/test_Poisson2D_EF_5272Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares16Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares4096Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares64Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares1024Cells.json', '../test_validation2DWaveSystemUpwindSquares/test_WaveSystem2DUpwind_meshSquareWithSquares256Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra55566Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra7986Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra750Cells.json', '../test_validation3DWaveSystemUpwindTetrahedra/test_WaveSystem3DUpwind_CubeWithTetrahedra105456Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_2560Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_160Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_640Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_10240Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_163840Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_40Cells.json', '../test_validation2DPoissonVF_loc_ref/test_Poisson2D_VF_40960Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_40Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_224Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_6422Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_25872Cells.json', '../test_validation2DPoissonEF/test_Poisson2D_EF_934Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_9261Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_68921Cells.json', '../test_validation3DPoissonVF_cubes/test_Poisson3D_VF_1331Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_256Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_1024Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_4096Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_16Cells.json', '../test_validation2DPoissonVF_squares/test_Poisson2D_VF_64Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares64Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares16Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares256Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares1024Cells.json', '../test_validation2DWaveSystemPStagSquares/test_WaveSystem2DPStag_SquareWithSquares4096Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes216Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes9261Cells.json', '../test_validation3DWaveSystemUpwindCubes/test_WaveSystem3DUpwind_CubeWithCubes1331Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_6422Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_224Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_25872Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_40Cells.json', '../test_validation2DPoissonVF_triangles/test_Poisson2D_VF_934Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles6422Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles224Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles40Cells.json', '../test_validation2DWaveSystemPStagTriangles/test_WaveSystem2DPStag_SquareWithTriangles934Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_16834Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_2081Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_7629Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_28561Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_63249Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_270Cells.json', '../test_validation3DPoissonEF/test_Poisson3D_EF_4077Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons255Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons25Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons1020Cells.json', '../test_validation2DWaveSystemUpwindHexagons/test_WaveSystem2DUpwind_squareWithHexagons2850Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_7629Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_16834Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_63249Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_4077Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_28561Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_270Cells.json', '../test_validation3DPoissonVF_tetrahedra/test_Poisson3D_VF_2081Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_36Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_18432Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_2304Cells.json', '../test_validation3DPoissonVF_checkerboard/test_Poisson3D_VF_288Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard40Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard2560Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard160Cells.json', '../test_validation2DWaveSystemUpwindCheckerboard/test_WaveSystem2DUpwind_squareWithCheckerboard640Cells.json']
    will be imported



```python
# Let's import all json files into a list of dictionaries
convergence_descriptions = []

for file_name in convergence_files:
    with open(file_name, 'r') as fd:
        convergence_descriptions.append(json.load(fd))

print("convergence json files have been imported")
```

    convergence json files have been imported



```python
# Let's create a pandas dataframe out of the convergence dictionary
df_convergence = pd.DataFrame(convergence_descriptions)
print("The convergence pandas dataframe has been created")

list_of_all_columns_convergence = df_convergence.columns
print("Printing the columns of the dataframe : these are the parameters of the database")
pp.pprint(list_of_all_columns_convergence)
```

    The convergence pandas dataframe has been created
    Printing the columns of the dataframe : these are the parameters of the database
    Index([u'Boundary_conditions', u'Color', u'Computational_time',
           u'Condition_numbers', u'Errors', u'Final_time', u'Final_time_step',
           u'Geometry', u'Global_comment', u'Global_name', u'Initial_data',
           u'Max_vel_norm', u'Mesh_cell_type', u'Mesh_description',
           u'Mesh_dimension', u'Mesh_is_unstructured', u'Mesh_names', u'Mesh_path',
           u'Mesh_sizes', u'Mesh_type', u'Numerical_error_pressure',
           u'Numerical_error_velocity', u'Numerical_method_name',
           u'Numerical_method_space_discretization',
           u'Numerical_method_time_discretization', u'Numerical_parameter_cfl',
           u'PDE_is_stationary', u'PDE_model',
           u'PDE_search_for_stationary_solution',
           u'Part_of_mesh_convergence_analysis', u'Scaling_preconditioner',
           u'Scheme_order', u'Scheme_order_press', u'Scheme_order_vel',
           u'Space_dimension', u'Test_color'],
          dtype='object')



```python
# Extract the most interesting column from df_convergence into a second dataframe df2_convergence
df2_convergence=df_convergence[['PDE_model','Numerical_method_name','Mesh_dimension','Mesh_type','Scheme_order','Mesh_cell_type','Test_color','Computational_time']]

# Let's create a jupyter table widget from the convergenc dataframe
qgrid_widget_convergence=qgrid.show_grid(df2_convergence, grid_options=gopt, show_toolbar=False)

# let's output this widget
qgrid_widget_convergence

```


    QgridWidget(grid_options={'defaultColumnWidth': 150, 'highlightSelectedRow': True, 'enableTextSelectionOnCells…



```python
#Now export convergence study table
df_convergence.to_csv('Convergence_table_all.csv')#Saving using csv format
output_file_name='Convergence_table_all.xlsx'
writer = pd.ExcelWriter(output_file_name)
df_convergence.to_excel(writer,'Sheet1')
writer.save()
print("Done writing file "+output_file_name)

df2_convergence.to_csv('Convergence_table_short.csv')#Saving using csv format
output_file_name='Convergence_table_short.xlsx'
writer = pd.ExcelWriter(output_file_name)
df2_convergence.to_excel(writer,'Sheet1')
writer.save()
print("Done writing file "+output_file_name)
```

    Done writing file Convergence_table_all.xlsx
    Done writing file Convergence_table_short.xlsx

